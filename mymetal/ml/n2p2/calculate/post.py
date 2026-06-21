"""Post-processing of n2p2 ``nnp-train`` output (data layer).

Parses ``learning-curve.out`` (physical-unit RMSE columns), reads the
``trainpoints.NNNNNN.out`` / ``trainforces.NNNNNN.out`` of a chosen epoch
(training units if ``normalize_data_set`` was used) and converts them back to
physical units with the ``mean_energy`` / ``conv_energy`` / ``conv_length``
constants that nnp-train writes into ``input.nn``, and aggregates per-tag
errors into a DataFrame.

Tags are NOT parsed here: the structure ``tag=`` labels come from
``mymetal.ml.n2p2.dataset.nnpdata`` (single source of truth for input.data).
Plotting lives in ``mymetal.universal.plot.n2p2``.

Unit conversion (n2p2 normalization):
    E_pu/atom = E_iu/atom / conv_energy + mean_energy
    F_pu      = F_iu * conv_length / conv_energy
RMSE-type quantities only need the conv factors (mean_energy cancels in
differences); absolute energies need the full formula.
"""

import numpy as np
import pandas as pd


def read_learning_curve(file_path: str = None) -> np.ndarray:
    """Parse an n2p2 ``learning-curve.out`` file.

    Args:
        file_path: Path to learning-curve.out.

    Returns:
        Array of shape ``(n_epochs, n_cols)`` with all numeric columns.
        For n2p2 2.3.0 (0-based indices): col 0 = epoch,
        col 1 = RMSEpa_Etrain_pu (energy RMSE per atom, physical units),
        col 9 = RMSE_Ftrain_pu (force RMSE, physical units).
        Test columns are NaN when test_fraction is 0 ('-NAN' in the file).
    """
    rows = []
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            fields = line.split()
            if not fields or fields[0].startswith('#'):
                continue
            # float('-NAN') -> nan，无测试集时的占位值
            rows.append([float(s) for s in fields])
    return np.array(rows)


def read_normalization(file_path: str = None) -> tuple:
    """Read data-set normalization constants from an ``input.nn`` file.

    nnp-train (keyword ``normalize_data_set``) writes ``mean_energy``,
    ``conv_energy`` and ``conv_length`` back into input.nn. For a
    non-normalized training these keywords are absent and the returned
    defaults make all conversions identities.

    Args:
        file_path: Path to input.nn (the copy inside the training run dir).

    Returns:
        Tuple ``(mean_energy, conv_energy, conv_length)``,
        defaults ``(0.0, 1.0, 1.0)``.
    """
    mean_energy, conv_energy, conv_length = 0.0, 1.0, 1.0
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            fields = line.split('#', 1)[0].split()
            if len(fields) < 2:
                continue
            if fields[0] == 'mean_energy':
                mean_energy = float(fields[1])
            elif fields[0] == 'conv_energy':
                conv_energy = float(fields[1])
            elif fields[0] == 'conv_length':
                conv_length = float(fields[1])
    return mean_energy, conv_energy, conv_length


def read_trainpoints(file_path: str = None) -> np.ndarray:
    """Parse ``trainpoints.NNNNNN.out`` (energy comparison).

    Args:
        file_path: Path to the trainpoints file.

    Returns:
        Array ``(n_structures, 3)``: structure index (= 0-based position in
        input.data), Eref/atom, Ennp/atom (training units if normalized).
    """
    return read_learning_curve(file_path)


def read_trainforces(file_path: str = None) -> np.ndarray:
    """Parse ``trainforces.NNNNNN.out`` (force comparison).

    Args:
        file_path: Path to the trainforces file.

    Returns:
        Array ``(n_components, 4)``: structure index, atom index (x, y, z on
        consecutive rows), Fref, Fnnp (training units if normalized).
    """
    return read_learning_curve(file_path)


def rmse_me(ref: np.ndarray = None, pred: np.ndarray = None) -> tuple:
    """Root-mean-square error and mean error (pred - ref).

    Args:
        ref: Reference values.
        pred: Predicted values.

    Returns:
        Tuple ``(rmse, me)``.
    """
    diff = np.asarray(pred) - np.asarray(ref)
    return float(np.sqrt((diff ** 2).mean())), float(diff.mean())


def build_rmse_by_tag_df(e_ref: np.ndarray = None, e_nnp: np.ndarray = None, tag_e: np.ndarray = None,
                         f_ref: np.ndarray = None, f_nnp: np.ndarray = None, tag_f: np.ndarray = None) -> pd.DataFrame:
    """Aggregate energy/force RMSE and ME per tag into a DataFrame.

    Energies are per structure (per atom), forces per component. All errors
    are returned in meV/atom and meV/Angstrom. A final ``TOTAL`` row spans the
    whole data set. ME = mean(NNP - DFT).

    Args:
        e_ref: DFT energies per atom (eV/atom).
        e_nnp: NNP energies per atom (eV/atom).
        tag_e: Tag label of each energy point.
        f_ref: DFT force components (eV/Angstrom).
        f_nnp: NNP force components (eV/Angstrom).
        tag_f: Tag label of each force component.

    Returns:
        DataFrame with columns ``tag, n_struct, E_RMSE_meV/at, E_ME_meV/at,
        n_fcomp, F_RMSE_meV/A, F_ME_meV/A`` (one row per tag + TOTAL).
    """
    e_ref, e_nnp, tag_e = np.asarray(e_ref), np.asarray(e_nnp), np.asarray(tag_e)
    f_ref, f_nnp, tag_f = np.asarray(f_ref), np.asarray(f_nnp), np.asarray(tag_f)

    def _row(label, me_mask, mf_mask):
        e_rmse, e_me = rmse_me(e_ref[me_mask], e_nnp[me_mask])
        f_rmse, f_me = rmse_me(f_ref[mf_mask], f_nnp[mf_mask]) if mf_mask.any() else (np.nan, np.nan)
        return {'tag': label, 'n_struct': int(me_mask.sum()),
                'E_RMSE_meV/at': e_rmse * 1e3, 'E_ME_meV/at': e_me * 1e3,
                'n_fcomp': int(mf_mask.sum()),
                'F_RMSE_meV/A': f_rmse * 1e3, 'F_ME_meV/A': f_me * 1e3}

    rows = [_row(t, tag_e == t, tag_f == t) for t in sorted(set(tag_e))]
    rows.append(_row('TOTAL', np.ones(len(e_ref), bool), np.ones(len(f_ref), bool)))
    return pd.DataFrame(rows)
