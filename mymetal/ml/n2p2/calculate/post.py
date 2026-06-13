"""Post-processing of n2p2 ``nnp-train`` output.

Workflow: parse ``learning-curve.out`` (physical-unit RMSE columns), pick an
epoch, read the corresponding ``trainpoints.NNNNNN.out`` /
``trainforces.NNNNNN.out`` (training units if ``normalize_data_set`` was used),
convert back to physical units with the ``mean_energy`` / ``conv_energy`` /
``conv_length`` constants that nnp-train writes into ``input.nn``, and map each
structure index back to the ``tag=`` label in the ``input.data`` comment line
for a per-tag error breakdown.

Unit conversion (n2p2 normalization):
    E_pu/atom = E_iu/atom / conv_energy + mean_energy
    F_pu      = F_iu * conv_length / conv_energy
RMSE-type quantities only need the conv factors (mean_energy cancels in
differences); absolute energies need the full formula.
"""

import re
import numpy as np


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


def read_tags_from_input_data(file_path: str = None, default_tag: str = 'unknown') -> list:
    """Extract one ``tag=`` label per structure from an n2p2 ``input.data``.

    The comment line written by mymetal looks like
    ``comment | tag=A11-1 | frame=1/1 | file=...``; the structure order in the
    file matches the structure index used in trainpoints/trainforces.

    Args:
        file_path: Path to input.data.
        default_tag: Tag used for structures without a tag= comment.

    Returns:
        List of tag strings, one per structure (in file order).
    """
    ltag = []
    tag = default_tag
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            fields = line.split()
            if not fields:
                continue
            if fields[0] == 'begin':
                tag = default_tag
            elif fields[0] == 'comment':
                m = re.search(r'tag=(\S+)', line)
                if m:
                    tag = m.group(1)
            elif fields[0] == 'end':
                ltag.append(tag)
    return ltag


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


def plot_learning_curve(epochs: np.ndarray = None, e_rmse_mev: np.ndarray = None,
                        f_rmse_mev: np.ndarray = None, file_path: str = None) -> None:
    """Plot energy/force training RMSE versus epoch (log-y).

    Args:
        epochs: Epoch numbers.
        e_rmse_mev: Energy RMSE in meV/atom.
        f_rmse_mev: Force RMSE in meV/Angstrom.
        file_path: Output figure path (e.g. learning_curve.pdf).
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(8, 3.2), constrained_layout=True)
    ax[0].semilogy(epochs, e_rmse_mev, lw=0.8, label='train')
    ax[0].set_xlabel('Epoch')
    ax[0].set_ylabel(r'Energy RMSE (meV/atom)')
    ax[1].semilogy(epochs, f_rmse_mev, lw=0.8, label='train')
    ax[1].set_xlabel('Epoch')
    ax[1].set_ylabel(r'Force RMSE (meV/$\mathrm{\AA}$)')
    for a in ax:
        a.grid(ls='--', lw=0.3)
        a.legend()
    fig.savefig(file_path)
    plt.close(fig)


def plot_compare(e_ref: np.ndarray = None, e_nnp: np.ndarray = None, tag_e: np.ndarray = None,
                 f_ref: np.ndarray = None, f_nnp: np.ndarray = None, tag_f: np.ndarray = None,
                 file_path: str = None, text_e: str = '', text_f: str = '') -> None:
    """DFT-vs-NNP scatter plots for energies and forces, colored by tag.

    Args:
        e_ref: DFT energies per atom (eV/atom, physical units).
        e_nnp: NNP energies per atom (eV/atom).
        tag_e: Tag label of each energy point.
        f_ref: DFT force components (eV/Angstrom).
        f_nnp: NNP force components (eV/Angstrom).
        tag_f: Tag label of each force component.
        file_path: Output figure path (e.g. dncompare.pdf).
        text_e: Annotation text for the energy panel (e.g. RMSE summary).
        text_f: Annotation text for the force panel.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    ltag = sorted(set(tag_e))
    cmap = plt.get_cmap('tab20')
    color = {t: cmap(i % 20) for i, t in enumerate(ltag)}

    # 上排全量程；下排按 DFT 值分位数裁剪（压缩二聚体等极端构型会把主体区域压成一个点）
    fig, ax = plt.subplots(2, 2, figsize=(9, 8.4), constrained_layout=True)
    panels = [(ax[0, 0], e_ref, e_nnp, tag_e, None),
              (ax[0, 1], f_ref, f_nnp, tag_f, None),
              (ax[1, 0], e_ref, e_nnp, tag_e, (0.0, 98.0)),
              (ax[1, 1], f_ref, f_nnp, tag_f, (0.5, 99.5))]
    for a, ref, nnp, tags, q in panels:
        if q is None:
            lo, hi = min(ref.min(), nnp.min()), max(ref.max(), nnp.max())
        else:
            lo, hi = np.percentile(ref, q)
        pad = 0.05 * (hi - lo)
        a.plot([lo - pad, hi + pad], [lo - pad, hi + pad], 'k-', lw=0.5)
        for t in ltag:
            mask = tags == t
            if mask.any():
                # rasterized: 几万个点的 pdf 不会爆体积
                a.scatter(ref[mask], nnp[mask], s=2, alpha=0.5, color=color[t],
                          label=t, rasterized=True)
        a.set_xlim(lo - pad, hi + pad)
        a.set_ylim(lo - pad, hi + pad)
        a.set_aspect('equal')
        a.grid(ls='--', lw=0.3)
    ax[0, 0].set_title('Energies')
    ax[0, 1].set_title('Force components')
    ax[1, 0].set_title('Energies (zoom)')
    ax[1, 1].set_title('Force components (zoom)')
    for a in [ax[0, 0], ax[1, 0]]:
        a.set_xlabel(r'DFT energy (eV/atom)')
        a.set_ylabel(r'NNP energy (eV/atom)')
    for a in [ax[0, 1], ax[1, 1]]:
        a.set_xlabel(r'DFT force (eV/$\mathrm{\AA}$)')
        a.set_ylabel(r'NNP force (eV/$\mathrm{\AA}$)')
    ax[0, 0].legend(fontsize=5, ncol=2, loc='lower right', markerscale=2)
    if text_e:
        ax[0, 0].text(0.03, 0.97, text_e, transform=ax[0, 0].transAxes, va='top', fontsize=6)
    if text_f:
        ax[0, 1].text(0.03, 0.97, text_f, transform=ax[0, 1].transAxes, va='top', fontsize=6)
    fig.savefig(file_path, dpi=300)
    plt.close(fig)
