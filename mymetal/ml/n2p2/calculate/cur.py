"""CUR-based symmetry function selection for n2p2.

Workflow: each candidate SF is evaluated separately by ``nnp-scaling`` in its
own directory (input.nn = no-SF header + one dummy zero-value SF + one real SF),
producing a ``function.data`` file. This module collects the per-atom SF values
from all those directories into feature matrices, filters out nearly-zero SFs,
and selects the most informative ones by CUR decomposition (skcosmo).
"""

import numpy as np
from pathlib import Path


def _read_function_data(file_path: str = None) -> tuple:
    """Parse one n2p2 ``function.data`` file (written by nnp-scaling).

    File format per frame: one line with the number of atoms, then one line per
    atom (``element_Z sf_value_1 sf_value_2 ...``), then one closing line with
    charge/energy information.

    Args:
        file_path: Path to the function.data file.

    Returns:
        A tuple ``(lsp_frm, lfeat_frm)`` where ``lsp_frm`` is a list (one entry
        per frame) of integer arrays with the atomic numbers, and ``lfeat_frm``
        is a list (one entry per frame) of per-atom lists of SF values.
    """
    # lsp_frm: [[79, 79], []]  every frame
    # lfeat_frm: [[[0.0, 1.0], [0.0, 1.0]], []]  every frame
    lsp_frm = []
    lfeat_frm = []
    sp = []
    feat = []
    nat = None
    counter = 0

    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            fields = line.split()
            if len(fields) == 1:
                # first line of a frame: number of atoms
                nat = int(fields[0])
                counter = 0
                sp = []
                feat = []
            elif counter < nat:
                sp.append(int(fields[0]))
                feat.append([float(s) for s in fields[1:]])
                counter += 1
            else:
                # closing line of a frame
                lsp_frm.append(np.array(sp))
                lfeat_frm.append(feat)

    return lsp_frm, lfeat_frm


def collect_sf_features(ldir: list = None, fname: str = 'function.data') -> tuple:
    """Collect SF feature matrices from a list of nnp-scaling directories.

    Each directory corresponds to one candidate SF. The per-atom SF value is
    taken as ``max`` over the SF columns: the dummy SF is always ~0 and real SF
    values are non-negative, so the max picks the real one regardless of the
    column order n2p2 uses internally. Atoms of elements that do not carry the
    SF columns are marked with -1 (only relevant for multi-element systems).

    Args:
        ldir: List of directories, ordered consistently with the SF list.
        fname: Name of the n2p2 function file inside each directory.

    Returns:
        A tuple ``(feat_atom, feat_av)`` where ``feat_atom`` has shape
        ``(n_atoms_total, n_sf)`` (per-atom SF values, -1 for other elements)
        and ``feat_av`` has shape ``(n_frames, n_sf)`` (per-frame averages).
    """
    lcol_atom = []
    lcol_av = []

    for d in ldir:
        lsp_frm, lfeat_frm = _read_function_data(Path(d) / fname)

        # find_sp 这段其实比较奇怪，以后对于多元素势函数或许需要更稳健的方法来筛选
        find_sp = None
        for sp_frm, feat_frm in zip(lsp_frm, lfeat_frm):
            for s, v in zip(sp_frm, feat_frm):
                if len(v) == 2:
                    find_sp = s
                    break
            if find_sp is not None:
                break
        if find_sp is None:
            raise ValueError(f"❌ No atom with 2 SF columns found in {Path(d) / fname}, "
                             f"cannot identify the center element of the real SF.")

        col_atom = []   # [1, 2, 3, 4], if 1, 2 is a frame, 3, 4 is another frame
        col_av = []     # [1.5,   3.5]
        for sp_frm, feat_frm in zip(lsp_frm, lfeat_frm):
            sf_sum = 0.0
            for s, v in zip(sp_frm, feat_frm):
                if s == find_sp:
                    sf = max(v)
                    col_atom.append(sf)
                    sf_sum += sf
                else:
                    col_atom.append(-1.0)
            # 每一个frame内做了一个平均
            col_av.append(sf_sum / len(sp_frm))

        lcol_atom.append(np.array(col_atom))
        lcol_av.append(np.array(col_av))

    # 转置操作：将lcol_atom 和 lcol_av 的每一行向量转换成列向量排列
    feat_atom = np.column_stack(lcol_atom)
    feat_av = np.column_stack(lcol_av)
    return feat_atom, feat_av


def filter_zero_columns(feat_atom: np.ndarray = None, max_zero_frac: float = 0.05,
                        zero_atol: float = 0.0, min_std: float = 0.0) -> tuple:
    """Drop SF columns that carry no information for this dataset.

    Two failure modes are removed here:

    1. **Too many zeros** — a SF that is zero for most atoms has an (almost)
       empty cutoff sphere on the dataset (``n_zero / n_valid > max_zero_frac``).
    2. **Near-constant / low variance** — a SF whose per-atom value barely
       changes across the whole dataset (``std <= min_std``) is *insensitive to
       the dataset*: it carries no discriminative information regardless of
       whether its constant value is zero or not. The zero-fraction test alone
       cannot catch a SF that is constant at a non-zero value, so this is the
       criterion that actually screens "对数据集不敏感" SFs.

    Args:
        feat_atom: Per-atom feature matrix ``(n_atoms_total, n_sf)``, -1 marks
            atoms of other elements and is excluded from the statistics.
        max_zero_frac: Maximum allowed fraction of zero values per column.
        zero_atol: Values with ``abs(v) <= zero_atol`` count as zero. The
            default 0.0 reproduces the reference behavior (exact zeros only,
            i.e. printed as 0.0000000000 in function.data); ~1e-9 additionally
            counts cutoff-tail values that are informationally zero. Must be
            < 1 so the -1 markers are never matched.
        min_std: Minimum per-column standard deviation (over valid, non -1
            entries) to keep a SF. 0.0 disables the variance screen (backward
            compatible); a small value such as 1e-4 removes near-constant SFs.

    Returns:
        A tuple ``(kept_idx, dropped_idx)`` of integer index arrays.
    """
    kept_idx = []
    dropped_idx = []
    for i in range(feat_atom.shape[1]):
        col = feat_atom[:, i]
        valid = col[col != -1]
        # 所有的SF数量
        n_valid = valid.size
        # 为 0 的SF数量
        n_zero = np.count_nonzero(np.abs(valid) <= zero_atol)
        std = float(valid.std()) if n_valid > 0 else 0.0
        if n_valid == 0 or n_zero / n_valid > max_zero_frac or std <= min_std:
            dropped_idx.append(i)
        else:
            kept_idx.append(i)
    return np.array(kept_idx, dtype=int), np.array(dropped_idx, dtype=int)


def cur_select(feat: np.ndarray = None, n_select: int = 48, standardize: bool = True) -> np.ndarray:
    """Select the most informative feature columns by CUR decomposition.

    Args:
        feat: Feature matrix ``(n_samples, n_sf)``, e.g. per-atom SF values.
            Rows that are entirely -1 (atoms of other elements) are dropped
            before the decomposition.
        n_select: Number of columns to select.
        standardize: If True (default) each column is z-scored (centered and
            scaled to unit std) before CUR. CUR ranks columns by their weight in
            the dominant singular subspace of the *uncentered* matrix, so on raw
            data a high-magnitude but low-variance (dataset-insensitive) SF can
            be favoured through its mean offset. Centering removes that bias so
            selection follows variation, not absolute level; scaling makes SFs
            of different magnitude comparable. Columns with zero std are left at
            zero (they should already be removed by ``filter_zero_columns``).

    Returns:
        Integer array with the indices of the selected columns.
    """
    # lazy import: keep the module usable without skcosmo installed
    from skcosmo.feature_selection import CUR

    X = np.asarray(feat, dtype=float)
    # 丢掉整行都是 -1 的原子（非中心元素），避免 -1 标记污染分解
    row_valid = ~np.all(X == -1, axis=1)
    X = X[row_valid]

    if standardize:
        mu = X.mean(axis=0)
        sd = X.std(axis=0)
        nz = sd > 0
        Xs = np.zeros_like(X)
        Xs[:, nz] = (X[:, nz] - mu[nz]) / sd[nz]
        X = Xs

    y = np.random.rand(X.shape[0])  # required by the fit API, not used by CUR
    selected_idx = np.asarray(CUR(n_to_select=int(n_select)).fit(X, y).selected_idx_, dtype=int)

    # 旧版 skcosmo 在特征退化时（重复列、n_select 超过有效秩）会重选同一索引；
    # 参考模板用 set() 静默去重（结果会少于 n_select），这里改为显式报错
    if len(np.unique(selected_idx)) != len(selected_idx):
        raise ValueError(f"❌ CUR returned duplicate indices: n_select={n_select} likely exceeds "
                         f"the number of linearly independent features. Reduce n_select.")
    return selected_idx
