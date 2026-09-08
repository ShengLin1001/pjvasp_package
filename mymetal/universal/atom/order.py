"""
order submodule

Local bond-order (Steinhardt) parameters used to ask "is this frame still the
phase I built it as?". An MD snapshot keeps its file name after the structure
underneath has changed, so the phase label of a trajectory is an assumption
until something measures it; ``q4``/``q6`` measure it from the positions alone,
without reference lattice sites or atom mapping.

Reference values for the 12-coordinated close-packed structures (ideal, 0 K):
FCC ``q4 = 0.191``, ``q6 = 0.575``; HCP ``q4 = 0.097``, ``q6 = 0.485``. FCC and
HCP are far apart in ``q4`` and close in ``q6``, so ``q4`` is the discriminating
one; a liquid drops both (``q6 < 0.3``).

Functions:
    - get_closepacked_cutoff: First-shell neighbour cutoff of a close-packed cell.
    - get_steinhardt_ql: Per-atom Steinhardt q_l, averaged over the structure.
    - get_closepacked_phase: Classify a close-packed frame as fcc / hcp / other.
"""

import numpy as np
from ase import Atoms
from ase.neighborlist import neighbor_list
from scipy.special import sph_harm

# 理想值（0 K、12 配位）；判据用 q4，因为 FCC/HCP 在 q4 上差一倍而在 q6 上很近
QL_IDEAL = {'fcc': {4: 0.191, 6: 0.575}, 'hcp': {4: 0.097, 6: 0.485}}


def get_closepacked_cutoff(atoms: Atoms = None, n_shell: int = 12,
                           margin: float = 1.15, rmax: float = 6.0) -> float:
    """Neighbour cutoff that keeps a close-packed frame 12-coordinated.

    Anchoring on the shortest bond in the cell is fragile: one hot pair of
    atoms drags the cutoff below the first shell and leaves other atoms with no
    neighbours at all. The 12th-nearest distance, taken as a median over atoms,
    is the first shell itself and does not care about the tail of the bond
    distribution.

    Args:
        atoms (Atoms): Structure with a periodic cell.
        n_shell (int): Coordination the cutoff should enclose (12 for FCC/HCP).
        margin (float): Multiplies that distance, leaving room between the first
            and second shells.
        rmax (float): Search radius for the underlying neighbour list.

    Returns:
        float: Cutoff in Angstrom.

    Raises:
        ValueError: If the cell has no neighbours within ``rmax``.
    """
    lindex, ldist = neighbor_list('id', atoms, rmax)
    if len(ldist) == 0:
        raise ValueError('no neighbours within %.1f A; is the cell sane?' % rmax)
    lshell = []
    for index in range(len(atoms)):
        lone = np.sort(ldist[lindex == index])
        if len(lone) >= n_shell:
            lshell.append(lone[n_shell - 1])
    if not lshell:
        raise ValueError('no atom has %d neighbours within %.1f A' % (n_shell, rmax))
    return float(margin * np.median(lshell))


def get_steinhardt_ql(atoms: Atoms = None, l: int = 4, cutoff: float = None,
                      if_average: bool = True, if_per_atom: bool = False):
    """Steinhardt bond-order parameter ``q_l`` of a structure.

    Args:
        atoms (Atoms): Structure with a periodic cell.
        l (int): Spherical-harmonic degree; 4 separates FCC from HCP, 6 measures
            crystallinity as such.
        cutoff (float): Neighbour cutoff in Angstrom. ``None`` puts it between
            the first and second close-packed shells (1.2 x the shortest bond),
            which is what keeps the coordination at 12 for a thermal snapshot.
        if_average (bool): Use the averaged form of Lechner and Dellago
            (J. Chem. Phys. 129, 114707), i.e. average ``q_lm`` over the
            neighbour shell *and* the atom itself before taking the modulus.
            The raw parameter smears badly once the frame is hot -- a 600 K FCC
            snapshot drifts far enough in raw ``q4`` to be mistaken for another
            structure -- while the averaged one stays on its 0 K value.
        if_per_atom (bool): Return the per-atom values instead of their mean.

    Returns:
        float or np.ndarray: Mean ``q_l`` over atoms, or the per-atom array.

    Raises:
        ValueError: If the structure has no atoms or no neighbours in range.
    """
    if atoms is None or len(atoms) == 0:
        raise ValueError('atoms is required and must be non-empty')
    if cutoff is None:
        cutoff = get_closepacked_cutoff(atoms)

    lindex, lneighbor, lvector = neighbor_list('ijD', atoms, cutoff)
    lqlm = np.zeros((len(atoms), 2 * l + 1), dtype=complex)
    for index in range(len(atoms)):
        lvec = lvector[lindex == index]
        if len(lvec) == 0:
            raise ValueError('atom %d has no neighbour within %.3f A' % (index, cutoff))
        lr = np.linalg.norm(lvec, axis=1)
        # scipy 的 sph_harm(m, l, azimuthal, polar) —— 方位角在前，别写反
        ltheta = np.arccos(np.clip(lvec[:, 2] / lr, -1.0, 1.0))
        lphi = np.arctan2(lvec[:, 1], lvec[:, 0])
        lqlm[index] = [sph_harm(m, l, lphi, ltheta).mean() for m in range(-l, l + 1)]

    if if_average:
        # Lechner-Dellago：先把 q_lm 在「自身 + 近邻」上平均，再取模。
        # 热噪声在这一步被邻域平均掉，判据才在几百 K 下还站得住。
        lqlm_avg = np.zeros_like(lqlm)
        for index in range(len(atoms)):
            lnb = lneighbor[lindex == index]
            lqlm_avg[index] = (lqlm[lnb].sum(axis=0) + lqlm[index]) / (len(lnb) + 1)
        lqlm = lqlm_avg
    lql = np.sqrt(4.0 * np.pi / (2 * l + 1) * (np.abs(lqlm) ** 2).sum(axis=1))
    return lql if if_per_atom else float(lql.mean())


def get_closepacked_phase(atoms: Atoms = None, cutoff: float = None,
                          tol: float = 0.25, if_average: bool = True) -> tuple:
    """Classify a close-packed frame as ``fcc`` / ``hcp`` / ``other``.

    The verdict is on ``q4`` alone (see module docstring); ``q6`` only guards
    against calling a molten or heavily defective frame a crystal.

    Args:
        atoms (Atoms): Structure to classify.
        cutoff (float): Neighbour cutoff passed to :func:`get_steinhardt_ql`.
        if_average (bool): Use the averaged (Lechner-Dellago) parameters; keep
            this on for thermal snapshots.
        tol (float): Half-width of each phase window, as a fraction of the
            ``q4`` gap between the ideal FCC and HCP values. 0.25 puts the
            decision boundary at the midpoint with a quarter-gap margin on
            either side, so a frame caught mid-transformation reads ``other``
            instead of being forced into one of the two.

    Returns:
        tuple: ``(phase, q4, q6)``.
    """
    q4 = get_steinhardt_ql(atoms, l=4, cutoff=cutoff, if_average=if_average)
    q6 = get_steinhardt_ql(atoms, l=6, cutoff=cutoff, if_average=if_average)
    q4_fcc, q4_hcp = QL_IDEAL['fcc'][4], QL_IDEAL['hcp'][4]
    gap = abs(q4_fcc - q4_hcp)
    if q6 < 0.3:                       # 晶体性都没有了，谈 FCC/HCP 无意义
        return 'other', q4, q6
    if abs(q4 - q4_fcc) <= tol * gap:
        return 'fcc', q4, q6
    if abs(q4 - q4_hcp) <= tol * gap:
        return 'hcp', q4, q6
    return 'other', q4, q6
