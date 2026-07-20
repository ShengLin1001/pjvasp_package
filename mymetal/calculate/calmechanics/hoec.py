"""
hoec module

Group-theory core for higher-order elastic constants (HOEC: 2nd + 3rd + 4th order)
by the energy-strain method.

Reference method: Hao Wang & Mo Li, "Ab initio calculations of second-, third-,
and fourth-order elastic constants for single crystals", Phys. Rev. B 79, 224102
(2009). That paper only treats cubic crystals; here the same energy-strain scheme
is generalized to arbitrary point-group symmetry by deriving the reduced strain-
energy density directly from the point group.

Physics:
    The strain-energy density is expanded in the engineering-Voigt Green-Lagrange
    strain eta = (eta1..eta6) with eta4=2*eps_yz, eta5=2*eps_xz, eta6=2*eps_xy::

        u(eta) = 1/2! C_IJ eta_I eta_J + 1/3! C_IJK eta_I eta_J eta_K
                 + 1/4! C_IJKL eta_I eta_J eta_K eta_L + ...

    with fully index-symmetric Brugger constants C. Point-group invariance
    u(M(g) eta) = u(eta) for every rotation g reduces the constants to the
    independent set (cubic: 3/6/11, hexagonal: 5/10/19).

    For a single-parameter deformation mode eta = xi * d the energy is a polynomial::

        [U(xi) - U(0)] / V0 = P2 xi^2 + P3 xi^3 + P4 xi^4 + O(xi^5),

    and P_n = sum_k B_n_k(d) * C_n_k, i.e. each P coefficient is a fixed linear
    combination of the independent order-n constants (the analog of Wang-Li Table I).
    This module produces those combinations for any mode and any symmetry, so the
    structure generator and the post-processing solve exactly the same linear system.

Functions:
    - rot_z: rotation matrix about z.
    - rot_axis: rotation matrix about an arbitrary axis.
    - voigt_strain_M: 6x6 engineering-Voigt strain transformation for a rotation.
    - close_group: close a set of generators into the full rotation group.
    - monomials: index multisets of the order-n strain monomials.
    - mono_eval_vec: evaluate all monomials at a strain vector.
    - voigt_name: index multiset -> Brugger Voigt name.
    - n_orderings: distinct orderings of an index multiset.
    - rationalize: float matrix -> exact Fraction matrix.
    - reduce_symmetry: reduce order-n strain energy to independent Brugger constants.
    - HOECModel: reduced 2nd/3rd/4th-order elastic model for one symmetry.
    - get_model: cached HOECModel.
    - get_strain_list: symmetric list of xi values including 0.
    - get_deformation_gradient: symmetric F with Green-Lagrange strain xi * d.
    - get_mode_severity: how hard a mode has actually deformed the cell at amplitude xi.
    - get_mode_window: largest |xi| a mode may take within a severity budget.
    - get_shear_window: largest |xi| keeping a mode's tensor shear within a cap.
    - scale_mode_shear: multiply a mode's engineering shear entries by a factor (small-shear).
    - get_shear_modes: names of the shear-carrying core modes (small-shear default target).
    - get_hoec_modes: core mode table, optionally shear-rescaled and/or with extra normal modes.
    - get_mode_strain_lists: per-mode xi lists of equal severity (or one shared window).
    - check_symmetry: validate / auto-detect 'cubic' or 'hex' from a cell.
    - selftest_hoec: verify cubic against Wang-Li Table I and check mode-set rank.

Change log:
    - Written by J. P., 2026, for the Cu/Ag/Au FCC+HCP higher-order elastic-constant
      workflow (``hoec_core.py``).
    - Merged into mymetal by J. P. on 2026-07-10; the deformation gradient, the
      strain list and the symmetry detection moved in from the generator script.
    - Revised by J. P. on 2026-07-14 with the per-mode severity budget
      (``get_mode_severity`` / ``get_mode_window`` / ``get_mode_strain_lists``). One shared
      xi window over-strains the multi-component modes badly enough to destroy the
      constants they alone determine -- for Au at xi=0.12, mode H reached -27.9% volume
      (4.1 eV/atom, above the cohesive energy) and returned C1266 = 1891 GPa against a
      literature 9402.
    - Fixed by Codex on 2026-07-14: project elastic coefficients with the dual monomial
      representation, add hex mode M20 so the fourth-order system has rank 19/19, and
      validate hex rotational invariance plus the analytic M04/M06 SOEC identity.
    - Revised by J. P. on 2026-07-18: added the opt-in hex small-shear option
      (``get_hoec_modes`` / ``scale_mode_shear`` / ``get_shear_modes`` / ``get_shear_window``) so a
      relaxed-ion HCP run can keep the shear directions without following the basal shuffle into
      FCC. It rescales the selected core modes' engineering shear IN PLACE (default 0.5, target =
      every shear-carrying mode), keeping their names, and caps their window at
      ``SMALL_SHEAR_SHEAR_CAP``. No mode is added or renamed, so mymetal.post.hoec_energy must
      read each mode's direction from the run manifest, not from the global ``MODES`` table.
"""

import numpy as np
from fractions import Fraction
from itertools import combinations_with_replacement
from math import factorial
from collections import Counter

VOIGT = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]  # Voigt 1..6 -> tensor ij


# ---------------------------------------------------------------------------
# rotations and the engineering-Voigt strain transformation
# ---------------------------------------------------------------------------
def rot_z(t: float = None) -> np.ndarray:
    """Rotation matrix about the z axis.

    Args:
        t (float): Rotation angle in radians.

    Returns:
        np.ndarray: The 3x3 rotation matrix.
    """
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])


def rot_axis(axis: list = None, t: float = None) -> np.ndarray:
    """Rotation matrix about an arbitrary axis (Rodrigues' formula).

    Args:
        axis (list): Rotation axis, normalized internally.
        t (float): Rotation angle in radians.

    Returns:
        np.ndarray: The 3x3 rotation matrix.
    """
    ax = np.array(axis, float)
    ax = ax / np.linalg.norm(ax)
    x, y, z = ax
    c, s = np.cos(t), np.sin(t)
    C = 1 - c
    return np.array([
        [c + x * x * C, x * y * C - z * s, x * z * C + y * s],
        [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
        [z * x * C - y * s, z * y * C + x * s, c + z * z * C]])


def voigt_strain_M(R: np.ndarray = None) -> np.ndarray:
    """6x6 matrix mapping the engineering-Voigt strain under eps -> R eps R^T.

    Args:
        R (np.ndarray): A 3x3 rotation matrix.

    Returns:
        np.ndarray: The 6x6 engineering-Voigt strain transformation.
    """
    M = np.zeros((6, 6))
    for j in range(6):
        E = np.zeros((3, 3))
        a, b = VOIGT[j]
        val = 1.0 if a == b else 0.5     # engineering shear carries a 1/2 in the tensor
        E[a, b] = val
        if a != b:
            E[b, a] = val
        Ep = R @ E @ R.T
        M[:, j] = [Ep[0, 0], Ep[1, 1], Ep[2, 2], 2 * Ep[1, 2], 2 * Ep[0, 2], 2 * Ep[0, 1]]
    return M


# generators of the proper rotation subgroups (enough for elastic-tensor symmetry)
GENERATORS = {
    'cubic': [rot_z(np.pi / 2), rot_axis([1, 1, 1], 2 * np.pi / 3)],       # O, order 24
    'hex':   [rot_z(np.pi / 3), rot_axis([1, 0, 0], np.pi)],               # D6, order 12
}


def close_group(lgen: list = None, tol: float = 1e-9, maxn: int = 2000) -> list:
    """Close a set of generators into the full rotation group.

    Args:
        lgen (list): List of 3x3 generator matrices.
        tol (float): Tolerance for deciding whether two matrices coincide.
        maxn (int): Safety cap on the group order.

    Returns:
        list: The group elements as 3x3 matrices.
    """
    lgroup = [np.eye(3)]

    def known(A):
        return any(np.allclose(A, G, atol=tol) for G in lgroup)

    for g in lgen:
        if not known(g):
            lgroup.append(g)
    changed = True
    while changed and len(lgroup) < maxn:
        changed = False
        for A in list(lgroup):
            for g in lgen:
                B = A @ g
                if not known(B):
                    lgroup.append(B)
                    changed = True
    return lgroup


# ---------------------------------------------------------------------------
# monomial basis and reduction to the independent Brugger constants
# ---------------------------------------------------------------------------
def monomials(order: int = None) -> list:
    """Index multisets of the order-n strain monomials.

    Args:
        order (int): Expansion order (2, 3 or 4).

    Returns:
        list: Index multisets, e.g. (0, 0) for eta1^2.
    """
    return list(combinations_with_replacement(range(6), order))


def mono_eval_vec(vec: np.ndarray = None, lmono: list = None) -> np.ndarray:
    """Evaluate all monomials at one strain vector.

    Args:
        vec (np.ndarray): Length-6 engineering-Voigt strain.
        lmono (list): Monomial index multisets from :func:`monomials`.

    Returns:
        np.ndarray: Monomial values, aligned with ``lmono``.
    """
    out = np.ones(len(lmono))
    for k, m in enumerate(lmono):
        p = 1.0
        for idx in m:
            p *= vec[idx]
        out[k] = p
    return out


def voigt_name(m: tuple = None) -> str:
    """Index multiset -> Brugger Voigt name, e.g. (0, 0) -> '11'.

    Args:
        m (tuple): Monomial index multiset (0-based).

    Returns:
        str: The 1-based Voigt name.
    """
    return ''.join(str(i + 1) for i in m)


def n_orderings(m: tuple = None) -> int:
    """Distinct orderings of an index multiset (multinomial coefficient).

    Args:
        m (tuple): Monomial index multiset.

    Returns:
        int: Number of distinct orderings.
    """
    num = factorial(len(m))
    for c in Counter(m).values():
        num //= factorial(c)
    return num


def rationalize(A: np.ndarray = None, maxden: int = 4620) -> np.ndarray:
    """Convert a float matrix to exact Fractions.

    The elastic reduction coefficients are small rationals; 4620 = lcm(1..11)/... is
    a denominator bound comfortably above every coefficient that occurs up to 4th order.

    Args:
        A (np.ndarray): Float matrix.
        maxden (int): Maximum denominator passed to ``limit_denominator``.

    Returns:
        np.ndarray: Object array of ``Fraction`` with the same shape.
    """
    out = np.empty(A.shape, dtype=object)
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            out[i, j] = Fraction(A[i, j]).limit_denominator(maxden)
    return out


def reduce_symmetry(symmetry: str = None, order: int = None,
                    lprefer: list = None, seed: int = 1) -> tuple:
    """Reduce the order-n strain-energy density to independent Brugger constants.

    Args:
        symmetry (str): 'cubic' or 'hex'.
        order (int): 2, 3, or 4.
        lprefer (list): Voigt names to select as the independent set first, so the
            naming matches a chosen textbook convention.
        seed (int): RNG seed for the sample points defining the monomial representation.

    Returns:
        tuple: ``(lmono, R, lname_indep, lpivot)`` where ``lmono`` are the monomial
        index multisets, ``R`` (nmono x d, Fractions) satisfies ``c = R @ theta`` with
        ``c_m`` the coefficient of monomial m in u and ``theta`` the independent
        *Brugger* constants, ``lname_indep`` are the d Voigt names, and ``lpivot``
        are the monomial indices chosen as the independent set.

    Raises:
        ValueError: If ``symmetry`` is not a supported point group.
    """
    if symmetry not in GENERATORS:
        raise ValueError("symmetry must be 'cubic' or 'hex', got %r" % symmetry)
    lgroup = close_group(GENERATORS[symmetry])
    lmono = monomials(order)
    n = len(lmono)
    rng = np.random.default_rng(seed)
    S = rng.standard_normal((n, 6))                # well-conditioned sample points
    V0 = np.array([mono_eval_vec(S[i], lmono) for i in range(n)])
    V0inv = np.linalg.inv(V0)

    ### average the induced representation to the group projector -----------------
    Pr = np.zeros((n, n))
    for R3 in lgroup:
        M = voigt_strain_M(R3)
        Vg = np.array([mono_eval_vec(M @ S[i], lmono) for i in range(n)])
        Tg = (V0inv @ Vg).T                        # mono_m(M s) = sum_m' T[m,m'] mono_m'(s)
        # Energy is q(s)^T c, so q(Ms)=Tq(s) requires T^T c=c. The monomials and
        # their coefficients belong to dual representations; projecting with T itself
        # happens to work for cubic here but violates hex invariance at every order.
        Pr += Tg.T
    Pr /= len(lgroup)
    # to here: Pr projects the monomial-coefficient space onto invariant vectors

    ### invariant space = null space of (Pr - I), robust via SVD ------------------
    U, sv, Vt = np.linalg.svd(Pr - np.eye(n))
    d = int(np.sum(sv < 1e-7 * max(sv[0], 1.0)))
    B = Vt[len(sv) - d:, :].T                      # (n x d) orthonormal invariant basis

    ### fix the naming of the independent set: preferred names first, then lex -----
    dict_name2idx = {voigt_name(m): i for i, m in enumerate(lmono)}
    ltry = []
    if lprefer:
        for nm in lprefer:
            if nm in dict_name2idx:
                ltry.append(dict_name2idx[nm])
    ltry += [k for k in range(n) if k not in ltry]
    lpivot = []
    for k in ltry:
        ltrial = lpivot + [k]
        s = np.linalg.svd(B[ltrial, :], compute_uv=False)
        if s[-1] > 1e-8 * s[0]:                    # rows still independent
            lpivot.append(k)
        if len(lpivot) == d:
            break

    ### express every coefficient via the pivots, rescale pivots to Brugger --------
    R = B @ np.linalg.inv(B[lpivot, :])            # c = R @ (c at pivots)
    alpha = np.array([n_orderings(lmono[p]) / factorial(order) for p in lpivot])
    R = R * alpha[None, :]                         # theta now = Brugger constants
    R = rationalize(R)
    lname_indep = [voigt_name(lmono[p]) for p in lpivot]
    return lmono, R, lname_indep, lpivot


# ---------------------------------------------------------------------------
# textbook independent sets (fix the pivot naming)
# ---------------------------------------------------------------------------
PREFER = {
    'cubic': {
        2: ['11', '12', '44'],
        3: ['111', '112', '123', '144', '155', '456'],
        4: ['1111', '1112', '1122', '1123', '1144', '1155',
            '1255', '1266', '1456', '4444', '4455'],
    },
    'hex': {
        2: ['11', '12', '13', '33', '44'],
        3: ['111', '112', '113', '123', '133', '144', '155', '222', '333', '344'],
        4: None,   # 19 independent, auto-named lexicographically
    },
}


# ---------------------------------------------------------------------------
# deformation-mode tables (single-parameter Voigt strain directions d)
# ---------------------------------------------------------------------------
# The following formulas come directly from the continuum expansion
#
#   u(xi) - u(0) = P2*xi^2 + P3*xi^3 + P4*xi^4 + O(xi^5),  eta = xi*d,
#
# after imposing the crystal symmetry on
#
#   u(eta) = C_IJ*eta_I*eta_J/2!
#          + C_IJK*eta_I*eta_J*eta_K/3!
#          + C_IJKL*eta_I*eta_J*eta_K*eta_L/4! + ... .
#
# Here d is the row vector [eta1, eta2, eta3, eta4, eta5, eta6]/xi; eta4--eta6
# are engineering shear strains. The formulas are kept beside the mode definitions
# so a generated strain path can be audited without reverse-engineering the
# symmetry-reduction matrices.
#
# Cubic modes (independent constants follow PREFER['cubic']):
#
# A: d = [1, 0, 0, 0, 0, 0]
#    P2 = 1/2*C11
#    P3 = 1/6*C111
#    P4 = 1/24*C1111
# B: d = [1, 1, 0, 0, 0, 0]
#    P2 = C11 + C12
#    P3 = 1/3*C111 + C112
#    P4 = 1/12*C1111 + 1/3*C1112 + 1/4*C1122
# C: d = [1, -1, 0, 0, 0, 0]
#    P2 = C11 - C12
#    P3 = 0
#    P4 = 1/12*C1111 - 1/3*C1112 + 1/4*C1122
# D: d = [1, 0, 0, 2, 0, 0]
#    P2 = 1/2*C11 + 2*C44
#    P3 = 1/6*C111 + 2*C144
#    P4 = 1/24*C1111 + C1144 + 2/3*C4444
# E: d = [1, 0, 0, 0, 0, 2]
#    P2 = 1/2*C11 + 2*C44
#    P3 = 1/6*C111 + 2*C155
#    P4 = 1/24*C1111 + C1155 + 2/3*C4444
# F: d = [0, 0, 0, 2, 2, 2]
#    P2 = 6*C44
#    P3 = 8*C456
#    P4 = 2*C4444 + 12*C4455
# G: d = [0, 0, 0, 2, 0, 0]
#    P2 = 2*C44
#    P3 = 0
#    P4 = 2/3*C4444
# H: d = [1, 1, 0, 0, 0, 2]
#    P2 = C11 + C12 + 2*C44
#    P3 = 1/3*C111 + C112 + 4*C155
#    P4 = 1/12*C1111 + 1/3*C1112 + 1/4*C1122 + 2*C1155
#         + 2*C1266 + 2/3*C4444
# I: d = [1, 1, 0, 2, 0, 0]
#    P2 = C11 + C12 + 2*C44
#    P3 = 1/3*C111 + C112 + 2*C144 + 2*C155
#    P4 = 1/12*C1111 + 1/3*C1112 + 1/4*C1122 + C1144 + C1155
#         + 2*C1255 + 2/3*C4444
# J: d = [1, 0, 0, 2, 2, 2]
#    P2 = 1/2*C11 + 6*C44
#    P3 = 1/6*C111 + 2*C144 + 4*C155 + 8*C456
#    P4 = 1/24*C1111 + C1144 + 2*C1155 + 8*C1456 + 2*C4444
#         + 12*C4455
# K: d = [1, 1, 1, 0, 0, 0]
#    P2 = 3/2*C11 + 3*C12
#    P3 = 1/2*C111 + 3*C112 + C123
#    P4 = 1/8*C1111 + C1112 + 3/4*C1122 + 3/2*C1123
#
# All 11 cubic modes are required by the 11-dimensional fourth-order system:
# removing any one of A--K lowers its rank.
# cubic: the 11 modes A..K of Wang-Li Table I.
MODES_CUBIC = {
    'A': (1, 0, 0, 0, 0, 0),
    'B': (1, 1, 0, 0, 0, 0),
    'C': (1, -1, 0, 0, 0, 0),
    'D': (1, 0, 0, 2, 0, 0),
    'E': (1, 0, 0, 0, 0, 2),
    'F': (0, 0, 0, 2, 2, 2),
    'G': (0, 0, 0, 2, 0, 0),
    'H': (1, 1, 0, 0, 0, 2),
    'I': (1, 1, 0, 2, 0, 0),
    'J': (1, 0, 0, 2, 2, 2),
    'K': (1, 1, 1, 0, 0, 0),
}
# Hexagonal modes (independent constants follow PREFER['hex']; the fourth-order
# names are the lexicographic independent set reported by HOECModel.names(4)):
#
# M01: d = [1, 0, 0, 0, 0, 0]
#      P2 = 1/2*C11
#      P3 = 1/6*C111
#      P4 = 1/24*C1111
# M02: d = [0, 0, 1, 0, 0, 0]
#      P2 = 1/2*C33
#      P3 = 1/6*C333
#      P4 = 1/24*C3333
# M03: d = [1, 1, 0, 0, 0, 0]
#      P2 = C11 + C12
#      P3 = 2/3*C111 + C112 - 1/3*C222
#      P4 = 11/108*C1111 + 17/54*C1112 + 1/4*C1122 - 1/9*C1166
# M04: d = [1, -1, 0, 0, 0, 0]
#      P2 = C11 - C12
#      P3 = 2/3*C111 - 2/3*C222
#      P4 = 1/36*C1111 - 5/18*C1112 + 1/4*C1122 + 1/3*C1166
# M05: d = [0, 0, 0, 2, 0, 0]
#      P2 = 2*C44
#      P3 = 0
#      P4 = 2/3*C4444
# M06: d = [0, 0, 0, 0, 0, 2]
#      P2 = C11 - C12
#      P3 = 0
#      P4 = 1/36*C1111 - 5/18*C1112 + 1/4*C1122 + 1/3*C1166
# M07: d = [1, 0, 1, 0, 0, 0]
#      P2 = 1/2*C11 + C13 + 1/2*C33
#      P3 = 1/6*C111 + 1/2*C113 + 1/2*C133 + 1/6*C333
#      P4 = 1/24*C1111 + 1/6*C1113 + 1/4*C1133 + 1/6*C1333
#           + 1/24*C3333
# M08: d = [1, 1, 1, 0, 0, 0]
#      P2 = C11 + C12 + 2*C13 + 1/2*C33
#      P3 = 2/3*C111 + C112 + C113 + C123 + C133 - 1/3*C222
#           + 1/6*C333
#      P4 = 11/108*C1111 + 17/54*C1112 + 1/3*C1113 + 1/4*C1122
#           + 2/3*C1123 + 1/2*C1133 - 1/9*C1166 + 1/3*C1223
#           + 1/2*C1233 + 1/3*C1333 + 1/24*C3333
# M09: d = [1, 0, 0, 2, 0, 0]
#      P2 = 1/2*C11 + 2*C44
#      P3 = 1/6*C111 + 2*C144
#      P4 = 1/24*C1111 + C1144 + 2/3*C4444
# M10: d = [0, 0, 1, 2, 0, 0]
#      P2 = 1/2*C33 + 2*C44
#      P3 = 1/6*C333 + 2*C344
#      P4 = 1/24*C3333 + C3344 + 2/3*C4444
# M11: d = [0, 0, 1, 0, 0, 2]
#      P2 = C11 - C12 + 1/2*C33
#      P3 = C113 - C123 + 1/6*C333
#      P4 = 1/36*C1111 - 5/18*C1112 + 1/4*C1122 + 1/2*C1133
#           + 1/3*C1166 - 1/2*C1233 + 1/24*C3333
# M12: d = [0, 0, 0, 2, 0, 2]
#      P2 = C11 - C12 + 2*C44
#      P3 = 0
#      P4 = 1/36*C1111 - 5/18*C1112 + 1/4*C1122 + C1144 + C1155
#           + 1/3*C1166 + C1244 - 3*C1255 + 2/3*C4444
# M13: d = [0, 0, 0, 2, 2, 2]
#      P2 = C11 - C12 + 4*C44
#      P3 = -4*C144 + 4*C155
#      P4 = 1/36*C1111 - 5/18*C1112 + 1/4*C1122 + 2*C1144
#           + 2*C1155 + 1/3*C1166 - 2*C1244 - 2*C1255 + 8/3*C4444
# M14: d = [1, 0, 0, 2, 2, 0]
#      P2 = 1/2*C11 + 4*C44
#      P3 = 1/6*C111 + 2*C144 + 2*C155
#      P4 = 1/24*C1111 + C1144 + C1155 + 8/3*C4444
# M15: d = [1, 0, 1, 2, 0, 0]
#      P2 = 1/2*C11 + C13 + 1/2*C33 + 2*C44
#      P3 = 1/6*C111 + 1/2*C113 + 1/2*C133 + 2*C144 + 1/6*C333
#           + 2*C344
#      P4 = 1/24*C1111 + 1/6*C1113 + 1/4*C1133 + C1144 + 1/6*C1333
#           + 2*C1344 + 1/24*C3333 + C3344 + 2/3*C4444
# M16: d = [1, 1, 1, 2, 0, 0]
#      P2 = C11 + C12 + 2*C13 + 1/2*C33 + 2*C44
#      P3 = 2/3*C111 + C112 + C113 + C123 + C133 + 2*C144 + 2*C155
#           - 1/3*C222 + 1/6*C333 + 2*C344
#      P4 = 11/108*C1111 + 17/54*C1112 + 1/3*C1113 + 1/4*C1122
#           + 2/3*C1123 + 1/2*C1133 + C1144 + C1155 - 1/9*C1166
#           + 1/3*C1223 + 1/2*C1233 + C1244 + C1255 + 1/3*C1333
#           + 2*C1344 + 2*C1355 + 1/24*C3333 + C3344 + 2/3*C4444
# M17: d = [1, 1, 1, 0, 0, 2]
#      P2 = 2*C11 + 2*C13 + 1/2*C33
#      P3 = 2/3*C111 + 2*C113 + C133 + 2/3*C222 + 1/6*C333
#      P4 = 14/27*C1111 + 4/27*C1112 + 4/3*C1113 + 2/3*C1123
#           + C1133 + 8/9*C1166 - 2/3*C1223 + 1/3*C1333 + 1/24*C3333
# M18: d = [1, 0, 1, 0, 0, 2]
#      P2 = 3/2*C11 - C12 + C13 + 1/2*C33
#      P3 = -5/6*C111 - 1/2*C112 + 3/2*C113 - C123 + 1/2*C133
#           + 3/2*C222 + 1/6*C333
#      P4 = 5/72*C1111 - 5/18*C1112 + 2/3*C1113 + 1/4*C1122
#           + C1123 + 3/4*C1133 + 4/3*C1166 - 3/2*C1223
#           - 1/2*C1233 + 1/6*C1333 + 1/24*C3333
# M19: d = [0, 1, -1, 0, 0, 0]
#      P2 = 1/2*C11 - C13 + 1/2*C33
#      P3 = -1/2*C113 + 1/2*C133 + 1/6*C222 - 1/6*C333
#      P4 = 5/216*C1111 + 1/54*C1112 - 1/6*C1113 - 1/6*C1123
#           + 1/4*C1133 + 1/9*C1166 + 1/6*C1223 - 1/6*C1333
#           + 1/24*C3333
# M20: d = [0, 1, 0, 0, 0, 0]
#      P2 = 1/2*C11
#      P3 = 1/6*C222
#      P4 = 5/216*C1111 + 1/54*C1112 + 1/9*C1166
#
# M04 and M06 have the same P2 and P4 but different P3, one classic redundant pair. Within the
# curated M01..M20 core, M01--M19 span only 18 of the 19 fourth-order constants after the
# physically correct dual projection; M20 supplies the missing direction. M21--M23 are three
# extra pure-normal directions, added to the default set so the pure-normal 4th-order block is
# fully determined by shear-free modes (see their comment below). The set is therefore
# deliberately over-determined; selftest_hoec() checks that every order stays full rank.
# hexagonal: 23 single-parameter modes (M01..M20 curated core + M21..M23 pure-normal).
MODES_HEX = {
    'M01': (1, 0, 0, 0, 0, 0),    # a-axis uniaxial
    'M02': (0, 0, 1, 0, 0, 0),    # c-axis uniaxial
    'M03': (1, 1, 0, 0, 0, 0),    # basal equibiaxial
    'M04': (1, -1, 0, 0, 0, 0),   # basal orthorhombic
    'M05': (0, 0, 0, 2, 0, 0),    # prism shear (C44)
    'M06': (0, 0, 0, 0, 0, 2),    # basal shear (C66)
    'M07': (1, 0, 1, 0, 0, 0),    # a+c biaxial
    'M08': (1, 1, 1, 0, 0, 0),    # volumetric-like
    'M09': (1, 0, 0, 2, 0, 0),
    'M10': (0, 0, 1, 2, 0, 0),
    'M11': (0, 0, 1, 0, 0, 2),
    'M12': (0, 0, 0, 2, 0, 2),
    'M13': (0, 0, 0, 2, 2, 2),
    'M14': (1, 0, 0, 2, 2, 0),
    'M15': (1, 0, 1, 2, 0, 0),
    'M16': (1, 1, 1, 2, 0, 0),
    'M17': (1, 1, 1, 0, 0, 2),
    'M18': (1, 0, 1, 0, 0, 2),
    'M19': (0, 1, -1, 0, 0, 0),
    'M20': (0, 1, 0, 0, 0, 0),    # second basal uniaxial direction; closes FOEC rank
    # Three extra pure-normal directions, now part of the default set. See
    # hcp_static_vs_relax_mode_zh.md sec.4: M01..M20's 8 pure-normal modes pin the 10 pure-normal
    # 4th-order constants only to rank 7/10, leaving C1113/C1333 in the underdetermined subspace,
    # so the fit had to lean on the shear modes (M17/M18) whose clamped-ion shear pollution leaks
    # into C1113/C1333. These raise the pure-normal C4 rank to 10/10 (verified: (1,0,2)->8,
    # (0,1,2)->9, (2,1,1)->10), fixing C1113/C1333 from shear-free modes alone. Pure normal keeps
    # the two atoms on their inversion centre (Dz==0.5) -- no FCC collapse, clamped-ion is exact.
    'M21': (1, 0, 2, 0, 0, 0),    # a1 + 2c   (pure-normal C4 rank 7 -> 8)
    'M22': (0, 1, 2, 0, 0, 0),    # a2 + 2c   (8 -> 9)
    'M23': (2, 1, 1, 0, 0, 0),    # anisotropic 1-2-3 mix, breaks eta1=eta2 (9 -> 10)
}
# Small-shear HCP option (opt-in, -small_shear). NOT a separate mode table.
#
# A metastable HCP branch survives a *static* (clamped-ion) scan at any strain, but once the
# ions are allowed to relax the basal shear opens an internal downhill channel and the cell
# shuffles into an FCC stacking. hcp_critical_window_and_runtime.md measured, per mode, the
# |xi| at which Au leaves the HCP branch: the shear-carrying modes transform first (safe
# window as small as +-0.03 for M15/M16, vs the uniaxial modes which never transform).
#
# small-shear lets a *relaxed-ion* run keep the shear directions without transforming: it
# multiplies the engineering shear entries (d4,d5,d6) of the selected core modes by
# SMALL_SHEAR_SCALE IN PLACE (get_hoec_modes / scale_mode_shear), so e.g.
# M09 [1 0 0 2 0 0] -> [1 0 0 1 0 0]: the tensor shear eps = xi*d/2 becomes HALF the
# accompanying tensor normal strain and is halved at any given xi. Modes without a shear entry
# are untouched, so the *default* target ("every shear-carrying mode") leaves the pure-normal
# modes exactly as they were. The mode NAMES are unchanged -- there is no M##s and no extra
# table -- so nothing is renamed and mymetal.post.hoec_energy must read each mode's direction
# from the run's manifest (it does), not from this global table.
#
# Scaling the direction is not enough on its own: get_mode_strain_lists' equal-severity scaling
# would just widen a now-milder mode's window until the shear is back where it started, and it
# does nothing at all for a pure-shear mode (its severity IS its shear). So the generator also
# caps each scaled mode's window at SMALL_SHEAR_SHEAR_CAP (get_shear_window). That cap also
# bounds the mode's normal-strain range, which is intended and harmless: these modes exist to
# carry the shear directions, while the normal constants come from the full-window pure-normal
# modes (M01, M02, ...).
SMALL_SHEAR_SCALE = 0.5

# Default cap on a scaled mode's largest tensor shear at its window edge. 0.03 keeps every Au
# shear mode inside the tightest safe window in hcp_critical_window_and_runtime.md
# (M15/M16, +-0.03); raise it for a stiffer metal, lower it if a relaxed run still transforms.
SMALL_SHEAR_SHEAR_CAP = 0.03

# MODES is the resolution fallback used by the post-processing when a manifest predates the
# per-mode ``direction`` field. small-shear adds no names (it rescales core modes in place), so
# nothing extra is folded in here. Generation selects/scales through get_hoec_modes.
MODES = {'cubic': MODES_CUBIC, 'hex': MODES_HEX}


def scale_mode_shear(d_dir: tuple = None, factor: float = None) -> tuple:
    """Multiply a mode's engineering shear entries (d4, d5, d6) by ``factor``; normals unchanged.

    This is the whole small-shear transform: ``M09 (1,0,0,2,0,0) -> (1,0,0,1,0,0)`` at
    ``factor=0.5``. A mode with no shear entry is returned unchanged.

    Args:
        d_dir (tuple): Length-6 engineering-Voigt mode direction.
        factor (float): Multiplier for the shear entries.

    Returns:
        tuple: The rescaled direction.
    """
    return tuple(v * factor if i >= 3 else v for i, v in enumerate(d_dir))


def get_shear_modes(symmetry: str = None) -> list:
    """Names of the core modes that carry an engineering shear entry (d4/d5/d6 != 0).

    This is the default target of ``-small_shear``: rescaling every shear-carrying mode leaves
    the pure-normal modes exactly as they were.

    Args:
        symmetry (str): 'cubic' or 'hex'.

    Returns:
        list: Mode names in the core-table order.

    Raises:
        ValueError: If ``symmetry`` is not a supported point group.
    """
    if symmetry not in MODES:
        raise ValueError("symmetry must be 'cubic' or 'hex', got %r" % symmetry)
    base = MODES_CUBIC if symmetry == 'cubic' else MODES_HEX
    return [name for name, d in base.items() if any(d[3:])]


def get_hoec_modes(symmetry: str = None, small_shear: bool = False,
                   shear_scale: float = SMALL_SHEAR_SCALE,
                   lsmall_shear_modes: list = None) -> dict:
    """Ordered mode table to generate for a symmetry, with the opt-in small-shear modifier.

    The base is the curated core set (cubic A..K, hex M01..M23). ``small_shear`` (hex only)
    rescales the engineering shear entries of the target modes by ``shear_scale`` IN PLACE,
    keeping their names -- so a relaxed-ion run gets the reduced-shear directions without adding
    or renaming any mode. The target defaults to every shear-carrying mode
    (:func:`get_shear_modes`); pass ``lsmall_shear_modes`` to rescale only some (the rest keep
    their original shear). With no modifier the returned dict is exactly the core set.

    Args:
        symmetry (str): 'cubic' or 'hex'.
        small_shear (bool): Rescale the target modes' shear entries (hex only).
        shear_scale (float): Shear multiplier for the small-shear target modes.
        lsmall_shear_modes (list): Modes to rescale; defaults to every shear-carrying mode.

    Returns:
        dict: Ordered mode name -> length-6 engineering-Voigt direction.

    Raises:
        ValueError: If ``symmetry`` is unsupported, ``shear_scale`` is not finite and positive,
            or a requested mode name is unknown.
    """
    if symmetry not in MODES:
        raise ValueError("symmetry must be 'cubic' or 'hex', got %r" % symmetry)
    dict_out = dict(MODES_CUBIC if symmetry == 'cubic' else MODES_HEX)
    if small_shear and symmetry == 'hex':
        if not np.isfinite(shear_scale) or shear_scale <= 0:
            raise ValueError("shear_scale must be finite and positive, got %r" % shear_scale)
        ltarget = (lsmall_shear_modes if lsmall_shear_modes is not None
                   else get_shear_modes(symmetry))
        lbad = [n for n in ltarget if n not in dict_out]
        if lbad:
            raise ValueError("unknown small-shear mode(s): %s" % ", ".join(lbad))
        for name in ltarget:
            dict_out[name] = scale_mode_shear(dict_out[name], shear_scale)
    return dict_out


# ---------------------------------------------------------------------------
# model
# ---------------------------------------------------------------------------
class HOECModel:
    """Reduced 2nd/3rd/4th-order elastic model for one crystal symmetry.

    Attributes:
        symmetry (str): 'cubic' or 'hex'.
        orders (dict): Per-order reduction data (monomials, R, names, pivots).
    """

    def __init__(self, symmetry: str = None):
        """Build the reduced model for all three orders.

        Args:
            symmetry (str): 'cubic' or 'hex'.

        Raises:
            ValueError: If ``symmetry`` is not a supported point group.
        """
        if symmetry not in GENERATORS:
            raise ValueError("symmetry must be 'cubic' or 'hex', got %r" % symmetry)
        self.symmetry = symmetry
        self.orders = {}
        for order in (2, 3, 4):
            lprefer = PREFER.get(symmetry, {}).get(order)
            lmono, R, lname, lpivot = reduce_symmetry(symmetry, order, lprefer=lprefer)
            # Rf: float view of the Fraction matrix, so P_coeffs stays a plain matmul
            self.orders[order] = dict(lmono=lmono, R=R, Rf=np.array(R, dtype=float),
                                      lname=lname, lpivot=lpivot)

    def names(self, order: int = None) -> list:
        """Voigt names of the independent order-n constants.

        Args:
            order (int): 2, 3 or 4.

        Returns:
            list: Voigt names, e.g. ['11', '12', '44'].
        """
        return list(self.orders[order]['lname'])

    def P_coeffs(self, order: int = None, d_dir: tuple = None) -> np.ndarray:
        """Coefficients over the independent order-n constants that give P_n.

        Args:
            order (int): 2, 3 or 4.
            d_dir (tuple): Length-6 engineering-Voigt mode direction.

        Returns:
            np.ndarray: Row of coefficients, aligned with :meth:`names`.
        """
        o = self.orders[order]
        mono = mono_eval_vec(np.array(d_dir, float), o['lmono'])
        return o['Rf'].T @ mono

    def system(self, order: int = None, lmode: list = None) -> np.ndarray:
        """Coefficient matrix A with ``P_order_vector = A @ constants``.

        Args:
            order (int): 2, 3 or 4.
            lmode (list): Mode directions (each length-6).

        Returns:
            np.ndarray: Matrix of shape (nmode, nconst).
        """
        return np.array([self.P_coeffs(order, d) for d in lmode])


_DICT_CACHE = {}


def get_model(symmetry: str = None) -> HOECModel:
    """Cached :class:`HOECModel` (the symmetry reduction is the expensive part).

    Args:
        symmetry (str): 'cubic' or 'hex'.

    Returns:
        HOECModel: The reduced model.
    """
    if symmetry not in _DICT_CACHE:
        _DICT_CACHE[symmetry] = HOECModel(symmetry)
    return _DICT_CACHE[symmetry]


# ---------------------------------------------------------------------------
# deformation
# ---------------------------------------------------------------------------
def get_strain_list(emax: float = 0.12, de: float = 0.01) -> list:
    """Symmetric list of xi in [-emax, emax] including 0.0, sorted ascending.

    Args:
        emax (float): Maximum |xi|.
        de (float): Step in xi.

    Returns:
        list: The xi values.

    Raises:
        ValueError: If ``emax`` or ``de`` is not positive.
    """
    if not (emax > 0 and de > 0):
        raise ValueError("emax and de must be positive, got emax=%r de=%r" % (emax, de))
    n = int(round(emax / de))
    return [round(k * de, 10) for k in range(-n, n + 1)]


def get_deformation_gradient(d_dir: tuple = None, xi: float = None) -> np.ndarray:
    """Symmetric F whose Green-Lagrange strain is eta = xi * d_dir.

    The engineering-Voigt direction is converted to the tensor strain (shears carry a
    1/2), then F = (I + 2 eta)^(1/2) is taken as the symmetric square root, so the mode
    carries no spurious rigid-body rotation.

    Args:
        d_dir (tuple): Length-6 engineering-Voigt mode direction.
        xi (float): Strain amplitude.

    Returns:
        np.ndarray: The symmetric 3x3 deformation gradient.

    Raises:
        ValueError: If I + 2 eta is not positive definite (non-physical deformation).
    """
    d1, d2, d3, d4, d5, d6 = d_dir
    eta = np.array([
        [xi * d1,       xi * d6 / 2.0, xi * d5 / 2.0],
        [xi * d6 / 2.0, xi * d2,       xi * d4 / 2.0],
        [xi * d5 / 2.0, xi * d4 / 2.0, xi * d3]])
    C = np.eye(3) + 2.0 * eta                       # C = F^T F for symmetric F
    w, V = np.linalg.eigh(C)
    if np.any(w <= 0):
        raise ValueError("non-physical deformation (F^T F not positive definite) "
                         "at xi=%.4f, d=%s" % (xi, d_dir))
    return (V * np.sqrt(w)) @ V.T                   # F = C^(1/2), symmetric


def get_mode_severity(d_dir: tuple = None, xi: float = None) -> float:
    """How far a mode has actually deformed the crystal at amplitude ``xi``.

    xi is a mode *amplitude*, not a strain: the Voigt directions carry shear entries of 2
    and up to three nonzero components, so one xi is a wildly different deformation from
    mode to mode. At xi=0.12 the cubic modes range from -12.8% volume (A) to -33.7% (K),
    and mode J reaches a 35% principal strain -- far outside the radius where a 4th-order
    Taylor expansion of the energy means anything.

    Severity is the worse of the two ways a cell can be pushed too far::

        max( max_i |lambda_i(F) - 1| ,  |det F - 1| )

    i.e. the largest principal stretch deviation and the volume change. Both are needed:
    the principal stretch alone misses the hydrostatic mode K (12.8% principal strain but
    -33.7% volume), and the volume change alone misses the pure shears (F, G).

    For the uniaxial reference mode A the two coincide, so ``get_mode_severity(A, xi)``
    is just A's own strain and can be used directly as the target for every other mode
    (see :func:`get_mode_window`).

    Args:
        d_dir (tuple): Length-6 engineering-Voigt mode direction.
        xi (float): Strain amplitude; both signs are tested and the worse one returned,
            since compression is always the stiffer, more nonlinear side.

    Returns:
        float: The severity, as a dimensionless strain-like number.
    """
    sev = 0.0
    for s in (abs(xi), -abs(xi)):
        F = get_deformation_gradient(d_dir, s)
        w = np.linalg.eigvalsh(F)
        sev = max(sev, np.max(np.abs(w - 1.0)), abs(np.linalg.det(F) - 1.0))
    return float(sev)


def get_mode_window(d_dir: tuple = None, target: float = None,
                    hi: float = 0.5, tol: float = 1e-6) -> float:
    """Largest |xi| at which a mode's severity still stays within ``target``.

    Severity increases monotonically with |xi|, so a bisection is exact.

    Args:
        d_dir (tuple): Length-6 engineering-Voigt mode direction.
        target (float): Severity budget, normally ``get_mode_severity(d_A, emax)``.
        hi (float): Upper bracket for the bisection.
        tol (float): Bracket width at which to stop.

    Returns:
        float: The mode's own maximum |xi|.

    Raises:
        ValueError: If ``target`` is not positive.
    """
    if not target > 0:
        raise ValueError("target severity must be positive, got %r" % target)
    lo = 0.0
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        try:
            ok = get_mode_severity(d_dir, mid) <= target
        except ValueError:                  # F^T F not positive definite this far out
            ok = False
        if ok:
            lo = mid
        else:
            hi = mid
    return lo


def get_shear_window(d_dir: tuple = None, shear_cap: float = None) -> float:
    """Largest |xi| at which a mode's tensor shear strain still stays within ``shear_cap``.

    The tensor shear grows linearly with |xi| (the engineering shear entry d4/d5/d6 carries the
    tensor value xi*d/2), so the bound is exact: ``|xi| <= shear_cap / (max|d4,d5,d6|/2)``. A mode
    with no shear entry is unbounded and returns +inf. This is the extra cap the small-shear
    modes take on top of the severity budget, because equal-severity scaling alone does not
    bound a pure-shear mode's shear (its severity IS its shear, so an equal-severity window keeps
    the shear unchanged).

    Args:
        d_dir (tuple): Length-6 engineering-Voigt mode direction.
        shear_cap (float): Maximum allowed tensor shear strain.

    Returns:
        float: The mode's own maximum |xi|, or ``float('inf')`` when the mode carries no shear.

    Raises:
        ValueError: If ``shear_cap`` is not finite and positive.
    """
    if not np.isfinite(shear_cap) or shear_cap <= 0:
        raise ValueError("shear_cap must be finite and positive, got %r" % shear_cap)
    shear_per_xi = max(abs(d_dir[3]), abs(d_dir[4]), abs(d_dir[5])) / 2.0
    return float('inf') if shear_per_xi == 0 else shear_cap / shear_per_xi


def get_mode_strain_lists(symmetry: str = None, emax: float = 0.12, de: float = 0.01,
                          scale_window: bool = True, dict_modes: dict = None,
                          shear_cap: float = None, lcap_modes: list = None) -> dict:
    """Per-mode xi lists, each mode confined to an equally severe deformation.

    The uniaxial mode A sets the budget: it keeps ``emax`` and ``de`` exactly as given, and
    every other mode is shrunk to the |xi| at which it deforms the crystal as hard as A does
    at ``emax`` (:func:`get_mode_severity`). The step is shrunk by the same factor, so every
    mode keeps the *same number of points* -- the fits stay equally conditioned, and the
    heavily-scaled modes automatically get the fine xi increment Wang-Li call for.

    With ``scale_window=False`` every mode gets the same ``[-emax, emax]``, which is the
    literal reading of Wang-Li and the original behaviour of this workflow.

    ``shear_cap`` adds a second bound on the modes named in ``lcap_modes`` (the small-shear
    modes): their window is the smaller of the severity window and the :func:`get_shear_window`
    shear cap, so a relaxed run's basal shear cannot reach the value that tips HCP into FCC.
    Modes not in ``lcap_modes`` are never capped. Ignored when either is unset.

    Args:
        symmetry (str): 'cubic' or 'hex'.
        emax (float): Maximum |xi| for the reference mode.
        de (float): xi step for the reference mode.
        scale_window (bool): Scale each mode's window by its severity.
        dict_modes (dict): Mode subset to build (name -> direction); defaults to the core set
            for this symmetry. Pass :func:`get_hoec_modes` output to include any modifiers.
        shear_cap (float): Max tensor shear for the ``lcap_modes``; ``None`` disables the cap.
        lcap_modes (list): Mode names to apply ``shear_cap`` to (the small-shear modes).

    Returns:
        dict: Mode name -> dict with ``xi`` (list), ``emax``, ``de`` and ``scale``
        (the mode's window divided by ``emax``; 1.0 for every mode when not scaling).

    Raises:
        ValueError: If ``symmetry`` is unsupported or ``shear_cap`` is not finite and positive.
    """
    if symmetry not in MODES:
        raise ValueError("symmetry must be 'cubic' or 'hex', got %r" % symmetry)
    if shear_cap is not None and (not np.isfinite(shear_cap) or shear_cap <= 0):
        raise ValueError("shear_cap must be finite and positive, got %r" % shear_cap)
    dict_modes = dict_modes if dict_modes is not None else get_hoec_modes(symmetry)
    lxi_ref = get_strain_list(emax, de)
    npt = (len(lxi_ref) - 1) // 2                  # points on each side of xi=0
    scap = set(lcap_modes or [])
    # severity budget: mode A is uniaxial, so its principal stretch and its volume change
    # coincide -- its own severity at emax is the natural target for every other mode
    target = get_mode_severity(MODES[symmetry]['A' if symmetry == 'cubic' else 'M01'], emax)

    dict_out = {}
    for name, d in dict_modes.items():
        # scale_window shrinks every mode to equal severity; without it all modes share emax
        xi_max = get_mode_window(d, target) if scale_window else emax
        if shear_cap is not None and name in scap:
            xi_max = min(xi_max, get_shear_window(d, shear_cap))   # keep the shear sub-critical
        scale = min(xi_max / emax, 1.0)           # never widen past the requested emax
        de_m = de * scale
        dict_out[name] = dict(xi=[round(k * de_m, 10) for k in range(-npt, npt + 1)],
                              emax=emax * scale, de=de_m, scale=scale)
    return dict_out


def check_symmetry(symmetry: str = 'auto', cell: np.ndarray = None) -> str:
    """Validate an explicit symmetry, or auto-detect it from the cell.

    Args:
        symmetry (str): 'auto', 'cubic' or 'hex'.
        cell (np.ndarray): 3x3 cell with lattice vectors as rows.

    Returns:
        str: 'cubic' or 'hex'.

    Raises:
        ValueError: If ``symmetry`` is invalid, or if 'auto' cannot classify the cell.
    """
    if symmetry in ('cubic', 'hex'):
        return symmetry
    if symmetry != 'auto':
        raise ValueError("symmetry must be auto|cubic|hex, got %r" % symmetry)
    cell = np.asarray(cell, float)
    a, b, c = np.linalg.norm(cell, axis=1)
    cosang = lambda u, v: np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
    gamma = np.degrees(np.arccos(cosang(cell[0], cell[1])))
    if abs(a - b) < 1e-3 * a and abs(a - c) < 1e-3 * a and abs(gamma - 90) < 1.0:
        return 'cubic'
    if abs(a - b) < 1e-3 * a and abs(gamma - 120) < 1.0:
        return 'hex'
    raise ValueError("could not auto-detect symmetry (a,b,c=%.4f,%.4f,%.4f gamma=%.1f); "
                     "pass symmetry='cubic' or 'hex' explicitly" % (a, b, c, gamma))


# ---------------------------------------------------------------------------
# self-test
# ---------------------------------------------------------------------------
def _fmt(lcoeff: np.ndarray = None, lname: list = None) -> str:
    """Format a P-coefficient row as a human-readable linear combination."""
    lterm = []
    for c, nm in zip(lcoeff, lname):
        f = Fraction(c).limit_denominator(4620)
        if f != 0:
            lterm.append('%s*C%s' % (str(f), nm))
    return ' + '.join(lterm) if lterm else '0'


def _get_individually_removable_modes(model: HOECModel = None,
                                      dict_modes: dict = None,
                                      orders: tuple = (2, 3, 4)) -> list:
    """Find modes whose individual removal preserves full rank at every order."""
    lremovable = []
    for name_drop in dict_modes:
        ldirections = [direction for name, direction in dict_modes.items()
                       if name != name_drop]
        if all(np.linalg.matrix_rank(model.system(order, ldirections), tol=1e-7)
               == len(model.names(order)) for order in orders):
            lremovable.append(name_drop)
    return lremovable


def selftest_hoec() -> bool:
    """Verify cubic formulas plus hex invariance, identities and mode-set rank.

    Returns:
        bool: True when the cubic P coefficients reproduce Table I, the hex
        coefficients are rotationally invariant, analytic SOEC identities hold,
        and every mode-set/order combination has full column rank.
    """
    ### cubic: compare against Wang-Li Table I for unambiguous modes ---------------
    dict_ref = {
        'A': ('1/2*C11', '1/6*C111', '1/24*C1111'),
        'F': ('6*C44', '8*C456', '2*C4444 + 12*C4455'),
        'G': ('2*C44', '0', '2/3*C4444'),
        'K': ('3/2*C11 + 3*C12', '1/2*C111 + 3*C112 + 1*C123',
              '1/8*C1111 + 1*C1112 + 3/4*C1122 + 3/2*C1123'),
    }
    print('================ ✅ cubic (FCC) — reduced model & Table I check')
    mc = get_model('cubic')
    for o in (2, 3, 4):
        print('  order %d (%d): %s' % (o, len(mc.names(o)), ['C' + n for n in mc.names(o)]))
    ok = True
    for name, d in MODES_CUBIC.items():
        got = tuple(_fmt(mc.P_coeffs(od, d), mc.names(od)) for od in (2, 3, 4))
        flag = ''
        if name in dict_ref:
            match = all(g == r for g, r in zip(got, dict_ref[name]))
            ok = ok and match
            flag = '  ✅' if match else '  ❌ EXPECTED ' + str(dict_ref[name])
        print('  mode %s d=%s' % (name, d))
        print('      P2=%s | P3=%s | P4=%s%s' % (got[0], got[1], got[2], flag))
    for o in (2, 3, 4):
        A = mc.system(o, list(MODES_CUBIC.values()))
        r = np.linalg.matrix_rank(A, tol=1e-7)
        full = (r == len(mc.names(o)))
        ok = ok and full
        print('  order %d rank %d/%d %s' % (o, r, len(mc.names(o)), '✅' if full else '❌'))
    lcubic_removable = _get_individually_removable_modes(mc, MODES_CUBIC)
    cubic_minimal = not lcubic_removable
    ok = ok and cubic_minimal
    print('  orders 2--4 jointly removable modes: %s %s'
          % (', '.join(lcubic_removable) if lcubic_removable else 'none',
             '✅' if cubic_minimal else '❌ EXPECTED none'))

    ### hexagonal: invariance, analytic SOEC identities and rank -------------------
    print('\n================ ✅ hexagonal (HCP) — invariance, identities & rank')
    mh = get_model('hex')
    for o in (2, 3, 4):
        print('  order %d (%d): %s' % (o, len(mh.names(o)), ['C' + n for n in mh.names(o)]))

    rng = np.random.default_rng(20260714)
    lprobe = rng.standard_normal((30, 6))
    lgroup = close_group(GENERATORS['hex'])
    for o in (2, 3, 4):
        err = max(np.max(np.abs(mh.P_coeffs(o, voigt_strain_M(R3) @ d)
                                - mh.P_coeffs(o, d)))
                  for d in lprobe for R3 in lgroup)
        invariant = err < 1e-10
        ok = ok and invariant
        print('  order %d rotational-invariance max error %.3e %s'
              % (o, err, '✅' if invariant else '❌'))

    expected_m04_m06 = np.array([1.0, -1.0, 0.0, 0.0, 0.0])
    row_m04 = mh.P_coeffs(2, MODES_HEX['M04'])
    row_m06 = mh.P_coeffs(2, MODES_HEX['M06'])
    identity = (np.allclose(row_m04, expected_m04_m06, atol=1e-12)
                and np.allclose(row_m06, expected_m04_m06, atol=1e-12))
    ok = ok and identity
    print('  M04/M06 P2 = C11-C12: %s / %s %s'
          % (_fmt(row_m04, mh.names(2)), _fmt(row_m06, mh.names(2)),
             '✅' if identity else '❌'))

    for o in (2, 3, 4):
        A = mh.system(o, list(MODES_HEX.values()))
        r = np.linalg.matrix_rank(A, tol=1e-7)
        full = (r == len(mh.names(o)))
        ok = ok and full
        over = A.shape[0] - r
        print('  order %d system %dx%d rank %d/%d cond=%.1f row-redundancy=%d %s'
              % (o, A.shape[0], A.shape[1], r, len(mh.names(o)), np.linalg.cond(A),
                 over, '✅' if full else '❌'))

    # M04/M06 identity: same P2 and P4, different P3 -- a symmetry fact independent of the set.
    m04_m06 = (np.allclose(mh.P_coeffs(2, MODES_HEX['M04']),
                           mh.P_coeffs(2, MODES_HEX['M06']), atol=1e-12)
               and not np.allclose(mh.P_coeffs(3, MODES_HEX['M04']),
                                   mh.P_coeffs(3, MODES_HEX['M06']), atol=1e-12)
               and np.allclose(mh.P_coeffs(4, MODES_HEX['M04']),
                               mh.P_coeffs(4, MODES_HEX['M06']), atol=1e-12))
    ok = ok and m04_m06
    print('  M04/M06 identity (P2,P4 equal, P3 differ): %s' % ('✅' if m04_m06 else '❌'))
    # the set is deliberately over-determined (M21..M23 add pure-normal rows), so several modes
    # are now individually removable; the redundancy is reported, not asserted to a fixed set.
    lremovable = _get_individually_removable_modes(mh, MODES_HEX)
    print('  orders 2--4 individually removable (over-determined, informational): %s'
          % (', '.join(lremovable) if lremovable else 'none'))

    ### hex small-shear: in-place shear rescale, shear cap, rank preserved --------
    print('\n================ ✅ hexagonal small-shear (-small_shear, in-place rescale)')
    dict_small = get_hoec_modes('hex', small_shear=True)      # default: every shear mode, x0.5
    target_hex = get_mode_severity(MODES_HEX['M01'], 0.12)
    small_ok = True
    lshear = get_shear_modes('hex')
    for name in lshear:
        d0, d = MODES_HEX[name], dict_small[name]
        # normals untouched; each shear entry is exactly halved
        rel_ok = (d[:3] == d0[:3]
                  and all(abs(d[i] - d0[i] * SMALL_SHEAR_SCALE) < 1e-12 for i in range(3, 6))
                  and any(d0[3:]))
        xi_max = min(get_mode_window(d, target_hex), get_shear_window(d, SMALL_SHEAR_SHEAR_CAP))
        shear_edge = xi_max * max(abs(d[3]), abs(d[4]), abs(d[5])) / 2.0
        cap_ok = shear_edge <= SMALL_SHEAR_SHEAR_CAP + 1e-9
        small_ok = small_ok and rel_ok and cap_ok
        print('  %-4s %-16s -> %-18s xi_max=%.4f shear_edge=%.4f %s'
              % (name, str(d0), str(d), xi_max, shear_edge, '✅' if rel_ok and cap_ok else '❌'))
    # rescaling the shear in place must not lower any order's rank (the run stays solvable)
    for o in (2, 3, 4):
        r = np.linalg.matrix_rank(mh.system(o, list(dict_small.values())), tol=1e-7)
        full = (r == len(mh.names(o)))
        small_ok = small_ok and full
        print('  order %d rescaled-set rank %d/%d %s' % (o, r, len(mh.names(o)), '✅' if full else '❌'))
    # non-shear modes must be left exactly as they were (untouched by the default target)
    untouched = all(dict_small[n] == MODES_HEX[n] for n in MODES_HEX if n not in lshear)
    small_ok = small_ok and untouched
    print('  pure-normal modes untouched by default target: %s' % ('✅' if untouched else '❌'))
    ok = ok and small_ok

    print('\n%s' % ('🎉 self-test PASSED' if ok else '❌ self-test FAILED'))
    return ok


if __name__ == '__main__':
    raise SystemExit(0 if selftest_hoec() else 1)
