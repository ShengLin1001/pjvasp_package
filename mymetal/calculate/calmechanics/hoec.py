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
    - check_symmetry: validate / auto-detect 'cubic' or 'hex' from a cell.
    - selftest_hoec: verify cubic against Wang-Li Table I and check mode-set rank.

Change log:
    - Written by J. P., 2026, for the Cu/Ag/Au FCC+HCP higher-order elastic-constant
      workflow (``hoec_core.py``).
    - Merged into mymetal by J. P. on 2026-07-10; the deformation gradient, the
      strain list and the symmetry detection moved in from the generator script.
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
        Pr += Tg
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
# hexagonal: 19 curated single-parameter modes, chosen so that the mode set keeps
# full column rank at every order (SOEC 5, TOEC 10, FOEC 19). Dropping or reordering
# them can make the 4th-order system singular -- selftest_hoec() checks the rank.
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
}
MODES = {'cubic': MODES_CUBIC, 'hex': MODES_HEX}


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


def selftest_hoec() -> bool:
    """Verify the reduced model: cubic vs Wang-Li Table I, and mode-set rank.

    Returns:
        bool: True when the cubic P coefficients reproduce Table I and every
        mode-set/order combination has full column rank.
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

    ### hexagonal: only a rank check, no published table to compare against --------
    print('\n================ ✅ hexagonal (HCP) — reduced model & rank check')
    mh = get_model('hex')
    for o in (2, 3, 4):
        print('  order %d (%d): %s' % (o, len(mh.names(o)), ['C' + n for n in mh.names(o)]))
    for o in (2, 3, 4):
        A = mh.system(o, list(MODES_HEX.values()))
        r = np.linalg.matrix_rank(A, tol=1e-7)
        full = (r == len(mh.names(o)))
        ok = ok and full
        print('  order %d rank %d/%d cond=%.0f %s' %
              (o, r, len(mh.names(o)), np.linalg.cond(A), '✅' if full else '❌'))
    print('\n%s' % ('🎉 self-test PASSED' if ok else '❌ self-test FAILED'))
    return ok


if __name__ == '__main__':
    raise SystemExit(0 if selftest_hoec() else 1)
