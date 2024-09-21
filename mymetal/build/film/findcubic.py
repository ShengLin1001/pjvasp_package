from ase import Atoms
from math import sqrt
from mymetal.build.film.extrfilm import cal_area
import numpy as np
import warnings
def find_cubic( prim: Atoms = None, type: str = 'hcp') -> Atoms:
    """
    Find the xy plane cubic cell of the primitive atoms(hcp/fcc film).\n
    input [a,b,c,any, any, 90]
    """
    atoms = prim.copy()
    cell = atoms.get_cell()
    area = cal_area(atoms)
    #volume = atoms.get_volume()
    if type == 'hcp' or 'fcc':
        cell[0,0] = sqrt(area/sqrt(3.0))
        cell[1,1] = sqrt(3.0)*cell[0,0]
        cell[0,1] = 0
        cell[1,0] = 0
    else:
        print('Wrong type!')
    atoms.set_cell(cell)
    return atoms

# take from ase.build.supercells
def find_optimal_cell_shape(
    cell,
    target_size,
    target_shape,
    lower_limit=-2,
    upper_limit=2,
    verbose=False,
):
    """Obtain the optimal transformation matrix for a supercell of target size
    and shape.

    Returns the transformation matrix that produces a supercell
    corresponding to *target_size* unit cells with metric *cell* that
    most closely approximates the shape defined by *target_shape*.

    Updated with code from the `doped` defect simulation package
    (https://doped.readthedocs.io) to be rotationally invariant and
    allow transformation matrices with negative determinants, boosting
    performance.

    Please note that the function is based on the ASE implementation.
    https://wiki.fysik.dtu.dk/ase/ase/build/tools.html#ase.build.find_optimal_cell_shape
    https://wiki.fysik.dtu.dk/ase/tutorials/defects/defects.html#supercell-creation

    Parameters:

    cell: 2D array of floats
        Metric given as a (3x3 matrix) of the input structure.
    target_size: integer
        Size of desired supercell in number of unit cells.
    target_shape: str
        Desired supercell shape. Can be 'sc' for simple cubic or
        'fcc' for face-centered cubic.
    lower_limit: int
        Lower limit of search range.
    upper_limit: int
        Upper limit of search range.
    verbose: bool
        Set to True to obtain additional information regarding
        construction of transformation matrix.

    Returns:
        2D array of integers: Transformation matrix that produces the
        optimal supercell.
    """
    cell = np.asarray(cell)

    # Set up target metric
    if target_shape == 'sc':
        target_metric = np.eye(3)
    elif target_shape == 'fcc':
        target_metric = 0.5 * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]],
                                       dtype=float)
    if verbose:
        print("target metric (h_target):")
        print(target_metric)

    # Normalize cell metric to reduce computation time during looping
    norm = (target_size * abs(np.linalg.det(cell)) /
            np.linalg.det(target_metric)) ** (-1.0 / 3)
    norm_cell = norm * cell
    if verbose:
        print("normalization factor (Q): %g" % norm)

    # Approximate initial P matrix
    ideal_P = np.dot(target_metric, np.linalg.inv(norm_cell))
    if verbose:
        print("idealized transformation matrix:")
        print(ideal_P)
    starting_P = np.array(np.around(ideal_P, 0), dtype=int)
    if verbose:
        print("closest integer transformation matrix (P_0):")
        print(starting_P)

    # Prepare run.
    from itertools import product

    best_score = 1e6
    optimal_P = None
    for dP in product(range(lower_limit, upper_limit + 1), repeat=9):
        dP = np.array(dP, dtype=int).reshape(3, 3)
        P = starting_P + dP
        if int(np.around(np.linalg.det(P), 0)) != target_size:
            continue
        score = get_deviation_from_optimal_cell_shape(P @ cell, target_shape)
        if score < best_score:
            best_score = score
            optimal_P = P

    if optimal_P is None:
        print("Failed to find a transformation matrix.")
        return None

    if np.linalg.det(optimal_P) <= 0:
        optimal_P *= -1  # flip signs if negative determinant

    # Finalize.
    if verbose:
        print("smallest score (|Q P h_p - h_target|_2): %f" % best_score)
        print("optimal transformation matrix (P_opt):")
        print(optimal_P)
        print("supercell metric:")
        print(np.round(np.dot(optimal_P, cell), 4))
        print("determinant of optimal transformation matrix: %g" %
              np.linalg.det(optimal_P))
    return optimal_P

def get_deviation_from_optimal_cell_shape(cell, target_shape="sc", norm=None):
    r"""Calculate the deviation from the target cell shape.

    Calculates the deviation of the given cell metric from the ideal
    cell metric defining a certain shape. Specifically, the function
    evaluates the expression `\Delta = || Q \mathbf{h} -
    \mathbf{h}_{target}||_2`, where `\mathbf{h}` is the input
    metric (*cell*) and `Q` is a normalization factor (*norm*)
    while the target metric `\mathbf{h}_{target}` (via
    *target_shape*) represent simple cubic ('sc') or face-centered
    cubic ('fcc') cell shapes.

    Replaced with code from the `doped` defect simulation package
    (https://doped.readthedocs.io) to be rotationally invariant,
    boosting performance.

    Parameters:

    cell: 2D array of floats
        Metric given as a (3x3 matrix) of the input structure.
    target_shape: str
        Desired supercell shape. Can be 'sc' for simple cubic or
        'fcc' for face-centered cubic.
    norm: float
        Specify the normalization factor. This is useful to avoid
        recomputing the normalization factor when computing the
        deviation for a series of P matrices.

    Returns:
        float: Cell metric (0 is perfect score)

    .. deprecated:: 3.24.0
        `norm` is unused in ASE 3.24.0 and removed in ASE 3.25.0.

    """
    if norm is not None:
        warnings.warn(
            '`norm` is unused in ASE 3.24.0 and removed in ASE 3.25.0',
            FutureWarning,
        )

    cell_lengths = np.linalg.norm(cell, axis=1)
    eff_cubic_length = float(abs(np.linalg.det(cell)) ** (1 / 3))  # 'a_0'

    if target_shape == 'sc':
        target_length = eff_cubic_length

    elif target_shape == 'fcc':
        # FCC is characterised by 60 degree angles & lattice vectors = 2**(1/6)
        # times the eff cubic length:
        target_length = eff_cubic_length * 2 ** (1 / 6)

    inv_target_length = 1.0 / target_length

    # rms difference to eff cubic/FCC length:
    return np.sqrt(np.sum((cell_lengths * inv_target_length - 1.0) ** 2))

