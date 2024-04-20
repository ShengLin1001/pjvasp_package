from ase import Atoms
from mymetal.build.film.findhetero import build_supercells

def compare_atoms(atoms1: Atoms = None, atoms2: Atoms = None, tolerance: float=1e-5, focus_z: bool = True) -> list:
    """
    Compare cell length, angle, cartesian position, scaled position.\n
    Return a list of bool
        >> `[True, Ture, True, True]`
    """
    flag = []
    if focus_z:
        # Compare cell lengths
        cell_lengths1 = atoms1.cell.cellpar()[:3]
        cell_lengths2 = atoms2.cell.cellpar()[:3]
        # Compare cell angles
        cell_angles1 = atoms1.cell.cellpar()[3:]
        cell_angles2 = atoms2.cell.cellpar()[3:]
        # Compare atom positions
        positions1 = atoms1.get_positions()
        positions2 = atoms2.get_positions()
        # Compare atom scaled positions
        scaled_positions1 = atoms1.get_scaled_positions()
        scaled_positions2 = atoms2.get_scaled_positions()
    else:
        # Compare cell lengths
        cell_lengths1 = atoms1.cell.cellpar()[:2]
        cell_lengths2 = atoms2.cell.cellpar()[:2]
        # Compare cell angles
        cell_angles1 = atoms1.cell.cellpar()[5]
        cell_angles2 = atoms2.cell.cellpar()[5]
        # Compare atom positions
        positions1 = atoms1.get_positions()[:,:2]
        positions2 = atoms2.get_positions()[:,:2]
        # Compare atom scaled positions
        scaled_positions1 = atoms1.get_scaled_positions()[:,:2]
        scaled_positions2 = atoms2.get_scaled_positions()[:,:2]

    if not all(abs(l1 - l2) < tolerance for l1, l2 in zip(cell_lengths1, cell_lengths2)):
        flag.append(False)
    else:
        flag.append(True)
    
    if not all(abs(a1 - a2) < tolerance for a1, a2 in zip(cell_angles1, cell_angles2)):
        flag.append(False)
    else:
        flag.append(True)
    
    if not all(all(abs(p1 - p2) < tolerance for p1, p2 in zip(pos1, pos2)) for pos1, pos2 in zip(positions1, positions2)):
        flag.append(False)
    else:
        flag.append(True)

    if not all(all(abs(p1 - p2) < tolerance for p1, p2 in zip(pos1, pos2)) for pos1, pos2 in zip(scaled_positions1, scaled_positions2)):
        flag.append(False)
    else:
        flag.append(True)

    return flag

def cal_mismatch(bottom: Atoms = None, top: Atoms = None, hetero: Atoms = None) -> float:
    """
    For calculating the max mismatch.
    ------
    return a list\n
        >> `[a_mismatch, b_mismatch, gamma_mismatch]`
    """
    a_bottom, b_bottom, c_bottom, alpha_bottom, beta_bottom, gamma_bottom = bottom.cell.cellpar()
    a_top, b_top, c_top, alpha_top, beta_top, gamma_top = top.cell.cellpar()
    a_hetero, b_hetero, c_hetero, alpha_hetero, beta_hetero, gamma_hetero = hetero.cell.cellpar()

    mismatch = []
    for bot, to, het in zip([a_bottom, b_bottom, gamma_bottom],[a_top, b_top, gamma_top], [a_hetero, b_hetero, gamma_hetero]):
        mismatch.append(relative_diff(bot, to, het))
    return mismatch

def cal_stretch(bottom: Atoms = None, top: Atoms = None, hetero: Atoms = None) -> float:
    """
    For calculating the stretch of bottom and top layer.
    -------
    return two list:
        >> `[a_bot_stretch, b_bot_stretch, gamma_bot_stretch], [a_to_stretch, b_to_stretch, gamma_to_stretch]`
    """
    a_bottom, b_bottom, c_bottom, alpha_bottom, beta_bottom, gamma_bottom = bottom.cell.cellpar()
    a_top, b_top, c_top, alpha_top, beta_top, gamma_top = top.cell.cellpar()
    a_hetero, b_hetero, c_hetero, alpha_hetero, beta_hetero, gamma_hetero = hetero.cell.cellpar()

    bot_stretch = []
    for bot, het in zip([a_bottom, b_bottom, gamma_bottom],[a_hetero, b_hetero, gamma_hetero]):
        bot_stretch.append((het-bot)/bot)

    to_stretch = []
    for to, het in zip([a_top, b_top, gamma_top],[a_hetero, b_hetero, gamma_hetero]):
        to_stretch.append((het-to)/to)
    return bot_stretch, to_stretch

def relative_diff(bottom: float = None, top: float = None, hetero: float = None) -> float:
    """
    For calculate: `max(abs(bottom - hetero)/bottom, abs(top - hetero)/top)`.
    """
    max_mismatch = max(abs(bottom - hetero)/bottom, abs(top - hetero)/top) 
    return max_mismatch

def cal_atom_num(atoms: Atoms = None) -> int:
    return len(atoms.get_positions()[:,0])

def filter_results(bottom: float = None, top: float = None, results: list = None, 
                   length_mismatch: float = 0.05, angle_mismatch: float = 0.05, sorting_func = lambda x: cal_atom_num(x.stack)):
    """
    Filter results of heterostructures based on specified mismatches in length and angle.
    ------
    Args:
        bottom (float): Parameter related to the bottom layer of the heterostructure.
        top (float): Parameter related to the top layer of the heterostructure.
        results (list): List of heterostructure results to filter.
        length_mismatch (float): Maximum allowable relative difference in length.
        angle_mismatch (float): Maximum allowable angle difference in degrees.
        sorting_func (function): Function to determine the key for sorting results.

    Returns:
        list: Filtered list of heterostructure results.
    """
    filted_results = []
    for result in results:
        bot, to = build_supercells(bottom, top, result.M, result.N, result.angle, if_stack = False)
        mystack = build_supercells(bottom, top, result.M, result.N, result.angle, result._weight)
        if max(cal_mismatch(bot, to, mystack)[0:2]) < length_mismatch and cal_mismatch(bot, to, mystack)[2] < angle_mismatch:
            filted_results.append(result)
    filted_results.sort(key=sorting_func)
    return filted_results