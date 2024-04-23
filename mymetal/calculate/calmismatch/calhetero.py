from ase import Atoms
from mymetal.build.film.findhetero import build_supercells
from numpy import array
import numpy as np
from mymetal.calculate.calmechanics.deformation import cal_deform_matrix
from mymetal.calculate.calmechanics.strain import cal_principal_and_shear_strain, cal_strain_matrix, cal_strain_matrix_root
from hetbuilder import Interface
from mymetal.universial.matrix.adjust import adjust_matrix

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


def cal_stretch_lattice(bottom: Atoms = None, top: Atoms = None, result: list = None, 
                   bot: Atoms = None, to: Atoms = None, stack: Atoms = None) -> float:
    """
    For calculating the stretch of bottom and top layer according to lattice constant.
    
    Returns:
        two list 
        >> `[a_bot_stretch, b_bot_stretch, gamma_bot_stretch], [a_to_stretch, b_to_stretch, gamma_to_stretch]`
    
    Raise:
        ValueError: The input [bottom, top, results] or [bot, to, stack] must not be None simultaneously.
    """
    if not (all(value != None for value in [bottom, top, result]) or all(value != None for value in [bot, to, stack])):
        raise ValueError('The input [bottom, top, results] or [bot, to, stack] must not be None simultaneously.')
    if bot == None and to == None and stack == None:
        bot, to = build_supercells(bottom, top, result.M, result.N, result.angle, if_stack = False)
        stack = build_supercells(bottom, top, result.M, result.N, result.angle, result._weight)
    a_bottom, b_bottom, c_bottom, alpha_bottom, beta_bottom, gamma_bottom = bot.cell.cellpar()
    a_top, b_top, c_top, alpha_top, beta_top, gamma_top = to.cell.cellpar()
    a_hetero, b_hetero, c_hetero, alpha_hetero, beta_hetero, gamma_hetero = stack.cell.cellpar()

    bot_stretch = []
    for bot, het in zip([a_bottom, b_bottom, gamma_bottom],[a_hetero, b_hetero, gamma_hetero]):
        bot_stretch.append((het-bot)/bot)

    to_stretch = []
    for to, het in zip([a_top, b_top, gamma_top],[a_hetero, b_hetero, gamma_hetero]):
        to_stretch.append((het-to)/to)
    return bot_stretch, to_stretch

def relative_diff(bottom: float = None, top: float = None, hetero: float = None) -> float:
    """For calculate: `max(abs(bottom - hetero)/bottom, abs(top - hetero)/top)`."""
    max_mismatch = max(abs(bottom - hetero)/bottom, abs(top - hetero)/top) 
    return max_mismatch

def cal_atom_num(atoms: Atoms = None) -> int:
    return len(atoms.get_positions()[:,0])

# def filter_results(bottom: float = None, top: float = None, results: list = None, 
#                    tag: str = 'percent',
#                    length_mismatch: float = 0.05, angle_mismatch: float = 0.05, 
#                    sorting_func = lambda x: cal_atom_num(x.stack)):
#     """Filter results of heterostructures based on specified mismatches in length and angle.
    
#     Args:
#         bottom (float): Parameter related to the bottom layer of the heterostructure.
#         top (float): Parameter related to the top layer of the heterostructure.
#         results (list): List of heterostructure results to filter.
#         length_mismatch (float): Maximum allowable relative difference in length.
#         angle_mismatch (float): Maximum allowable angle difference in degrees.
#         sorting_func (function): Function to determine the key for sorting results.

#     Returns:
#         list: Filtered list of heterostructure results.
#     """
#     filted_results = []
#     for result in results:
#         bot, to = build_supercells(bottom, top, result.M, result.N, result.angle, if_stack = False)
#         mystack = build_supercells(bottom, top, result.M, result.N, result.angle, result._weight)
#         if max(cal_mismatch(bot, to, mystack)[0:2]) < length_mismatch and cal_mismatch(bot, to, mystack)[2] < angle_mismatch:
#             filted_results.append(result)
#     filted_results.sort(key=sorting_func)
#     return filted_results

def filter_results(bottom: Atoms = None, top: Atoms = None, results: list = None, 
                   bot: Atoms = None, to: Atoms = None, stack: Atoms = None,
                   tag_mismatch: str = 'strain', filter_tag: str = 'normal strain' ,tag_z: bool = False, tag_value = 50, tag_type: str='left',
                   length_mismatch: float = 0.05, angle_mismatch: float = 0.05, 
                   normal_strain_mismatch: float = 0.05, shear_strain_mismatch: float = 0.05, stretch_mismatch: float =0.05,
                   sorting_func = lambda x: cal_atom_num(x.stack)) -> list:
    """
    Filters heterostructure results based on specified mismatches (e.g., strain, length, angle) and provides 
    options for additional adjustments (e.g., modification of z-components). Results can be sorted using a custom function.

    Args:
        bottom (Atoms, optional): Reference atoms object for the bottom layer (not supercell).
        top (Atoms, optional): Reference atoms object for the top layer (not supercell).
        results (list of Interface, optional): List of heterostructure interfaces to filter.
        bot (Atoms, optional): Supercell of the bottom layer. Used if not None.
        to (Atoms, optional): Supercell of the top layer. Used if not None.
        stack (Atoms, optional): Combined supercell of top and bottom. Used if not None.
        tag_mismatch (str, optional): Type of mismatch to filter ('strain' or 'percent'). Defaults to 'strain'.
        filter_tag (str, optional): Specific strain mismatch criterion ('normal strain', 'stretch'). Defaults to 'normal strain'.
        tag_z (bool, optional): If True, modifies the z-component of the lattice cell. Defaults to False.
        tag_value (float, optional): Value to set for the z-component if tag_z is True. Defaults to 50.
        tag_type (str, optional): Specifies which eigenvectors to consider ('left' or 'right'). Defaults to 'left'.
        length_mismatch (float, optional): Maximum allowable relative difference in length. Defaults to 0.05.
        angle_mismatch (float, optional): Maximum allowable angle difference in degrees. Defaults to 0.05.
        normal_strain_mismatch (float, optional): Threshold for normal strain mismatch. Defaults to 0.05.
        shear_strain_mismatch (float, optional): Threshold for shear strain mismatch. Defaults to 0.05.
        stretch_mismatch (float, optional): Threshold for eigenvalue stretch mismatch. Defaults to 0.05.
        sorting_func (callable, optional): Function to sort the filtered results based on a custom property. 
                                           Defaults to sorting by atom number in the stack.

    Raises:
        ValueError: The input [bottom, top, results] or [bot, to, stack] must not be None simultaneously.
        ValueError: If the tag_mismatch is neither 'strain' nor 'lattice constant'.
        ValueError: If the tag_type is neither 'left' nor 'right'.
        ValueError: If the filter_tag is neither 'normal strain' nor 'stretch'
        
    Returns:
        list: A list of filtered Interface objects, sorted based on the specified sorting function.
    
    Example:
        # Example of using filter_results to filter and sort heterostructure interfaces
        filtered_results = filter_results(bottom=Atoms(...), top=Atoms(...), results=my_results,
                                          sorting_func=lambda x: cal_atom_num(x.stack))
        for result in filtered_results:
            print(result)

    Note:
        This function assumes the 'Interface' type objects contain all necessary information 
        to build supercells and perform filtering based on the provided criteria.

    Todo:
        shear_strain_matrix: doesn't work well.

    """
    filted_results = []
    if not (all(value != None for value in [bottom, top, results]) or all(value != None for value in [bot, to, stack])):
        raise ValueError('The input [bottom, top, results] or [bot, to, stack] must not be None simultaneously.')
    # set this tag to avoid being infulenced by for loop.
    if bot == None and to == None and stack == None:
        tag_build = True
    for result in results:
        if tag_build:
            bot, to = build_supercells(bottom, top, result.M, result.N, result.angle, if_stack = False)
            stack = build_supercells(bottom, top, result.M, result.N, result.angle, result._weight)
        if tag_mismatch == 'strain':
            top_cell = array(to.get_cell()).T
            bottom_cell = array(bot.get_cell()).T
            stack_cell = array(stack.get_cell()).T
            if not tag_z:
                top_cell = adjust_matrix(top_cell, 2, 2, 0)
                top_cell[2,2] = tag_value
                bottom_cell = adjust_matrix(bottom_cell, 2, 2, 0)
                bottom_cell[2,2] = tag_value
                stack_cell = adjust_matrix(stack_cell, 2, 2, 0)
                stack_cell[2,2] = tag_value
            if tag_type == 'left':
                num = 0
            elif tag_type == 'right':
                num = 2
            else:
                raise ValueError(f'Unsupported tag_type {tag_type}.')
            defor = cal_deform_matrix(top_cell, stack_cell)
            relative_stretch = np.sqrt(cal_strain_matrix_root(bot = bot, to =to,stack =stack, tag_z =tag_z, tag_value = tag_value)[2+num][0]) - 1
            #strain_matrix = cal_strain_matrix(defor) bot = bot, to = to, stack = stack,
            principal_strain = cal_strain_matrix_root(bot = bot, to = to, stack = stack, tag_z = tag_z, tag_value = tag_value)[3+num][2]
            #strain_list = cal_principal_and_shear_strain(strain_matrix[num])

            # For eigenvalue, normal strain
            if filter_tag == 'normal strain':
                if  all(abs(value) < normal_strain_mismatch for value in principal_strain):
                    filted_results.append(result)
            elif filter_tag == 'stretch':
                if all(abs(value) < stretch_mismatch for value in relative_stretch):
                    filted_results.append(result)
            else:
                print(f'Unsupported filter_tag type {filter_tag}')
        elif tag_mismatch == 'lattice constant':
            if max(cal_mismatch(bot, to, stack)[0:2]) < length_mismatch and cal_mismatch(bot, to, stack)[2] < angle_mismatch:
                filted_results.append(result)
        else:
            print(f'Unsupported tag_mismatch type {tag_mismatch}')
    filted_results.sort(key=sorting_func)
    return filted_results


