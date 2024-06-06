from mymetal.calculate.calmechanics.strain import cal_strain_matrix_root
from ase import Atoms
from hetbuilder import Interface
import numpy as np

def cal_relative_stretch(bottom: Atoms = None, top: Atoms = None, result: Interface = None, 
                            bot: Atoms = None, to: Atoms = None, stack: Atoms = None,
                            tag_z: bool = False, tag_value = 50, tag_type: str = 'left') -> list:
    """
    Calculate the relative stretch factor and its direction based on the deformation between two configurations.

    Args:
        bottom (Atoms, optional): Reference bottom atoms object for initial configuration (not supercell).
        top (Atoms, optional): Reference top atoms object for initial configuration (not supercell).
        result (Interface, optional): Interface object containing deformation information.
        bot (Atoms, optional): Supercell of the bottom layer, already calculated.
        to (Atoms, optional): Supercell of the top layer, already calculated.
        stack (Atoms, optional): Combined supercell, already calculated.
        tag_z (bool, optional): If True, modify the z-component of the lattice cell. Defaults to False.
        tag_value (float, optional): Value to set for the z-component if tag_z is True. Defaults to 50.
        tag_type (str, optional): Specifies the orientation ('left' or 'right') for the stretch calculation. Defaults to 'left'.

    Returns:
        list: [relative stretch factor, direction of stretch]
              relative stretch factor is the stretch relative to the original configuration minus one.
              direction of stretch is a vector indicating the principal direction of stretch.

    Raises:
        ValueError: If neither input configuration sets are provided, or an unsupported tag_type is given.
        ValueError: The input [bottom, top, result] or [bot, to, stack] must not be None simultaneously.
    """
    if not (all(value != None for value in [bottom, top, result]) or all(value != None for value in [bot, to, stack])):
        raise ValueError('The input [bottom, top, result] or [bot, to, stack] must not be None simultaneously.')
    stretch = cal_stretch(bottom, top, result,bot = bot, to = to, stack = stack, tag_z = tag_z, tag_value=tag_value, tag_type = tag_type)
    relative_stretch_factor = [stretch[0][0] - 1, stretch[0][1] - 1]
    relative_stretch_direction = stretch[1]
    return [relative_stretch_factor, relative_stretch_direction]

def cal_stretch(bottom: Atoms = None, top: Atoms = None, result: Interface = None, 
                            bot: Atoms = None, to: Atoms = None, stack: Atoms = None,
                            tag_z: bool = False, tag_value = 50, tag_type: str = 'left') -> list:
    """
    Calculate the absolute stretch factor and direction based on the deformation between initial and final configurations.

    Args:
        bottom (Atoms, optional): Reference bottom atoms object for initial configuration (not supercell).
        top (Atoms, optional): Reference top atoms object for initial configuration (not supercell).
        result (Interface, optional): Interface object containing deformation information.
        bot (Atoms, optional): Supercell of the bottom layer, already calculated.
        to (Atoms, optional): Supercell of the top layer, already calculated.
        stack (Atoms, optional): Combined supercell, already calculated.
        tag_z (bool, optional): If True, modify the z-component of the lattice cell. Defaults to False.
        tag_value (float, optional): Value to set for the z-component if tag_z is True. Defaults to 50.
        tag_type (str, optional): Specifies the orientation ('left' or 'right') for the stretch calculation. Defaults to 'left'.

    Returns:
        list: [stretch factor, direction of stretch]
              stretch factor is a scalar representing the degree of stretching.
              direction of stretch is a vector indicating the principal direction of stretch.

    Raises:
        ValueError: If neither input configuration sets are provided, or an unsupported tag_type is given.
        ValueError: The input [bottom, top, result] or [bot, to, stack] must not be None simultaneously.
    """
    if not (all(value != None for value in [bottom, top, result]) or all(value != None for value in [bot, to, stack])):
        raise ValueError('The input [bottom, top, result] or [bot, to, stack] must not be None simultaneously.')

    if tag_type in ['left', 'l', 'Left', 'L']:
        num = 0
    elif tag_type in ['right', 'r', 'Right', 'R']:
        num = 2
    else:
        raise ValueError('Unsupported type.')
    temp = cal_strain_matrix_root(bottom, top, result,bot = bot, to = to, stack = stack, tag_z = tag_z, tag_value=tag_value)
    stretch_all = [strain[2+num] for strain in temp]
    stretch_factor = [np.sqrt(stretch[0]) for stretch in stretch_all]
    stretch_direction = [stretch[1] for stretch in stretch_all]
    return [stretch_factor, stretch_direction]


            
