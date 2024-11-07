"""
strain module

This module provides functions for calculating various types of strain, including 
principal strain, shear strain, and von Mises strain. Strain analysis is vital 
in the study of material plasticity and mechanical properties.

Functions:
    - cal_principal_and_shear_strain: Calculate the principal and shear strains.
    - cal_principal_and_shear_strain_root: Calculate the square root of principal and shear strains.
    - cal_strain_matrix: Calculate the strain matrix for a given structure.
    - cal_strain_matrix_root: Calculate the square root of the strain matrix.
    - cal_von_mises_strain: Calculate the von Mises strain.
"""

from ase import Atoms
import numpy as np
from numpy import array
from mymetal.build.film.findhetero import build_supercells
from mymetal.calculate.calmechanics.deformation import cal_deform_matrix
from mymetal.universial.matrix.adjust import adjust_matrix
from hetbuilder import Interface

def cal_principal_and_shear_strain_root(bottom: Atoms = None, top: Atoms = None, result: Interface = None, 
                            bot: Atoms = None, to: Atoms = None, stack: Atoms = None,
                            tag_z: bool = False, tag_value = 50, tag_type: str = 'left') -> list:
    """
    Calculate the principal and shear strains for a given heterostructure configuration.

    This function computes the principal and shear strains based on the deformation matrix
    between a reference configuration and a deformed configuration. It uses the deformation
    matrix to determine the eigenvalues and eigenvectors, which represent the principal
    strains and their directions, as well as the shear strain matrix.

    Args:
        bottom (Atoms): Atoms object for the bottom layer in its initial configuration.
        top (Atoms): Atoms object for the top layer in its initial configuration.
        result (Interface): Interface object that contains information about the heterostructure.
        bot (Atoms): Atoms object for the bottom layer in its supercell configuration.
        to (Atoms): Atoms object for the top layer in its supercell configuration.
        stack (Atoms): Atoms object representing the combined supercell after stacking.
        tag_z (bool): Flag to indicate whether to adjust the z-component of the cell matrix.
        tag_value (float): The value to set for the z-component if tag_z is True.
        tag_type (str): Specifies which type of strain (left or right) to consider. Valid options are 'left' or 'right'.

    Returns:
        list: A list containing the principal strain values, principal strain directions, and the shear strain matrix.

    Raises:
        ValueError: If an unsupported tag_type is provided.
        ValueError: The input [bottom, top, result] or [bot, to, stack] must not be None simultaneously.

    Examples:
        >>> # Assuming that `bottom`, `top`, `result`, `bot`, `to`, and `stack` are defined elsewhere
        >>> strains = cal_principal_and_shear_strain_root(bottom, top, result, bot, to, stack, tag_z=False, tag_value=50, tag_type='left')
        >>> print(strains)
        # Output: [principal_strain_values, principal_strain_directions, shear_strain_matrix]
    """
    if not (all(value != None for value in [bottom, top, result]) or all(value != None for value in [bot, to, stack])):
        raise ValueError('The input [bottom, top, results] or [bot, to, stack] must not be None simultaneously.')
    if tag_type == 'left':
        num = 0
    elif tag_type == 'right':
        num = 2
    else:
        raise ValueError(f'Unsupported tag_type {tag_type}.')
    temp = cal_strain_matrix_root(bottom, top, result, bot = bot, to = to, stack = stack, tag_z = tag_z, 
                                                    tag_value = tag_value)
    principal_shear_strain = [strain[3+num][2] for strain in temp]
    principal_strain_matrix = [strain[0] for strain in principal_shear_strain]
    shear_strain_matrix = [strain[1] for strain in principal_shear_strain]
    priciple_strain = [strain[2] for strain in principal_shear_strain]
    priciple_direction = [strain[3] for strain in principal_shear_strain]
    return [priciple_strain, priciple_direction, shear_strain_matrix]

def cal_strain_matrix_root(bottom: Atoms = None, top: Atoms = None, result: Interface = None, 
                            bot: Atoms = None, to: Atoms = None, stack: Atoms = None,
                            tag_z: bool = False, tag_value = 50) -> list:
    """
    Calculate deformation and strain matrices for a given pair of bottom and top 
    Atoms objects along with their transformed supercell versions. Optionally adjusts
    the z-component of the cell matrices.

    Args:
        bottom (Atoms, optional): The bottom Atoms object from which the supercells are built.
        top (Atoms, optional): The top Atoms object from which the supercells are built.
        result (Interface, optional): The Interface object containing transformation matrices M and N.
        bot (Atoms, optional): The transformed bottom Atoms object.
        to (Atoms, optional): The transformed top Atoms object.
        stack (Atoms, optional): The transformed stack Atoms object.
        tag_z (bool, optional): If True, modify the z-component of the cell matrices. Defaults to False.
        tag_value (float, optional): The value to set for the z-components if tag_z is True. Defaults to 50.

    Returns:
        Two lists [top, bottom]
        list: A list containing the deformation matrix, strain matrix, eigenvalues and eigenvectors for left and right multiplication, 
        and principal and shear strains for both left and right matrices:
        - defor: The deformation matrix calculated as the dot product of the top cell and inverse of the stack cell.
        - strain_matrix: The strain matrix calculated from the deformation matrix.
        - [eigenvalues_left, eigenvectors_left]: Eigenvalues and eigenvectors of C = Ft* F
        - strain_list_left: Principal and shear strains calculated from the left strain matrix.
        - [eigenvalues_right, eigenvectors_right]: Eigenvalues and eigenvectors of F* Ft
        - strain_list_right: Principal and shear strains calculated from the right strain matrix.

    Example:
        # Define initial and transformed atoms
        bottom = Atoms(...)
        top = Atoms(...)
        result = Interface(...)
        # Calculate strain matrices
        strains = cal_strain_matrix_root(bottom, top, result)
        print(strains)

    Note:
        If neither `bot`, `to`, nor `stack` are provided, the function will build the supercells
        using the given `bottom`, `top`, and `result` objects.
    """
    if not (all(value != None for value in [bottom, top, result]) or all(value != None for value in [bot, to, stack])):
        raise ValueError('The input [bottom, top, results] or [bot, to, stack] must not be None simultaneously.')
    if bot == None and to == None and stack == None:
        bot, to = build_supercells(bottom, top, result.M, result.N, result.angle, if_stack = False)
        stack = build_supercells(bottom, top, result.M, result.N, result.angle, result._weight)
    top_cell = array(to.get_cell()).T
    bottom_cell = array(bot.get_cell()).T
    stack_cell = array(stack.get_cell()).T
    if not tag_z:
        top_cell = adjust_matrix(top_cell, 2, 2, 0)
        bottom_cell = adjust_matrix(bottom_cell, 2, 2, 0)
        stack_cell = adjust_matrix(stack_cell, 2, 2, 0)
        top_cell[2,2] = tag_value
        bottom_cell[2,2] = tag_value
        stack_cell[2,2] = tag_value
    defor_all = [cal_deform_matrix(top_cell, stack_cell), cal_deform_matrix(bottom_cell, stack_cell)]
    content = []
    for defor in defor_all:
        eigenvalues_left, eigenvectors_left = np.linalg.eig(np.dot(defor.T,defor))
        eigenvalues_right, eigenvectors_right = np.linalg.eig(np.dot(defor,defor.T))
        strain_matrix = cal_strain_matrix(defor)
        strain_list_left = cal_principal_and_shear_strain(strain_matrix[0])
        strain_list_right = cal_principal_and_shear_strain(strain_matrix[1])
        content.append([defor, strain_matrix, [eigenvalues_left, eigenvectors_left],strain_list_left, [eigenvalues_right, eigenvectors_right], strain_list_right])
    return content

def cal_strain_matrix(deformation_matrix: np.array = None) -> list:
    """To calculate the lagrangian and euler strain matrix.

    Args:
        initial_atoms (Atoms): initial configuration
        deformed_atoms (Atoms): deformed configuration

    Raises:
        ValueError: Input matrix cannot be None.
        ValueError: Input matrix must be `n * n (0 < n <= 3)`

    Returns:
        list: [lagrangian strain, euler strain]
    """

    if deformation_matrix is None:
        raise ValueError("Input matrix cannot be None.")
    if deformation_matrix.shape[0] > 3 or deformation_matrix.shape[0] != deformation_matrix.shape[1] or deformation_matrix.size == 0:
        raise ValueError(f"Expected a n*n np.array, but it's {deformation_matrix.shape[0]} * {deformation_matrix.shape[1]} if not empty.")

    F = deformation_matrix

    # For small defoemation
    # left Cauchy-Green deformation tensor
    # C = F^T * F
    lc = np.dot(F.T, F)
    # Lagrangian strain tensor
    # ε = 1/2 * (F^T * F - I)
    E = 0.5 * (lc - np.eye(F.shape[0]))


    # For large deformation
    # right Cauchy-Green deformation tensor
    rc = np.dot(F, F.T)
    # Euler strain tensor
    # ε = 1/2 * (I - F * F^T)
    e = 0.5 * (np.eye(F.shape[0]) - np.linalg.inv(rc))
    
    return [E, e]



def cal_principal_and_shear_strain(strain_matrix: np.array = None) -> list:
    """To calculate the principal and shear strain

    Args:
        strain_matrix (np.array): a 3*3 matrix

    Raises:
        ValueError: Input matrix cannot be None.
        ValueError: Input matrix must be `n * n (0 < n <= 3)`

    Returns:
        list: [principal_strain_matrix, shear_strain_matrix, eigenvalues, eigenvectors]
    """

    if strain_matrix is None:
        raise ValueError("Input matrix cannot be None.")
    if strain_matrix.shape[0] > 3 or strain_matrix.shape[0] != strain_matrix.shape[1] or strain_matrix.size == 0:
        raise ValueError(f"Expected a n*n np.array, but it's {strain_matrix.shape[0]} * {strain_matrix.shape[1]} if not empty.")

    # Principal strain, direction
    eigenvalues, eigenvectors = np.linalg.eig(strain_matrix)
    principal_strain_matrix = np.diag(eigenvalues)
    
    # devoratic strain
    shear_strain_matrix = strain_matrix - np.diag(np.diag(strain_matrix))
    
    return principal_strain_matrix, shear_strain_matrix, eigenvalues, eigenvectors

def cal_von_mises_strain(strain):
    """Equivalent strain to Von Mises Stress."""
    eps = strain - 1 / 3 * np.trace(strain) * np.identity(3)

    return np.sqrt(np.sum(eps * eps) * 2 / 3)

