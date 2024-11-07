"""
moveatom submodule

This submodule provides functions for moving atoms in a structure. It includes functions for translating atoms based on
a translation matrix. These functions are designed to streamline common tasks in materials science simulations and data handling.

Functions:
    - move_atoms: Move atoms in a structure by applying a translation matrix.
"""
from numpy import array, ndarray
from ase import Atoms
import numpy as np

# move atoms
# def move_atoms(atoms: Atoms = None,
#                translate_matrix: list = [0.1, 0.1, 0.0]) -> Atoms :
#     translate_matrix = array(translate_matrix)
#     scaled = atoms.get_scaled_positions()
#     scaled += translate_matrix
#     atoms.set_scaled_positions(scaled)
#     atoms.wrap()
#     return atoms

# move atoms
def move_atoms(atoms: Atoms = None,
               translate_matrix: ndarray = array([0.1, 0.1, 0.0]),
               if_scale_position: bool = True,
               atom_type: list = None,
               position_range: tuple = ((-np.inf, np.inf), (-np.inf, np.inf), (-np.inf, np.inf))) -> Atoms :
    """
    Move atoms by a translation matrix(in cartesian or direct coordination).
    The atoms to move can be filtered by atom type and position range(always in cartesian coordination).

    Args:
        atoms (Atoms, optional): The atomic structure to be translated, represented as an ASE `Atoms` object.
            This is the system of atoms to which the translation will be applied. Defaults to None.
        translate_matrix (ndarray, optional): The translation vector to apply to each selected atom. The translation
            can be applied in either Cartesian or scaled (fractional) coordinates, depending on the value of `if_scale_position`.
            Defaults to np.array([0.1, 0.1, 0.0]).
        if_scale_position (bool, optional): Whether to apply the translation in scaled (fractional) coordinates. If True, 
            the translation is applied to the scaled positions. If False, the translation is applied to Cartesian coordinates. 
            Defaults to True.
        atom_type (list, optional): A list of atomic symbols specifying which types of atoms should be moved (e.g., ['O', 'H']).
            If None, all atom types are considered. Defaults to None.
        position_range (tuple, optional): A tuple ((x_min, x_max), (y_min, y_max), (z_min, z_max)) specifying the range of positions 
            in Cartesian coordinates to filter which atoms should be moved. Only atoms within the specified position range will be moved. 
            Defaults to ((-np.inf, np.inf), (-np.inf, np.inf), (0.0, 1.0)) to restrict movement in the z-axis.

    Returns:
        Atoms: A new ASE `Atoms` object with translated positions for the selected atoms.

    Example:
        >>> from ase.build import molecule
        >>> atoms = molecule('H2O')
        >>> translated_atoms = move_atoms(atoms, translate_matrix=np.array([0.5, 0.5, 0.0]), 
                                          atom_type=['O'], position_range=((-np.inf, np.inf), (-np.inf, np.inf), (0.0, 1.0)))
        >>> print(translated_atoms.get_positions())

    """
    translate_matrix = array(translate_matrix)
    atoms_copy = atoms.copy()
    
    if if_scale_position:
        scaled = atoms_copy.get_scaled_positions()
        positions_to_modify = scaled
    else:
        positions = atoms_copy.get_positions()
        positions_to_modify = positions
    
    # Filter by atom type if specified
    
    if atom_type:
        indices1_to_move = [i for i, atom in enumerate(atoms) if atom.symbol in atom_type]
    else:
        indices1_to_move = list(range(len(atoms)))
    #print(indices1_to_move)

    # Filter by position range if specified
    if position_range:
        (x_min, x_max), (y_min, y_max), (z_min, z_max) = position_range
        indices_to_move2 = [
            i for i, atom in enumerate(atoms) if 
            (x_min <= atom.position[0] <= x_max) and
            (y_min <= atom.position[1] <= y_max) and
            (z_min <= atom.position[2] <= z_max)
        ]
    else:
        indices_to_move2 = list(range(len(atoms)))
    #print(indices_to_move2)

    # And operator
    indices_to_move = list(set(indices1_to_move) & set(indices_to_move2))
    #print(indices_to_move)

    # Apply translation to the selected atoms
    for i in indices_to_move:
        positions_to_modify[i] += translate_matrix
    
    if if_scale_position:
        atoms_copy.set_scaled_positions(positions_to_modify)
    else:
        atoms_copy.set_positions(positions_to_modify)
    
    #print(f'Move {len(indices_to_move)} atoms')
    atoms_copy.wrap()
    return atoms_copy


