"""
findprim submodule

This submodule provides functions for finding the primitive cell of a thin film
structure.

Functions:
    - my_find_prim: Finds the primitive cell of a thin film structure.
    - check_direction: Checks the direction of the cell for a 2D material.
"""

from numpy import array
from ase import Atoms
import spglib
from mymetal.universal.atom.moveatom import *

def my_find_prim(atoms: Atoms = None, move_list = [0, 0, 0], check_direction_tag = False, scale_atoms = False, to_primitive = 1,
                 if_fix_c: bool = False) -> Atoms:
    """Find the primitive or conventional cell of an ASE Atoms object using spglib.

    This function standardizes the given structure using spglib, returning either
    a primitive or conventional cell depending on `to_primitive`. For 2D systems or
    slabs or specified (ijk) plane-based bulk, it can fix the c-axis (z-direction) during symmetry reduction by temporarily
    adding vacuum and restoring the original height and zpositions afterward.

    Args:
        atoms (ase.Atoms): Input atomic structure.
        move_list (list[float], optional): Translation to apply after transformation
            (fractional coordinates). Default is [0, 0, 0].
        check_direction_tag (bool, optional): Whether to check and correct lattice
            orientation (recommended for 2D materials). Default is True. (Be carefully
            when setting to True!)
        scale_atoms (bool, optional): If True, scales atomic positions when changing
            cell dimensions. Default is False.
        to_primitive (int, optional): If 1, find the primitive cell; if 0, find the
            conventional cell. Default is 1.
        if_fix_c (bool, optional): If True, fix the c-axis (z-direction) for slabs or
            2D materials or specified (ijk) plane-based bulk by preserving the original cell height. Default is False.

    Returns:
        ase.Atoms: The standardized primitive or conventional cell.

    Notes:
        - When `if_fix_c=True`, the function:
            1. Temporarily centers the atoms with added vacuum along z.
            2. Runs spglib standardization.
            3. Restores the original z positions and c-axis length.
            4. This setting should only be used for specified (ijk) plane-based bulk.
        - When `check_direction_tag=True`, the lattice orientation is verified via
          `check_direction()`.
        - The final structure is sorted by atomic number, translated by `move_list`,
          and wrapped into the periodic cell.

    Example:
        >>> prim = my_find_prim(atoms, to_primitive=1)
        >>> conv = my_find_prim(atoms, to_primitive=0, if_fix_c=True)
    """

    

    # only find prim in xy plane for bulk
    # This part first add vacuum, then find prim, finally set back the cell[:,2] and move the atoms to original z positions
    if if_fix_c:
        # Find primitive cell of the slab
        # Save original cell/original positions
        old_cell = np.array(atoms.get_cell())
        old_positions = np.array(atoms.get_positions())
        old_min_z = np.min(old_positions[:, 2])
        atoms.center(vacuum=7.5, axis=2)


    # pass to spglib to find primitive cell
    lattice = array(atoms.get_cell())
    points = array(atoms.get_scaled_positions())
    numbers = array(atoms.get_atomic_numbers())
    pbc = array(atoms.get_pbc())
    cell = (lattice, points, numbers)
    # Important: pbc = [True, True, True], else spglib may treat it as 0D or 1D system, and raise some unexpected error!

    primitive_cell = spglib.standardize_cell(cell, to_primitive=to_primitive, no_idealize=1)
    # Convert the spglib output back to an ASE Atoms object
    primitive_atoms = Atoms(numbers = primitive_cell[2],
                            scaled_positions = primitive_cell[1],
                            cell = primitive_cell[0],
                            pbc=pbc)
    
    if if_fix_c:
        # Move the slab to original z positions
        new_cell = np.array(primitive_atoms.get_cell())
        new_positions = np.array(primitive_atoms.get_positions())
        new_min_z = np.min(new_positions[:, 2])
        shift_z = new_min_z - old_min_z
        new_positions[:, 2] -= shift_z

        # set cell, a, b unchanged, c changed to original c
        new_cell[2, :] = old_cell[2, :]
        primitive_atoms.set_cell(new_cell, scale_atoms=False)
        primitive_atoms.set_positions(new_positions)


    if check_direction_tag:
        primitive_atoms = check_direction(primitive_atoms, scale_atoms)
    primitive_atoms = move_atoms(primitive_atoms, move_list)
    primitive_atoms = primitive_atoms[primitive_atoms.numbers.argsort()]
    primitive_atoms.wrap()
    return primitive_atoms


# check the z-direction of cell is positive for 2D material primitive cell
def check_direction(atom: Atoms = None, scale_atoms: bool = False, align_row: int = 0, align_to: list = [1, 0, 0]) -> Atoms:
    """
    Check and adjust the z-direction of a 2D material's primitive cell to ensure it is positive. 
    Additionally, this function aligns the first lattice vector (a1) to the x-direction by rotating the structure.

    Args:
        atom (Atoms): An ASE Atoms object representing the crystal structure. 
                       It should be a valid Atoms object with a defined lattice.
        scale_atoms (bool, optional): Whether to scale atoms after modifying the lattice. 
                             If True, the atomic positions will be scaled to match the new lattice parameters. Default is False.
        align_row (int, optional): The index of the row (lattice vector) to be aligned to the x-direction. 
                                   Default is 0 (i.e., align the first lattice vector).
        align_to (list, optional): The target direction for alignment as a list of 3 values representing a vector (e.g., [1, 0, 0] for the x-direction). 
                                   Default is [1, 0, 0], which aligns the selected lattice vector to the x-axis.


    Returns:
        Atoms: The modified Atoms object with the adjusted lattice and rotated structure.

    Raises:
        ValueError: If the provided `atom` is None or not a valid ASE Atoms object.
        TypeError: If `scale_atoms` is not a boolean.

    Example:
        >>> atoms = check_direction(my_atoms, scale_atoms=True)
    """
    # Parameter checks
    if atom is None:
        raise ValueError("Input atom must be a valid ASE Atoms object.")
    if not isinstance(scale_atoms, bool):
        raise TypeError("scale_atoms must be a boolean value.")
    atoms = atom.copy()
    lattice = np.array(atoms.cell.copy())
    # Check and ensure the determinant of the lattice is positive
    determinant = np.linalg.det(lattice)
    if determinant < 0:
        lattice = lattice * -1
        atoms.set_cell(lattice, scale_atoms=scale_atoms)

    # Ensure z-direction of the lattice is positive
    if lattice[2, 2] < 0:
        lattice[0, :] *= -1  # Invert first row
        lattice[2, :] *= -1  # Invert third row
        atoms.set_cell(lattice, scale_atoms=scale_atoms)
    
    # Align the a1 vector to the x-direction
    a1 = lattice[align_row]
    b1 = align_to  
    a1 = np.asarray(a1, dtype=float) / np.linalg.norm(a1)
    b1 = np.asarray(b1, dtype=float) / np.linalg.norm(b1)
    # Calculate cross-product to determine rotation axis
    c1 = np.cross(a1, b1)
    if np.linalg.norm(c1) > 1e-10:
        c1 /= np.linalg.norm(c1)
        # Rotate the atom structure
        # it's important to rotate atoms, not atom
        atoms.rotate(np.arccos(np.dot(a1, b1)) / np.pi * 180, c1, rotate_cell=True)
    atoms.wrap()
    return atoms

