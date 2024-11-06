from numpy import array
from ase import Atoms
import spglib
from mymetal.universial.atom.moveatom import *

def my_find_prim(atoms: Atoms = None, move_list = [0, 0, 0], check_direction_tag = True, scale_atoms = False, to_primitive = 1) -> Atoms:
    """
    find primitive cell using spglib\n
    Convert to a format suitable for spglib\n
    if the material is not 2D material, please turn off the check_direction tag\n
    this function could be used to find conventional cell by controlling the to_primitive tag\n
    """
    
    lattice = array(atoms.get_cell())
    points = array(atoms.get_scaled_positions())
    numbers = array(atoms.get_atomic_numbers())
    pbc = array(atoms.get_pbc())
    cell = (lattice, points, numbers)

    primitive_cell = spglib.standardize_cell(cell, to_primitive=to_primitive, no_idealize=1)
    # Convert the spglib output back to an ASE Atoms object
    primitive_atoms = Atoms(numbers = primitive_cell[2],
                            scaled_positions = primitive_cell[1],
                            cell = primitive_cell[0],
                            pbc=pbc)
    if check_direction_tag:
        primitive_atoms = check_direction(primitive_atoms, scale_atoms)
    primitive_atoms = move_atoms(primitive_atoms, move_list)
    primitive_atoms = primitive_atoms[primitive_atoms.numbers.argsort()]
    return primitive_atoms


# check the z-direction of cell is positive for 2D material primitive cell
def check_direction(atoms: Atoms = None, scale_atoms: bool = False, align_row: int = 0, align_to: list = [1, 0, 0]) -> Atoms:
    """
    Check and adjust the z-direction of a 2D material's primitive cell to ensure it is positive. 
    Additionally, this function aligns a specified lattice vector to the x-direction by rotating the structure.
    It also ensures the lattice determinant is positive.

    Args:
        atoms (Atoms): An ASE Atoms object representing the crystal structure. 
                       The object should contain a valid lattice with the corresponding atomic positions.
        scale_atoms (bool, optional): Whether to scale atoms after modifying the lattice. 
                                       If True, the atomic positions will be scaled to match the new lattice parameters. 
                                       Default is False.
        align_row (int, optional): The index of the row (lattice vector) to be aligned to the x-direction. 
                                   Default is 0 (i.e., align the first lattice vector).
        align_to (list, optional): The target direction for alignment as a list of 3 values representing a vector (e.g., [1, 0, 0] for the x-direction). 
                                   Default is [1, 0, 0], which aligns the selected lattice vector to the x-axis.

    Returns:
        Atoms: The modified ASE Atoms object with the adjusted lattice and rotated structure.

    Raises:
        ValueError: If the `atom` is None or not a valid ASE Atoms object.
        TypeError: If `scale_atoms` is not a boolean or `align_to` is not a list of 3 numerical values.

    Example:
        >>> atoms = check_direction(my_atoms, scale_atoms=True, align_row=0, align_to=[1, 0, 0])
        >>> print(atoms.cell)
    """
    if atoms is None:
        raise ValueError("Input atoms must be a valid ASE Atoms object.")
    if not isinstance(scale_atoms, bool):
        raise TypeError("scale_atoms must be a boolean value.")
    atom = atoms.copy()
    lattice = array(atom.cell.copy())
    determinant = np.linalg.det(lattice)  # 计算行列式
    if determinant < 0:
        lattice = lattice * -1
        atom.set_cell(lattice, scale_atoms= scale_atoms)
    if lattice[2,2] < 0:
        lattice[0, :] *= -1  # keep the determinant positive
        lattice[2, :] *= -1  
        atom.set_cell(lattice, scale_atoms= scale_atoms)
    # rotate the a1 to the x-direction
    a1 = lattice[align_row]
    b1 = align_to
    a1 = np.asarray(a1, dtype=float) / np.linalg.norm(a1)
    b1 = np.asarray(b1, dtype=float) / np.linalg.norm(b1)
    c1 = np.cross(a1, b1)
    c1 /= np.linalg.norm(c1)
    #print(np.arccos(np.dot(a1, b1))/np.pi*180)
    atom.rotate(np.arccos(np.dot(a1, b1))/np.pi*180, c1, rotate_cell=True)
    return atom


## usage
# str1 = read_vasp('POSCAR1')
# str2 = read_vasp('POSCAR2')
# #str2 *= [2, 2, 1]
# #print(str1)
# #print(str2)
# prim1 = my_find_prim(str1, [0.001, 0.001, 0.01])
# prim2 = my_find_prim(str2, [0.001, 0.001, 0.01])

# print(prim1)
# print(prim2)

# write_vasp('./POSCAR1_prim', prim1, label = 'SiO', direct= False, vasp5= True)
# write_vasp('./POSCAR2_prim', prim2, label = 'Au', direct= False, vasp5= True)

# super1 = read_vasp('SUPERLATTICE_0001_18.4_0036.vasp')
# print(super1)

# print(cal_area(prim1))
# print(cal_area(prim2))
# print(cal_area(super1))