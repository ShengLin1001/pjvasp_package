from numpy import array, ndarray, cross
from ase import Atoms
from ase.io.vasp import read_vasp, write_vasp
from ase.visualize import view
import spglib

def my_find_prim(atoms: Atoms = None, move_list = [0, 0, 0], check_direction_tag = True, scale_atoms = False) -> Atoms:
    """
    find primitive cell using spglib\n
    Convert to a format suitable for spglib
    """
    
    lattice = array(atoms.get_cell())
    points = array(atoms.get_scaled_positions())
    numbers = array(atoms.get_atomic_numbers())
    pbc = array(atoms.get_pbc())
    cell = (lattice, points, numbers)

    primitive_cell = spglib.standardize_cell(cell, to_primitive=1, no_idealize=1)
    # Convert the spglib output back to an ASE Atoms object
    primitive_atoms = Atoms(numbers = primitive_cell[2],
                            scaled_positions = primitive_cell[1],
                            cell = primitive_cell[0],
                            pbc=pbc)
    if check_direction_tag:
        primitive_atoms = check_direction(primitive_atoms, scale_atoms)
    primitive_atoms = move_atoms(primitive_atoms, move_list)
    return primitive_atoms

# move atoms
def move_atoms(atoms: Atoms = None,
               translate_matrix: ndarray = array([0.1, 0.1, 0.0])) -> Atoms :
    scaled = atoms.get_scaled_positions()
    scaled += translate_matrix
    atoms.set_scaled_positions(scaled)
    atoms.wrap()
    return atoms

# check the z-direction of cell is positive
def check_direction(atoms: Atoms = None, scale_atoms = False) -> Atoms:
    lattice = array(atoms.get_cell())
    
    points = array(atoms.get_positions())
    #print(points)
    if lattice[2,2] < 0:
        lattice *= -1
    else:
        lattice = lattice
    atoms.set_cell(lattice, scale_atoms= scale_atoms)
    #print(atoms.get_positions())
    atoms.set_positions(points)
    #print(atoms.get_positions())
    return atoms

# calculate the surface area in xy plane
def cal_area(atoms: Atoms = None) -> float:
    lattice = array(atoms.get_cell())
    a = lattice[0,:]
    b = lattice[1,:]
    area = abs(cross(a, b)[2])
    return area

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