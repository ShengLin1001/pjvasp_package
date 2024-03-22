from numpy import array, ndarray
from ase import Atoms

# move atoms
def move_atoms(atoms: Atoms = None,
               translate_matrix: ndarray = array([0.1, 0.1, 0.0])) -> Atoms :
    scaled = atoms.get_scaled_positions()
    scaled += translate_matrix
    atoms.set_scaled_positions(scaled)
    atoms.wrap()
    return atoms
