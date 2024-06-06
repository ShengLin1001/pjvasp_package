from numpy import array, ndarray
from ase import Atoms

# move atoms
def move_atoms(atoms: Atoms = None,
               translate_matrix: list = [0.1, 0.1, 0.0]) -> Atoms :
    translate_matrix = array(translate_matrix)
    scaled = atoms.get_scaled_positions()
    scaled += translate_matrix
    atoms.set_scaled_positions(scaled)
    atoms.wrap()
    return atoms


