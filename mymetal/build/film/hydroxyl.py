from ase import Atoms, Atom
from numpy import array, ndarray

def add_hydro_atoms(atoms:  Atoms = None,
                    add_symbol: chr = 'H', 
                    added_symbol: chr = 'O', 
                    surf_range: ndarray = array([0.1, 11.1]), 
                    shift_distance: ndarray = array([0.0, 0.0, 2.5]),
                    surf_direction: int = 2) -> Atoms:
    """
    Add hydrogen atoms to the surface of the atoms object.
    """
    atoms = atoms.copy()
    #collect = atoms.copy()
    for atom in atoms:
        if atom.position[surf_direction] > surf_range[0] and atom.position[surf_direction] < surf_range[1] and atom.symbol == added_symbol:
            #print(atom.index, atom.symbol, atom.position)
            atoms.append(Atom(symbol=add_symbol, position=atom.position + shift_distance))
    atoms = atoms[atoms.numbers.argsort()]
    return atoms


## usage
# # INPUT the SiO2 model
# atoms = read_vasp('SiO2.poscar')
# print(atoms)
# # view(atoms)


# # SLAB model
# temp = atoms.repeat([1, 1, 2])
# print(temp)

# del temp[[atom.index for atom in temp if atom.position[2] < 0.9 or atom.position[2] > 10.18 ] ]
# temp = temp[temp.numbers.argsort()]
# temp.center(vacuum = 5, axis = 2)
# print(temp)
# #view(temp)

# # ADD H atoms
# print(temp)
# added_atoms = add_hydro_atoms(temp, add_symbol = 'H', added_symbol = 'O', surf_range=[14.2, 15], shift_distance = [0, 0, 0.9584])
# added_atoms = add_hydro_atoms(added_atoms, add_symbol = 'H', added_symbol = 'O', surf_range=[4.9, 5.1], shift_distance = [0, 0, -0.9584])
# added_atoms.center(vacuum=5, axis=2)
# print(added_atoms)
# view(added_atoms)
