from ase import Atoms, Atom
from numpy import array, ndarray, argsort
from ase.build import add_adsorbate

# def add_hydro_atoms(atoms:  Atoms = None,
#                     add_symbol: chr = 'H', 
#                     added_symbol: chr = 'O', 
#                     surf_range: ndarray = array([0.1, 11.1]), 
#                     shift_distance: ndarray = array([0.0, 0.0, 2.5]),
#                     surf_direction: int = 2) -> Atoms:
#     """
#     Add hydrogen atoms to the surface of the atoms object.
#     """
#     atoms = atoms.copy()
#     #collect = atoms.copy()
#     for atom in atoms:
#         if atom.position[surf_direction] > surf_range[0] and atom.position[surf_direction] < surf_range[1] and atom.symbol == added_symbol:
#             #print(atom.index, atom.symbol, atom.position)
#             atoms.append(Atom(symbol=add_symbol, position=atom.position + shift_distance))
#     atoms = atoms[atoms.numbers.argsort()]
#     return atoms

def add_hydro_atoms(my_atoms:  Atoms = None,
                    add_symbol: chr = 'H', 
                    added_symbol: chr = 'O', 
                    surf_range: ndarray = array([0.1, 11.1]), 
                    shift_distance: ndarray = array([0.0, 0.0, 2.5]),
                    surf_direction: int = 2,
                    reverse_every: bool = False,
                    reverse_length: int = 1,
                    fix_directions: list = [2],
                    sort_direction: int = 0,
                    offset = None,
                    mol_index = 0,
                    manual: bool = False) -> Atoms:
    """
    Adds hydrogen atoms to the surface of the specified atomic structure.
    ------
    Args:
        my_atoms (Atoms, optional): The atomic structure to which hydrogen atoms will be added.
        add_symbol (str, optional): The symbol of the atom to add, defaults to 'H' for hydrogen. (new)
        added_symbol (str, optional): The symbol of the surface atom to which hydrogen atoms are added, defaults to 'O' for oxygen. (old)
        surf_range (ndarray, optional): A 1D array specifying the range along the surface direction within which atoms are considered for hydrogen addition.
        shift_distance (ndarray, optional): A 3D vector specifying the shift in position to place the hydrogen atom relative to the surface atom.
        surf_direction (int, optional): The axis index (0, 1, or 2 corresponding to x, y, or z) along which the surface is oriented.
        reverse_every (bool, optional): A flag to reverse the shift direction after each addition, defaults to False.
        fix_directions (list, optional): A list of indices of the axes where the shift distance should not be reversed.
        sort_direction (int, optional): The axis index to sort the surface atoms by their coordinates.
        offset (optional): Offset for placing the added atom, can be used for complex surface structures.
        mol_index (int, optional): Molecular index to specify which molecule to modify, useful in multi-molecular systems.

    Returns:
        Atoms: The modified atomic structure with added hydrogen atoms.

    Raises:
        TypeError: If `my_atoms` is not an instance of `Atoms`.
        TypeError: If `add_symbol` or `added_symbol` are not instances of `str`.
        ValueError: If `surf_range` or `shift_distance` do not contain exactly two/three elements.
        ValueError: If `surf_direction` or `sort_direction` are not one of the valid indices (0, 1, or 2).

    Usages:
        added_atoms = my_add_hydro_atoms(temp, add_symbol = 'H', added_symbol = 'O', surf_range=[16.7, 17], shift_distance = [0.9441063, 0.2422, 0], 
                                 reverse_every=True)
        added_atoms = add_hydro_atoms(temp, add_symbol = 'H', added_symbol = 'O', surf_range=[16.7, 17], shift_distance = [0, 0, 0.9584])
        added_atoms = add_hydro_atoms(added_atoms, add_symbol = 'H', added_symbol = 'O', surf_range=[7.4, 7.6], shift_distance = [0, 0, -0.9584])
    """
    if not isinstance(my_atoms, Atoms):
        raise TypeError("Expected my_atoms to be an instance of Atoms")
    if not isinstance(add_symbol, str) or not isinstance(added_symbol, str):
        raise TypeError("Expected add_symbol, added_symbol to be the instance of chr")
    if len(surf_range) != 2:
        raise ValueError("surf_range must contain exactly two elements")
    if len(shift_distance) != 3:
        raise ValueError("shift_distance must contain exactly three elements")
    if surf_direction not in [0, 1, 2] or sort_direction not in [0, 1, 2]:
        raise ValueError("surf_direction and sort_direction must be 0, 1, or 2")
    #print(surf_range)
    surf_range = array(surf_range)
    shift_distance = array(shift_distance)
    atoms = my_atoms.copy()
    atoms.set_pbc(True)
    shift = shift_distance
    collect = []
    # test top position
    positions = atoms.get_positions()
    position = positions[:,fix_directions]
    position.sort()
    max_position = float(max(position))

    for atom in atoms:
        if atom.position[surf_direction] > surf_range[0] and atom.position[surf_direction] < surf_range[1] and atom.symbol == added_symbol:
            collect.append((atom.index, atom.position))

    sorted_indices = argsort([item[sort_direction] for index,item in collect])
    sorted_collect = [collect[i] for i in sorted_indices]
    #print(sorted_collect)
    count = 0
    for index,item in sorted_collect:
        item = array(item)
        count = count + 1
        if manual:
            # Copy atoms for visualization and highlight the specific atom
            view(atoms)
            # Manual input for direction change
            # Request user input for direction modifiers
            direction_modifiers = input(f"For this adsorbate of {index}\nEnter direction modifiers as a list (e.g., [1, 1, -1]): ")
            if direction_modifiers.strip() == '':
                direction_modifiers = [1, 1, 1]
            else:
                direction_modifiers = eval(direction_modifiers)  # Converts input string to a list

            shift = shift_distance * array(direction_modifiers)
            add_adsorbate(slab = atoms, adsorbate = add_symbol, height = shift[2], position = tuple(shift[:2] +item[:2]),
                        offset = offset, mol_index = mol_index)
        else:
            #print(shift[2])
            #print(tuple(shift[:2] +item[:2]))
            add_adsorbate(slab = atoms, adsorbate = add_symbol, height = item[2] - max_position + shift[2], position = tuple(shift[:2] +item[:2]),
                        offset = offset, mol_index = mol_index)
            if reverse_every and count % reverse_length == 0:
                shift = -1. * shift
                for fix_direction in fix_directions:
                    shift[fix_direction] = shift_distance[fix_direction]
        #atoms.append(Atom(symbol=add_symbol, position=atom.position + shift))
    atoms = atoms[atoms.numbers.argsort()]
    atoms.wrap()
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
