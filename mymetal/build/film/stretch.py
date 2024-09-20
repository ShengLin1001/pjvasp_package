from numpy import sqrt, ndarray, array
from ase import Atoms
from ase.build import bulk, surface
# For version independency, it could be instead by 'from ase.io.vasp import *'
import spglib
from mymetal.build.film.findprim import my_find_prim
import os

def generate_film(   symbols: str = None,                 # str
                structure: str = None,              # str
                num_layers: int = None,            # int
                replic_z: int = None,              # int
                my_vacuum: float = 7.5,           # float
                slice_plane: tuple[int] = (0, 0, 1),  # Tuple of ints
                my_tol: float = 1e-6,              # float
                my_periodic: bool = False,          # bool
                a_fcc: float = 2.95*sqrt(2.0),     # float
                a_hcp: float = 2.95,               # float
                my_covera: float = sqrt(8.0/3.0),  # float
                move_atom: list = [0.1, 0.1, 0.0],
                before_move_atom: list = [0.05, 0.05, 0.05],
                number_per_layer: float = 1.0,
                bulk_atoms: Atoms = None
                ) -> Atoms:
    """
    please see the ASE/cut function\n
    parameters : 1-8 line: general setting\n
                  9 10-11 line: fcc, hcp parameters\n
    when we use surface() function, my_bulk must be a conventional cell, not a primitive cell, so set cubic=True
    Generates a thin film structure (slab) from bulk material and applies modifications. Please see the ASE/cut function\n

    The function allows the user to create a thin film from an bulk structure from `bulk_atoms` or specify the lattice constant of only supported hcp/fcc structure
    , adjust the number of layers, apply vacuum padding, and move atoms before or after generating the slab. The slab is created along the given slice plane.

    Args:
        symbols (str): The chemical symbol of the material to generate the bulk (e.g., 'Si'). So far, only one element is supported.
        structure (str): The crystal structure type ('fcc' or 'hcp').
        num_layers (int): The number of layers in the film. Either this or `replic_z` must be specified, superior to replic_z.
        replic_z (int): The number of repetitions in the z direction (along the surface normal). Either this or `num_layers` must be specified.
        my_vacuum (float): The amount of vacuum to add to the slab, two sizes of the vacuum are added to the top and bottom of the slab (default is 7.5).
        slice_plane (tuple[int]): The slicing plane defined as a tuple of Miller indices (default is (0, 0, 1)).
        my_tol (float): The tolerance for detecting layers in the slab (default is 1e-6).
        my_periodic (bool): Whether to apply periodic boundary conditions to the slab (default is False).
        a_fcc (float): The lattice constant for FCC structures (default is 2.95 * sqrt(2.0)).
        a_hcp (float): The lattice constant for HCP structures (default is 2.95). 
        my_covera (float): The c/a ratio for HCP structures (default is sqrt(8.0/3.0)).
        move_atom (list): A list of x, y, z translations to apply to the atoms after slab generation (default is [0.1, 0.1, 0.0]), in lattice unit.
        before_move_atom (list): A list of x, y, z translations to apply to the atoms before slab generation (default is [0.05, 0.05, 0.05]), in lattice unit.
        number_per_layer (float): The number of atoms per layer (default is 1.0).
        bulk_atoms (Atoms, optional): A pre-existing bulk structure to use. If not provided, a new bulk is generated.

    Returns:
        Atoms: The final slab structure with applied modifications.

    Raises:
        ValueError: If the `structure` is not 'fcc' or 'hcp', or if neither `num_layers` nor `replic_z` is provided.

    Example:
        >>> slab = generate_film(symbols='Cu', structure='fcc', num_layers=10, my_vacuum=10)
        >>> print(slab)
    """

    if bulk_atoms:
        my_bulk = my_find_prim(bulk_atoms, to_primitive=0)
    else:
        if structure == 'fcc':
            my_bulk = bulk(symbols, structure, a=a_fcc, cubic=True)
        elif structure == 'hcp':
            my_bulk = bulk(symbols, structure, a = a_hcp, covera = my_covera, cubic=True)
        else:
            raise ValueError('%s is an invalid structure' % structure)
    my_bulk = move_atoms(my_bulk, before_move_atom)

    layer_number_per_slab = my_find_num_per_slab(my_bulk, slice_plane, my_tol, my_periodic, number_per_layer)
    # print('layer_number_per_slab: %s' % layer_number_per_slab)

    if num_layers:
        num_rep_z = int(num_layers/layer_number_per_slab)
    elif replic_z:
        num_rep_z = replic_z
    else:
        raise ValueError('%s or %s is an invalid value' % num_layers % num_rep_z)

    # print('rep_z: %s' %num_rep_z)
    my_slab = surface(my_bulk, slice_plane , num_rep_z, vacuum = my_vacuum, tol=my_tol, periodic=my_periodic)
    print(my_slab)
    my_slab = my_find_prim(my_slab)

    my_slab = move_atoms(my_slab, move_atom)
    my_slab.wrap()
    return my_slab

# stretch unit cell
def stretch_list_along_direction_to_cell(atoms: Atoms = None , stretch_factor_list: list = [0.997, 0.998, 0.999, 1.000, 1.001, 1.002, 1.003],
                                         stretch_direction_list: list = ['x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz'],
                                        stretch_3_direction_lists: list = None, my_scale_atoms: bool = True) -> list:
    """
    Stretches the unit cell of the provided atomic structure along specified directions by a set of stretch factors.
    The function supports stretching in multiple directions and provides an option to scale atomic positions accordingly.

    Args:
        atoms (Atoms): The atomic structure to be stretched, represented as an ASE Atoms object.
        stretch_factor_list (list of float): A list of stretch factors to apply to the unit cell in the specified directions.
        stretch_direction_list (list of str): A list of directions (e.g., 'x', 'y', 'z', 'xy', etc.) corresponding to stretch factors.
        stretch_3_direction_lists (list of list of float, optional): A list of lists of stretch factors for x, y, z directions 
            (e.g., [[0.997, 1.000, 1.000], [1.000, 0.997, 1.000], [1.000, 1.000, 0.997]]). If provided, it overrides `stretch_factor_list` and `stretch_direction_list`.
        my_scale_atoms (bool): Whether to scale the atomic positions along with the unit cell (default is True).

    Returns:
        list of Atoms: A list of atomic structures, each stretched by the corresponding stretch factor in the given direction.

    Raises:
        ValueError: If `atoms` is not provided or not an instance of ASE Atoms.
        ValueError: If the lengths of `stretch_factor_list` and `stretch_direction_list` are not equal when `stretch_3_direction_lists` is not provided.
    
    Example:
        >>> from ase.build import bulk
        >>> atoms = bulk('Cu', 'fcc', a=3.6)
        >>> stretched_atoms = stretch_list_along_direction_to_cell(atoms, stretch_factor_list=[0.99, 1.01], stretch_direction_list=['x', 'y'])
        >>> for stretched in stretched_atoms:
        >>>     print(stretched.get_cell())
    """

    # Input validation
    if atoms is None or not isinstance(atoms, Atoms):
        raise ValueError("`atoms` must be an instance of ASE Atoms.")
    
    if stretch_3_direction_lists is not None:
        if not isinstance(stretch_3_direction_lists, list) or not all(isinstance(direction_list, list) and len(direction_list) == 3 for direction_list in stretch_3_direction_lists):
            raise ValueError("`stretch_3_direction_lists` must be a list of lists, and each list must contain three stretch factors (for x, y, z directions).")
    else:
        if not isinstance(stretch_factor_list, list) or not all(isinstance(factor, (float, int)) for factor in stretch_factor_list):
            raise ValueError("`stretch_factor_list` must be a list of floats or integers representing stretch factors.")
        
        if not isinstance(stretch_direction_list, list) or not all(isinstance(direction, str) for direction in stretch_direction_list):
            raise ValueError("`stretch_direction_list` must be a list of strings representing stretch directions (e.g., 'x', 'y', 'z', 'xy').")

    atoms_copy = atoms.copy()
    # storing atoms after being stretched using list
    stretched_atoms = []
    if stretch_3_direction_lists:
        for stretch_3_direction_list in stretch_3_direction_lists:
            stretch_along_direction_to_cell(atoms_copy, stretch_3_direction_list = stretch_3_direction_list, my_scale_atoms = my_scale_atoms, stretched_atoms_list = stretched_atoms)
    elif len(stretch_factor_list) == len(stretch_direction_list):
        for index, stretch_factor in enumerate(stretch_factor_list):
            stretch_direction = stretch_direction_list[index]
            stretch_along_direction_to_cell(atoms_copy, stretch_factor, stretch_direction, my_scale_atoms = my_scale_atoms, stretched_atoms_list = stretched_atoms)
    elif len(stretch_direction_list) == 1:
        stretch_direction = stretch_direction_list[0]
        for stretch_factor in stretch_factor_list:
            stretch_along_direction_to_cell(atoms_copy, stretch_factor, stretch_direction, my_scale_atoms = my_scale_atoms, stretched_atoms_list = stretched_atoms)
    return stretched_atoms

def stretch_along_direction_to_cell(atoms: Atoms = None ,stretch_factor: float = None, stretch_direction: chr = None, stretch_3_direction_list: list = None,
                                my_scale_atoms: bool = True, stretched_atoms_list: list = None):
    """
    Stretches the unit cell of the provided atomic structure along specified directions or by specific stretch factors for each axis. 
    Allows for scaling of atomic positions along with the unit cell.

    Args:
        atoms (Atoms): The atomic structure to be stretched, represented as an ASE Atoms object.
        stretch_factor (float, optional): The factor by which to stretch the unit cell in the specified direction.
        stretch_direction (str, optional): The direction along which to apply the stretch (e.g., 'x', 'y', 'z', 'xy', 'xyz').
        stretch_3_direction_list (list of float, optional): A list of three stretch factors to apply along the x, y, and z axes. 
            If provided, this will override `stretch_factor` and `stretch_direction`.
        my_scale_atoms (bool, optional): Whether to scale the atomic positions along with the unit cell (default is True).
        stretched_atoms_list (list of Atoms, optional): A list to store the stretched atomic structure. If provided, the stretched
            structure will be appended to this list.

    Returns:
        Atoms: The atomic structure with the stretched unit cell.

    Raises:
        ValueError: If `atoms` is not provided or is not an instance of ASE Atoms.
        ValueError: If neither `stretch_factor` nor `stretch_3_direction_list` is provided.
        ValueError: If `stretch_3_direction_list` is provided but is not a list of three elements.
        ValueError: If `stretch_direction` is provided and is not a valid direction (e.g., 'x', 'y', 'z', 'xy', 'xyz').

    Example:
        >>> atoms = bulk('Cu', 'fcc', a=3.6)
        >>> stretched_atoms = stretch_along_direction_to_cell(atoms, stretch_factor=1.01, stretch_direction='x')
        >>> print(stretched_atoms.get_cell())
    """

    # Input validation
    if atoms is None or not isinstance(atoms, Atoms):
        raise ValueError("`atoms` must be an instance of ASE Atoms.")
    
    if stretch_3_direction_list is None:
        if stretch_factor is None or stretch_direction is None:
            raise ValueError("Either `stretch_3_direction_list` or both `stretch_factor` and `stretch_direction` must be provided.")
    else:
        if not isinstance(stretch_3_direction_list, list) or len(stretch_3_direction_list) != 3:
            raise ValueError("`stretch_3_direction_list` must be a list of three elements representing stretch factors for x, y, and z axes.")
    
    if stretch_3_direction_list is None and stretch_direction not in ['x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz', 'yx', 'zy', 'zx', 'xzy', 'yzx', 'yxz', 'zxy', 'zyx']:
        raise ValueError("`stretch_direction` must be one of the following: 'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz', or combinations like 'xzy'.")

    atoms_copy = atoms.copy()
    my_cell = array(atoms_copy.get_cell())
    temp = my_cell.copy()
    if stretch_3_direction_list and len(stretch_3_direction_list) == 3:
        for i, factor in enumerate(stretch_3_direction_list):
            temp[:,i] = my_cell[:,i] * factor
    else:
        if stretch_direction == 'x':
            temp[:,0] = my_cell[:,0] * stretch_factor
        elif stretch_direction == 'y':
            temp[:,1] = my_cell[:,1] * stretch_factor
        elif stretch_direction == 'z':
            temp[:,2] = my_cell[:,2] * stretch_factor
        elif stretch_direction == 'xy' or 'yx':
            temp[:,:2] = my_cell[:,:2] * stretch_factor
        elif stretch_direction == 'yz' or 'zy':
            temp[:,1:] = my_cell[:,1:] * stretch_factor
        elif stretch_direction == 'zx' or 'xz':
            temp[:,2] = my_cell[:,2] * stretch_factor
            temp[:,0] = my_cell[:,0] * stretch_factor
        elif stretch_direction == 'xyz' or 'xzy' or 'yzx' or 'yxz' or 'zxy' or 'zyx':
            temp[:,2] = my_cell[:,2] * stretch_factor
    atoms_copy.set_cell(temp, scale_atoms=my_scale_atoms)
    if stretched_atoms_list is not None:
        stretched_atoms_list.append(atoms_copy)
    return atoms_copy

def print_after_what(char_name = '', variable = None,calling_function = '', specified_blank = '', character = '||', after_what = 'none'):
    # print(f"the value is None! - {calling_function}")
    print(f"the {character} {char_name} {variable} {character} is {after_what}! - {calling_function}")
    return specified_blank

# move atoms
def move_atoms(atoms: Atoms = None,
               translate_matrix: ndarray = array([0.1, 0.1, 0.0]),
               if_scale_position: bool = True) -> Atoms :
    """
    Translates the atomic positions in the provided structure by a given translation matrix. 
    The translation can be applied either in scaled (fractional) coordinates or in Cartesian coordinates.

    Args:
        atoms (Atoms): The atomic structure to be translated, represented as an ASE Atoms object.
        translate_matrix (ndarray): The translation vector to apply to each atom, either in fractional or Cartesian coordinates, depending on `if_scale_position`. 
            The default is [0.1, 0.1, 0.0].
        if_scale_position (bool): Whether to apply the translation in scaled (fractional) coordinates. If `True`, the translation is applied to the scaled 
            positions (default). If `False`, the translation is applied to the Cartesian positions.

    Returns:
        Atoms: The translated atomic structure.

    Raises:
        ValueError: If the `atoms` argument is not provided or is not an instance of ASE Atoms.
        ValueError: If `translate_matrix` is not a 1D array of length 3.
        ValueError: If `if_scale_position` is not a boolean.

    Example:
        >>> from ase.build import molecule
        >>> atoms = molecule('H2O')
        >>> translated_atoms = move_atoms(atoms, translate_matrix=np.array([0.5, 0.5, 0.0]), if_scale_position=False)
        >>> print(translated_atoms.get_positions())
    """
    translate_matrix = array(translate_matrix)
    # Input validation
    if atoms is None or not isinstance(atoms, Atoms):
        raise ValueError("`atoms` must be an instance of ASE Atoms.")
    
    if not isinstance(translate_matrix, ndarray) or translate_matrix.shape != (3,):
        raise ValueError("`translate_matrix` must be a 1D numpy array of length 3.")
    
    if not isinstance(if_scale_position, bool):
        raise ValueError("`if_scale_position` must be a boolean value.")
    
    atoms_copy = atoms.copy()
    if if_scale_position:
        scaled = atoms_copy.get_scaled_positions()
        scaled += translate_matrix
        atoms_copy.set_scaled_positions(scaled)
    else:
        positions = atoms_copy.get_positions()
        positions += translate_matrix
        atoms_copy.set_positions(positions)
    atoms_copy.wrap()
    return atoms_copy

# useless
def file_name(dir_path: str = None, base_name: str = None, special_name: str = None, format_type: str = None) -> str:
    formatted_factor = format_type % special_name
    filename = f'{dir_path}{base_name}{special_name}'
    return file_name


# def my_find_prim(atoms: Atoms = None) -> Atoms:
#     """
#     find primitive cell using spglib\n
#     Convert to a format suitable for spglib
#     """
    
#     lattice = array(atoms.get_cell())
#     points = array(atoms.get_scaled_positions())
#     numbers = array(atoms.get_atomic_numbers())
#     pbc = array(atoms.get_pbc())
#     cell = (lattice, points, numbers)

#     primitive_cell = spglib.standardize_cell(cell, to_primitive=1, no_idealize=1)
#     # Convert the spglib output back to an ASE Atoms object
#     primitive_atoms = Atoms(numbers = primitive_cell[2],
#                             scaled_positions = primitive_cell[1],
#                             cell = primitive_cell[0],
#                             pbc=pbc)
#     return primitive_atoms

# find number of atoms per slab
def my_find_num_per_slab(my_bulk: Atoms = None,
                         slice_plane:  tuple[int] = (0, 0, 1),
                         my_tol: float = 1e-6,          
                         my_periodic: bool = False,
                         number_per_layer: float = None) -> float:
    my_one_slab = surface(my_bulk, slice_plane , 1, vacuum = 20, tol=my_tol, periodic=my_periodic)
    prim_one_slab = my_find_prim(my_one_slab)
    atom_number_per_slab = len(prim_one_slab.get_positions())
    layer_number_per_slab = atom_number_per_slab/number_per_layer
    return layer_number_per_slab
    

# # usage
# def uni_axial_stretch():
#     film = generate_film(symbols = 'Au', structure = 'fcc', num_layers = 12, my_vacuum = 20, slice_plane = (1,1,1), a_fcc = 2.95*sqrt(2.0))
#     stretch_factor_list = [0.997 + i * 0.001 for i in range(7)]
#     #[0.997, 0.998, 0.999, 1.000, 1.001, 1.002, 1.003]
#     films_stretch = stretch_list_along_direction_to_cell(film , stretch_factor_list = stretch_factor_list, stretch_direction_list = ['x'])
    
#     format_type = '%.3f'
#     for i, film_stretch in enumerate(films_stretch):
#         formatted_i = format_type % stretch_factor_list[i]
#         #print(formatted_i)
#         #print(array(film_stretch.get_cell()))
#         filename = f'./y_dir/{formatted_i}/POSCAR' 
#         makedirs(path.dirname(filename), exist_ok=True)   
#         my_write_vasp(filename, film_stretch, label = f'Au thin film {formatted_i}')
        
# uni_axial_stretch()

