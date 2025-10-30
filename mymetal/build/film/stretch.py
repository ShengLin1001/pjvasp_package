"""
stretch submodule

This submodule provides functions for stretching the unit cell of a thin film structure along specified directions.

Functions:
    - generate_film: Generates a thin film structure (slab) from bulk material and applies modifications.
    - stretch_list_along_direction_to_cell: Stretches the unit cell of the provided atomic structure along specified directions by a set of stretch factors.
    - stretch_along_direction_to_cell: Stretches the unit cell of the provided atomic structure along specified directions or by specific stretch factors for each axis.
    - file_name: Generates a file name based on the provided directory path, base name, special name, and format type.
    - find_num_per_slab: Finds the number of atoms per slab in a thin film structure.
    - adjust_lattice_for_volume_conservation: Adjusts the lattice vectors to conserve volume after stretching.
    - generate_bulk_from_film: Generates a bulk structure from a given thin film structure.
"""
# For version independency, it could be instead by 'from ase.io.vasp import *'
# The strict has been removed.
import spglib

import os

import numpy as np
from numpy import sqrt, ndarray, array

from ase import Atoms
from ase.build import bulk, surface

from mymetal.build.film.findprim import my_find_prim
from mymetal.build.film.extrfilm import my_extr_thick
from mymetal.universal.atom.moveatom import move_atoms

####################################################################################################
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
####################################################################################################
# # Example usage:
# # Lattice vectors before stretching
# lattice_before = np.array([[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]])

# # Lattice vectors after stretching
# lattice_after = np.array([[2.5, 0.0, 0.0], [0.0, 2.5, 0.0], [0.0, 0.0, 2.5]])

# # Adjust lattice to conserve volume
# adjusted_lattice = adjust_lattice_for_volume_conservation(lattice_before, lattice_after, change_axis=[0])
# print(adjusted_lattice)
####################################################################################################

def generate_film(   symbols: str = None,                 # str
                structure: str = None,              # str
                num_layers: int = None,            # int
                replic_z: int = None,              # int
                my_vacuum: float = 7.5,           # float
                slice_plane: tuple[int] = (0, 0, 1),  # Tuple of ints
                my_tol: float = 1e-6,              # float
                my_periodic: bool = True,          # bool
                a_fcc: float = 2.95*sqrt(2.0),     # float
                a_hcp: float = 2.95,               # float
                my_covera: float = sqrt(8.0/3.0),  # float
                move_atom: list = [1e-6, 1e-6, 1e-6],
                before_move_atom: list = [1e-6, 1e-6, 1e-6],
                number_per_layer: float = 1.0,
                bulk_atoms: Atoms = None,
                if_find_prim: bool = True,
                if_fix_c: bool = False
                ) -> Atoms:
    """Generate a slab OR specified (ijk) plane-based bulk from a conventional ideal bulk crystal (supports FCC/HCP or a provided bulk).

    This utility wraps ASE-like operations (e.g. `bulk`, `surface`, `cut`) and custom
    helpers (`my_find_prim`, `my_find_num_per_slab`, `move_atoms`) to produce a slab
    oriented along a Miller plane, optionally with vacuum padding and atom translations.

    The function supports two main modes:
      * Build a conventional bulk cell internally from `symbols` + `structure` ('fcc' or 'hcp'),
        then slice it.
      * Use an externally supplied `bulk_atoms` ASE Atoms object.
    It also supports generating a bulk cell that is aligned with an arbitrary (h,k,l)
    plane by setting `vacuum=None`, `my_periodic=True`, `if_fix_c=True` (i.e., produce a periodic
    bulk rather than a non-periodic slab with vacuum).

    Args:
        symbols (str, optional):
            Element symbol for internally-generated bulk (e.g. "Cu"). Only single-element
            systems are supported when generating bulk internally. Required if `bulk_atoms`
            is not provided.
        structure (str, optional):
            Crystal type for internal bulk generation. One of {'fcc', 'hcp'}. Required if
            `bulk_atoms` is None.
        num_layers (int, optional):
            Desired number of atomic layers in the final slab. If provided, it takes precedence
            over `replic_z`. The function determines the number of layers per primitive repeat
            using `my_find_num_per_slab()` and computes the required z-repetitions.
        replic_z (int, optional):
            Direct number of repetitions along the surface normal (z-like direction). Used
            if `num_layers` is not provided.
        my_vacuum (float, optional):
            Vacuum thickness (Å) to add above and below the slab. Default 7.5 Å.
            If `my_vacuum` is `None`, no vacuum is added. Note: the implementation adds
            *two* vacuum regions (top and bottom) equal to this value when creating non-periodic slabs.
        slice_plane (tuple[int], optional):
            Miller indices (h, k, l) defining the slicing plane / surface normal. Default (0,0,1).
        my_tol (float, optional):
            Numerical tolerance used for detecting layers and geometry checks (units: Å). Default 1e-6.
        my_periodic (bool, optional):
            If True (default), the produced slab/cell will keep periodic boundary conditions along the surface
            normal (useful to produce a bulk-like (h,k,l)-based cell). If False, the slab
            will be non-periodic along the surface normal with vacuum added (unless `my_vacuum is None`).
            Be carefully when setting to False!
        a_fcc (float, optional):
            Lattice constant used when `structure == 'fcc'` and `bulk_atoms` not provided.
            Default uses `2.95 * sqrt(2.0)` to form a conventional FCC cell.
        a_hcp (float, optional):
            Lattice `a` constant used when `structure == 'hcp'` and `bulk_atoms` not provided.
        my_covera (float, optional):
            c/a ratio for HCP crystals when building a bulk cell internally.
        move_atom (list[float], optional):
            Translation applied *after* slab generation. Given as fractional translation
            in lattice (or internal coordinate) units (x, y, z). Default [1e-6, 1e-6, 1e-6].
        before_move_atom (list[float], optional):
            Translation applied *before* any slab creation (i.e., applied to the initial bulk).
            Also fractional in lattice units. Default [1e-6, 1e-6, 1e-6].
            Purpose: break symmetries / avoid atoms lying exactly on cutting planes.
        number_per_layer (float, optional):
            Expected number of atoms per layer used by `my_find_num_per_slab()` to determine
            how many layers constitute one repeat. Default 1.0.
        bulk_atoms (ase.Atoms, optional):
            If provided, this Atoms object is used as the bulk from which to slice. The
            function will call `my_find_prim(bulk_atoms, to_primitive=0)` to make a suitable
            conventional/working cell before slicing.
        if_find_prim (bool, optional):
            If True (default), the resulting slab will be reduced to a primitive cell via
            `my_find_prim()` after slab creation.
        if_fix_c (bool, optional):
            If True, when reducing to primitive cell, the c-axis and zpositions (surface normal) 
            will be preserved/fixed. Default False.

    Returns:
        ase.Atoms:
            The final slab / thin-film as an ASE Atoms object. The returned Atoms will have
            been optionally reduced to primitive cell, translated by `move_atom`, and wrapped
            (`Atoms.wrap()`) before return.

    Raises:
        ValueError:
            * If `structure` is not 'fcc' or 'hcp' while `bulk_atoms` is None.
            * If neither `num_layers` nor `replic_z` is provided.
            The raised message contains the invalid argument(s).

    Notes:
        * Layer counting: `my_find_num_per_slab(my_bulk, slice_plane, my_tol, my_periodic, number_per_layer)`
          is used to determine how many atomic layers appear per slab repeat. `num_layers` is converted
          to `replic_z` with `num_rep_z = int(num_layers / layer_number_per_slab)`. If this division does not
          produce an integer, the integer cast truncates — ensure `num_layers` is compatible with the
          slab repeat or supply `replic_z` directly.
        * `vacuum=None` + `my_periodic=True` + `if_fix_c=True`:
            This produces a periodic (h,k,l)-based bulk cell (no vacuum added) suitable for calculations
            where periodicity along the surface-normal is required (i.e., not a slab).
        * `move_atom` and `before_move_atom` are applied in fractional/lattice units (not Å). They are used
          to avoid atoms sitting exactly on the slicing plane or to intentionally displace atoms.
        * The function expects helper utilities to exist in the same codebase:
            - `my_find_prim(atoms, to_primitive=0)`
            - `my_find_num_per_slab(...)`
            - `move_atoms(atoms, translation_list)`
            - ASE functions: `bulk()`, `surface()`, and `Atoms.wrap()`.
        * If `bulk_atoms` is already primitive vs conventional, `my_find_prim(..., to_primitive=0)` is called
          to ensure the correct conventional/predictable cell prior to slicing (except for HCP where cubic
          flag is not applicable).
        * It's importan to use `move_atom` and `before_move_atom` to avoid atoms exactly on the cutting plane,
            which can lead to unpredictable behavior.

    Examples:
        # 1) Generate a 10-layer FCC Cu slab with vacuum (default behavior)
        >>> slab = generate_film(symbols='Cu', structure='fcc', num_layers=10, my_vacuum=10.0)

        # 2) Use an externally constructed bulk and make a (1,1,1) periodic cell (no vacuum)
        >>> bulk_atoms = bulk('Au', 'fcc', a=4.08, cubic=True)
        >>> surface_atoms = generate_film(bulk_atoms = bulk_atoms, slice_plane = (1, 1, 1), num_layers = 12,
        ...                   my_vacuum=None, my_periodic = True, if_find_prim=True, if_fix_c=True,
        ...                   before_move_atom = [1e-6, 1e-6, 1e-6], move_atom = [1e-6, 1e-6, 1e-6])

        # 3) Shift atoms before cutting to avoid atoms on the cutting plane,
        #    then apply a small post-generation shift and reduce to primitive
        >>> slab = generate_film(symbols='Al', structure='fcc',
        ...                      num_layers=6,
        ...                      before_move_atom=[0.01, 0.01, 0.02],
        ...                      move_atom=[0.05, 0.05, 0.0],
        ...                      if_find_prim=True)

    See Also:
        ASE: `ase.build.bulk`, `ase.build.surface`, `ase.geometry.cut`
        Local helpers: `my_find_prim`, `my_find_num_per_slab`, `move_atoms`

    Implementation detail:
        * The function will:
            1. Build or accept a bulk (`my_bulk`).
            2. Apply `before_move_atom` via `move_atoms(my_bulk, before_move_atom)`.
            3. Compute `layer_number_per_slab = my_find_num_per_slab(...)`.
            4. Convert `num_layers` -> `num_rep_z` or use `replic_z`.
            5. Call `surface(my_bulk, slice_plane, num_rep_z, vacuum=my_vacuum, tol=my_tol, periodic=my_periodic)`.
            6. Optionally call `my_find_prim(my_slab)`.
            7. Apply `move_atoms(my_slab, move_atom)`, call `my_slab.wrap()`, and return.

    """
    if bulk_atoms:
        my_bulk = my_find_prim(bulk_atoms, to_primitive=0)
    else:
        if structure == 'fcc':
            my_bulk = bulk(symbols, structure, a=a_fcc, cubic=True)
        elif structure == 'hcp':
            my_bulk = bulk(symbols, structure, a = a_hcp, covera = my_covera)
        else:
            raise ValueError('%s is an invalid structure' % structure)
    my_bulk = move_atoms(my_bulk, before_move_atom)

    print('number_per_layer: %s' % number_per_layer)
    layer_number_per_slab = my_find_num_per_slab(my_bulk, slice_plane, my_tol, my_periodic, number_per_layer)

    if num_layers:
        num_rep_z = int(num_layers/layer_number_per_slab)
    elif replic_z:
        num_rep_z = replic_z
    else:
        raise ValueError('%s or %s is an invalid value' % num_layers % num_rep_z)

    print('num_rep_z: %s' %num_rep_z)
    my_slab = surface(my_bulk, slice_plane , num_rep_z, vacuum = my_vacuum, tol=my_tol, periodic=my_periodic)
    # print(my_slab)

    if if_find_prim:
        my_slab = my_find_prim(my_slab, check_direction_tag=False, if_fix_c = if_fix_c)
    print('len(my_slab): %s' % len(my_slab))

    my_slab = move_atoms(my_slab, move_atom)
    my_slab.wrap()
    return my_slab

# stretch unit cell
def stretch_list_along_direction_to_cell(atoms: Atoms = None , stretch_factor_list: list = [0.997, 0.998, 0.999, 1.000, 1.001, 1.002, 1.003],
                                         stretch_direction_list: list = ['x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz', '1', '2', '3', '12', '123'],
                                        stretch_3_direction_lists: list = None, my_scale_atoms: bool = True, stretch_3_direction_symbol: str = 'xyz') -> list:
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
            stretch_along_direction_to_cell(atoms_copy, stretch_3_direction_list = stretch_3_direction_list, my_scale_atoms = my_scale_atoms, stretched_atoms_list = stretched_atoms,
                                            stretch_3_direction_symbol = stretch_3_direction_symbol)
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
                                my_scale_atoms: bool = True, stretched_atoms_list: list = None, stretch_3_direction_symbol: str = 'xyz'):
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
            if stretch_3_direction_symbol == 'xyz':
                temp[:,i] = my_cell[:,i] * factor
            elif stretch_3_direction_symbol == 'a1a2a3':
                temp[i,:] = my_cell[i,:] * factor
    else:
        if stretch_direction == 'x':
            temp[:, 0] = my_cell[:, 0] * stretch_factor
        elif stretch_direction == 'y':
            temp[:, 1] = my_cell[:, 1] * stretch_factor
        elif stretch_direction == 'z':
            temp[:, 2] = my_cell[:, 2] * stretch_factor
        elif stretch_direction in ['xy', 'yx']:
            temp[:, :2] = my_cell[:, :2] * stretch_factor
        elif stretch_direction in ['yz', 'zy']:
            temp[:, 1:] = my_cell[:, 1:] * stretch_factor
        elif stretch_direction in ['xz', 'zx']:
            temp[:, [0,2]] = my_cell[:, [0,2]] * stretch_factor
        elif stretch_direction in ['xyz', 'xzy', 'yzx', 'yxz', 'zxy', 'zyx']:
            temp[:, :] = my_cell[:, :] * stretch_factor
        elif stretch_direction == '1':
            temp[0, :] = my_cell[0, :] * stretch_factor
        elif stretch_direction == '2':
            temp[1, :] = my_cell[1, :] * stretch_factor
        elif stretch_direction == '3':
            temp[2, :] = my_cell[2, :] * stretch_factor
        elif stretch_direction in ['12', '21']:
            temp[:2, :] = my_cell[:2, :] * stretch_factor
        elif stretch_direction in ['23', '32']:
            temp[1:, :] = my_cell[1:, :] * stretch_factor
        elif stretch_direction in ['13', '31']:
            temp[[0,2], :] = my_cell[[0,2], :] * stretch_factor
        elif stretch_direction in ['123', '132', '213', '231', '312', '321']:
            temp[:, :] = my_cell[:, :] * stretch_factor
        else:
            raise ValueError(f"Unknown stretch_direction: {stretch_direction}")
    atoms_copy.set_cell(temp, scale_atoms=my_scale_atoms)
    if stretched_atoms_list is not None:
        stretched_atoms_list.append(atoms_copy)
    return atoms_copy

def print_after_what(char_name = '', variable = None,calling_function = '', specified_blank = '', character = '||', after_what = 'none'):
    # print(f"the value is None! - {calling_function}")
    print(f"the {character} {char_name} {variable} {character} is {after_what}! - {calling_function}")
    return specified_blank

# useless
def file_name(dir_path: str = None, base_name: str = None, special_name: str = None, format_type: str = None) -> str:
    formatted_factor = format_type % special_name
    filename = f'{dir_path}{base_name}{special_name}'
    return file_name


# find number of atoms per slab
def my_find_num_per_slab(my_bulk: Atoms = None,
                         slice_plane:  tuple[int] = (0, 0, 1),
                         my_tol: float = 1e-6,          
                         my_periodic: bool = False,
                         number_per_layer: float = 1) -> float:
    my_one_slab = surface(my_bulk, slice_plane , 1, vacuum = 20, tol=my_tol, periodic=my_periodic)
    prim_one_slab = my_find_prim(my_one_slab)
    #print('prim_one_slab: %s' % prim_one_slab)
    atom_number_per_slab = len(prim_one_slab.get_positions())
    print('atom_number_per_slab: %s' % atom_number_per_slab)
    layer_number_per_slab = atom_number_per_slab/number_per_layer
    print('layer_number_per_slab: %s' % layer_number_per_slab)
    return layer_number_per_slab


def adjust_lattice_for_volume_conservation(lattice_before, lattice_after, change_axis: list = [0, 1, 2]):
    """
    Adjust the 'lattice_after' so that the volume remains the same as 'lattice_before'.
    
    Parameters:
    lattice_before: np.array
        3x3 array representing the lattice vectors before stretching.
    lattice_after: np.array
        3x3 array representing the lattice vectors after stretching.
        
    Returns:
    adjusted_lattice_after: np.array
        Adjusted 3x3 array representing the lattice vectors after stretching, 
        ensuring volume conservation.
    """
    
    # Compute the volume before and after stretching
    volume_before = np.abs(np.linalg.det(lattice_before))
    volume_after = np.abs(np.linalg.det(lattice_after))
    
    # Scaling factor to adjust the lattice_after to conserve volume
    scaling_factor = None
    if len(change_axis) == 3:
        scaling_factor = (volume_before / volume_after) ** (1/3)
    elif len(change_axis) == 2:
        scaling_factor = (volume_before / volume_after) ** (1/2)
    elif len(change_axis) == 1:
        scaling_factor = volume_before / volume_after
    
    for axis in change_axis:
        lattice_after[axis, :] = lattice_after[axis, :] * scaling_factor

    adjusted_lattice_after = lattice_after
    
    return adjusted_lattice_after, scaling_factor


def generate_bulk_from_film(film: Atoms=None, if_find_prim: bool = False, vacuum: float=None, number_per_layer: int = 1) -> Atoms:
    """
    [Deprecated] Convert a thin film back into an approximate bulk structure.

    ⚠️ **Deprecated:** This function has been replaced by the more general
    `generate_film()` function with the setting:
        `vacuum=None` and `my_periodic=True`.

    The older approach attempted to reconstruct a bulk by removing the vacuum
    regions from a given slab. It is kept for backward compatibility and testing
    but is no longer recommended for production use.

    Args:
        film (Atoms, optional): The Atoms object representing the thin film. Defaults to None.
        if_find_prim (bool, optional): Whether to find the primitive cell of the bulk structure. Defaults to False.

    Returns:
        Atoms: A new Atoms object representing the bulk structure.

    Raises:
        ValueError: If the number of layers in the film is less than or equal to 1.
        
    Example:
        >>> bulk_atoms = generate_bulk_from_film(film_atoms)
    """
    atoms = film.copy()
    numatom = len(atoms)
    numlayer = numatom / number_per_layer
    if numlayer <= 1:
        raise ValueError("The number of layers in the film must be greater than 1.")
    thick = my_extr_thick(atoms)
    if vacuum is not None:
        vacuum = vacuum
    else:
        vacuum = thick/(numlayer-1)/2
    atoms.center(vacuum=vacuum, axis=2)

    if if_find_prim:
        atoms = my_find_prim(atoms)

    return atoms
