from numpy import sqrt, ndarray, array, int32
from os    import makedirs, path
from inspect import stack
from ase.cell import Cell
from ase.symbols import Symbols
import numpy as np

from ase import Atoms
from ase.data import chemical_symbols
from ase.build import bulk, surface
from ase.io import write
from ase.visualize import view
# For version independency, it could be instead by 'from ase.io.vasp import *'
from ase.utils import reader, writer

import re
from pathlib import Path
from typing import List, Optional, TextIO, Tuple
import spglib

import sys
print(sys.path)

from mymetal.io.vasp import my_write_vasp

def generate_film(   symbols: str = None,                 # str
                structure: str = None,              # str
                num_layers: int = None,            # int
                replic_z: int = None,              # int
                my_vacuum: float = None,           # float
                slice_plane: tuple[int] = (0, 0, 1),  # Tuple of ints
                my_tol: float = 1e-6,              # float
                my_periodic: bool = False,          # bool
                a_fcc: float = 2.95*sqrt(2.0),     # float
                a_hcp: float = 2.95,               # float
                my_covera: float = sqrt(8.0/3.0),  # float
                move_atom: list = [0.1, 0.1, 0.0],
                number_per_layer: float = 1.0
                ) -> Atoms:
    # parameters : 1-8 line: general setting
    #              9 10-11 line: fcc, hcp parameters
    # when we use surface() function, my_bulk must be a conventional cell, not a primitive cell, so set cubic=True

    if structure == 'fcc':
        my_bulk = bulk(symbols, structure, a=a_fcc, cubic=True)
    elif structure == 'hcp':
        my_bulk = bulk(symbols, structure, a = a_hcp, covera = my_covera, cubic=True)
    else:
        raise ValueError('%s is an invalid structure' % structure)
    
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
    my_slab = my_find_prim(my_slab)

    my_slab = move_atoms(my_slab, move_atom)
    my_slab.wrap()
    return my_slab

# stretch unit cell
def stretch_list_along_direction_to_cell(atoms: Atoms = None , stretch_factor_list: list = [0.997, 0.998, 0.999, 1.000, 1.001, 1.002, 1.003],
                                         stretch_direction_list: list = ['x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz'],
                                        stretch_3_direction_lists: list = None, my_scale_atoms: bool = True):
    # stretch again loop
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
               translate_matrix: ndarray = array([0.1, 0.1, 0.0])) -> Atoms :
    scaled = atoms.get_scaled_positions()
    scaled += translate_matrix
    atoms.set_scaled_positions(scaled)
    atoms.wrap()
    return atoms


# useless
def file_name(dir_path: str = None, base_name: str = None, special_name: str = None, format_type: str = None) -> str:
    formatted_factor = format_type % special_name
    filename = f'{dir_path}{base_name}{special_name}'
    return file_name


def my_find_prim(atoms: Atoms = None) -> Atoms:
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
    return primitive_atoms

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
    

# usage
def uni_axial_stretch():
    film = generate_film(symbols = 'Au', structure = 'fcc', num_layers = 12, my_vacuum = 20, slice_plane = (1,1,1), a_fcc = 2.95*sqrt(2.0))
    stretch_factor_list = [0.997 + i * 0.001 for i in range(7)]
    #[0.997, 0.998, 0.999, 1.000, 1.001, 1.002, 1.003]
    films_stretch = stretch_list_along_direction_to_cell(film , stretch_factor_list = stretch_factor_list, stretch_direction_list = ['x'])
    
    format_type = '%.3f'
    for i, film_stretch in enumerate(films_stretch):
        formatted_i = format_type % stretch_factor_list[i]
        #print(formatted_i)
        #print(array(film_stretch.get_cell()))
        filename = f'./y_dir/{formatted_i}/POSCAR' 
        makedirs(path.dirname(filename), exist_ok=True)   
        my_write_vasp(filename, film_stretch, label = f'Au thin film {formatted_i}')
        
uni_axial_stretch()

