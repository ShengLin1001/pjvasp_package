"""
mymetal.build

This subpackage provides functions for building structures, including
creating surfaces, adding atoms, and stretching structures. The subpackage
is organized into several modules that focus on different aspects of
building structures.

Modules:
    - film: Contains functions for building thin films, including adding
        hydroxyl groups, stretching films, and finding the primitive cell.
"""

from mymetal.build.film.findprim import (my_find_prim, move_atoms, check_direction)
from mymetal.build.film.hydroxyl import (add_hydro_atoms, get_neighbors, find_matching_atom_in_bulk, 
                                         find_unsaturated_neighbors, passivate_surface_custom)
from mymetal.build.film.stretch import (generate_film, stretch_list_along_direction_to_cell)
#from mymetal.build.film.findhetero import (find_hetero, my_plot_results, set_length, magnitude, build_supercells, 
#                                      stack_atoms, scale_cell_xy, split_model)
from mymetal.build.film.extrfilm import (my_extr_lattice, my_extr_thick, my_extr_etot, cal_area)
from mymetal.build.film.findcubic import (find_cubic)


__all__ = ['my_find_prim', 'move_atoms', 'check_direction',
           'add_hydro_atoms', 'get_neighbors', 'find_matching_atom_in_bulk', 'find_unsaturated_neighbors', 'passivate_surface_custom',
            'generate_film', 'stretch_list_along_direction_to_cell',
#            'find_hetero', 'my_plot_results', 'set_length', 'magnitude', 'build_supercells', 'stack_atoms', 'scale_cell_xy', 'split_model',
#            'my_extr_lattice', 'my_extr_thick', 'my_extr_etot', 'cal_area',
            'find_cubic',
]
