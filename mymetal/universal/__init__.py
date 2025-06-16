"""
mymetal.universal

This subpackage provides a collection of utility functions and tools that are commonly used across various 
modules in the `mymetal` package. It includes functions for input validation, data adjustments, atom manipulation, 
data transformations, matrix adjustments, and more. These utilities are designed to streamline common tasks 
in materials science simulations and data handling.

Modules:
    - check: Module for input validation and error checking.
    - data: Module for data adjustments and transformations.
    - atom: Module for atom manipulation and analysis.
    - search: Module for searching and extracting data from files.
    - print: Module for printing messages and warnings.
    - matrix: Module for matrix adjustments and transformations.
    - index: Module for indexing transformations and adjustments.
    - plot: Module for plotting data and results.
"""

from mymetal.universal.check.checkinput import check_input
from mymetal.universal.data.dataadjust import rm_blank, my_add_list, my_down_up, myjust, normalize_float_int
from mymetal.universal.atom.moveatom import move_atoms
from mymetal.universal.data.datachange import list_to_char, char_to_dic, dic_to_char
from mymetal.universal.search.find import find_line_position, extract_line_at_position
from mymetal.universal.data.patterntrans import my_pattern_trans
from mymetal.universal.print.printafter import print_after_blank, print_after_cant_read, print_after_read, print_after_not_supported
from mymetal.universal.atom.delatom import mydel_pos_type, check_position
from mymetal.universal.atom.fixatom import fixatoms
from mymetal.universal.matrix.adjust import adjust_matrix
from mymetal.universal.index.indextrans import three_index_to_four_index
from mymetal.universal.plot.plot import my_plot, my_plot_brokenaxed, my_plot_energy_components, my_plot_interlayer_distance, my_plot_zpositions
from mymetal.universal.plot.oldplotdos import my_plot_orientation, my_plot_complete_dos, my_plot_horizontal_vertical, my_plot_complete_dos, my_plot_idos
from mymetal.universal.math.operations import get_integration


__all__ = ['check_input',
           'rm_blank', 'my_add_list', 'my_down_up', 'myjust', 'normalize_float_int',
              'move_atoms',
                'list_to_char', 'char_to_dic', 'dic_to_char',
                    'find_line_position', 'extract_line_at_position',
                    'my_pattern_trans',
                    'print_after_blank', 'print_after_cant_read', 'print_after_read', 'print_after_not_supported',
                    'mydel_pos_type', 'check_position',
                    'fixatoms',
                    'adjust_matrix',
                    'three_index_to_four_index',
                    'my_plot', 'my_plot_brokenaxed', 'my_plot_energy_components', 'my_plot_interlayer_distance', 'my_plot_zpositions'      ,
                    'my_plot_orientation', 'my_plot_complete_dos', 'my_plot_horizontal_vertical', 'my_plot_complete_dos', 'my_plot_idos',
                    'get_integration'
]