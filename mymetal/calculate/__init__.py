"""
mymetal.calculate

This subpackage provides various calculation functions for material analysis, 
including surface energy calculations, mismatch analysis, deformation, strain 
analysis, and quantum mechanics (QM) calculations. The subpackage is organized 
into several modules that focus on different types of calculations.

Modules:
    - calenergy: Contains functions for calculating energy-related properties, 
      such as surface energy.
    - calmismatch: Contains functions for comparing structures, calculating 
      atom numbers, mismatches, and lattice stretching.
    - calmechanics: Includes functions for mechanical deformation and strain 
      analysis, such as calculating strain matrices, von Mises strain, and 
      deformation matrices.
    - calmath: Provides matrix-related operations, including calculating the 
      Hermite normal form.
    - calqm: Contains functions for quantum mechanics calculations, such as 
      generating k-points for QM calculations.
    - electronic_structure: Contains functions for plotting electronic structure
      data and other universal electronic structure-related functions.
    - material_science: Contains functions for material science calculations,
"""


#from mymetal.calculate.calmismatch.calhetero import (compare_atoms, cal_atom_num, relative_diff, cal_mismatch, filter_results, cal_stretch_lattice)
#from mymetal.calculate.calmechanics.deformation import (cal_deform_matrix)
#from mymetal.calculate.calmechanics.strain import (cal_principal_and_shear_strain, cal_principal_and_shear_strain_root, 
#                                                   cal_strain_matrix, cal_strain_matrix_root, cal_von_mises_strain)
#from mymetal.calculate.calmechanics.stretch import (cal_relative_stretch, cal_stretch)
#from mymetal.calculate.calmath.matrix import (hermite_normal_form)
from mymetal.calculate.calenergy.surfenergy import cal_surface_energy
from mymetal.calculate.electronic_structure.plotter import *
from mymetal.calculate.electronic_structure.universal import *
from mymetal.calculate.material_science.schmid import *
from mymetal.calculate.calqm.kpoints import get_kpoints_by_size, get_size_by_distance, get_distance_by_size, cal_reciprocal_matrix, cal_reciprocal_matrix2, get_lattice_information
__all__ = [
            'compare_atoms', 'cal_atom_num', 'relative_diff', 'cal_mismatch', 'filter_results', 'cal_stretch_lattice',
            'cal_deform_matrix',

            'cal_principal_and_shear_strain', 'cal_principal_and_shear_strain_root', 
            'cal_strain_matrix', 'cal_strain_matrix_root', 'cal_von_mises_strain',

            'cal_stretch', 'cal_relative_stretch',

            'hermite_normal_form',

            'cal_surface_energy',

            'cal_fcc_schmid_factors',

            'get_kpoints_by_size', 'get_size_by_distance', 'get_distance_by_size','cal_reciprocal_matrix', 'cal_reciprocal_matrix2', 'get_lattice_information',

]
