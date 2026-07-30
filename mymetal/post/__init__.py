"""
mymetal.post

This subpackage provides post-processing utilities for materials simulations and analysis.

Modules:
    - oldmain: Legacy post-processing entry point for materials simulations.
    - newmain: Newer post-processing entry point for materials simulations.
    - general: Shared helpers (sorting, polynomial fitting, CONTCAR reading,
      structural information extraction) used by the other post modules.
    - convergence: Post-processes ENCUT/k-point convergence tests into plots
      and summary text files.
    - relax_convergence: Plots ionic-relaxation energy/force convergence from
      ``pei_vasp_univ_extract_convergence`` data files.
    - Cij_energy: Second-order elastic constants from the energy-strain method.
    - hoec_energy: Second-, third- and fourth-order elastic constants from the
      energy-strain method (the higher-order sibling of Cij_energy).
    - gsfe: Post-processes generalized stacking fault energy data from VASP or
      LAMMPS into plots and summary files.
    - neb: Post-processes NEB calculations into spline plots and summary files.
    - stretch: Post-processes VASP/LAMMPS stretch calculations into plots and
      summary files.
    - E_in_1_2_bulk: Analyzes and plots E_in_1/2 bulk contour/profile results.
    - kpar_ncore: Two-panel KPAR/NCORE timing and energy benchmark
      post-processing.

"""
