"""
mymetal.build

This subpackage provides functions for building structures, including bulk
cells, thin films (surfaces, heterostructures, GSFE models), and VASP workflow
input generators. It is organized into subpackages that focus on different
aspects of structure building.

Subpackages:
    - bulk: (hkl)-oriented crystal cells and generalized stacking fault energy
      (GSFE) models for FCC and HCP structures.
    - film: Thin-film/slab construction, primitive-cell and cubic-cell finding,
      heterostructure building, hydroxyl passivation, and film stretching.
    - workflow: Directory/input generators that drive the VASP workflow scripts
      in ``vasp_utils`` (higher-order elastic constants, KPAR/NCORE timing).
"""


