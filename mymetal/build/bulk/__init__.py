"""
mymetal.build.bulk

This subpackage provides functions for creating (hkl)-oriented crystal cells
(replacing the older ``generate_bulk_from_film()`` approach) and generalized
stacking fault energy (GSFE) models.

Modules:
    - create: Creates FCC (111), HCP (0001)/(10-10)/(10-11) and other oriented
      crystal cells for dislocation and slip-system analysis.
    - gsfe: Builds GSFE supercells for FCC and HCP from the cells produced by
      ``create``.
"""
