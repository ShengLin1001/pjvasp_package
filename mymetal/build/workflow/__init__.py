"""
mymetal.build.workflow

This subpackage provides directory/input generators that drive the VASP workflow
scripts in ``vasp_utils``. They start from a relaxed reference directory
(``y_full_relax``) and lay out the deformed or stretched calculations underneath it.

Modules:
    - general: Shared helpers for copying VASP inputs and comparing lattices.
    - hoec: Generates the deformed input directories for higher-order (2nd/3rd/4th)
      elastic constants, plus the mode manifest read by ``mymetal.post.hoec_energy``.
"""
