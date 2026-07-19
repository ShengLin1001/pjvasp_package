"""
mymetal.slurm

This subpackage provides a tool-agnostic engine for generating Slurm job
scripts and submitting batch jobs across many subdirectories (e.g. VASP,
n2p2, LAMMPS workflows). The standalone ``pei_slurm_univ_submit.py`` command
is its thin argparse/preset CLI and must keep parameter names aligned with the
library entry point.

Modules:
    - submit: Builds the base/sequential Slurm scripts, splits job directories
        into scheduling lanes, and drives submission in the parallel /
        each-subdir / single-alloc modes via ``pei_slurm_univ_submit``.
"""
