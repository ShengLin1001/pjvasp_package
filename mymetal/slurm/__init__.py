"""
mymetal.slurm

This subpackage provides a tool-agnostic engine for generating Slurm job
scripts and submitting batch jobs across many subdirectories (e.g. VASP,
n2p2, LAMMPS workflows). It is the Python twin of the standalone
``pei_slurm_univ_submit`` shell tool and must stay output-compatible with it.

Modules:
    - submit: Builds the base/sequential Slurm scripts, splits job directories
        into chunks, and drives submission in the parallel / each-subdir /
        single-alloc modes via ``pei_slurm_univ_submit``.
"""
