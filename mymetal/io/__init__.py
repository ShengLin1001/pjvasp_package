"""
mymetal.io

This subpackage provides input and output utilities for interacting with VASP files 
and post-processing functions used in material simulations.

It includes functions for reading and writing VASP input/output files, as well 
as creating and writing content to files for post-processing tasks.

Modules:
    - vasp: Reads and writes VASP POSCAR/CONTCAR files (adapted from ase.io.vasp).
    - extxyz: Converts Extended XYZ trajectory files into a list of ASE Atoms.
    - general: Reads delimited files into a DataFrame and writes DataFrames to
      formatted plain text.
    - post: Constructs and writes post-processing content to files.
"""
