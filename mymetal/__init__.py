"""
mymetal package

The ``mymetal`` package provides a comprehensive set of tools and utilities for
materials science and computational modeling. It includes subpackages that cover
structure building, property calculation, I/O, post-processing, machine-learning
potential tooling, Slurm submission, contact-resonance analysis, and shared
universal helpers. It is designed to support the development, analysis, and
optimization of materials, particularly for thin film and other advanced material
systems.

Subpackages:
    - build: Functions for constructing bulk and film structures (surfaces,
      heterostructures, GSFE models, hydroxyl passivation, stretching) and
      workflow input generators for VASP.
    - calculate: Mathematical and computational utilities for material property
      calculations (energy, mechanics, mismatch, QM k-points, electronic
      structure, material-science quantities such as Schmid factors).
    - cr: Contact-resonance (CR) post-processing scripts, including stiffness
      fitting from force-distance curves and k-contact plotting.
    - io: Input/output for VASP files, Extended XYZ, and general delimited
      files, plus post-processing file construction and writing.
    - ml: Machine-learning tools, in particular the n2p2 neural-network
      potential wrapper (dataset, training workflow, symmetry-function
      parameter generation and CUR selection).
    - post: Post-processing functions for analyzing and visualizing VASP/LAMMPS
      simulation results (elastic constants, GSFE, NEB, stretch, convergence,
      etc.).
    - slurm: Tool-agnostic engine for generating Slurm job scripts and
      submitting batch jobs across many subdirectories.
    - universal: General-purpose utilities (atom manipulation, input checks,
      data transforms, math/matrix helpers, plotting, printing, file search)
      shared across the package.

Usage:
    To use this package, import the relevant subpackage and access the necessary
    functions and classes:

    Example:
    >>> from mymetal.build.film.stretch import generate_film
    >>> generate_film('Au', 'fcc', num_layers = 12, a_fcc = 4.08, slice_plane = (1, 1, 1))

    Each subpackage includes detailed documentation and example usage to guide
    users through its functionalities. You can see the examples and documentation
    for each subpackage by importing it and using the ``help()`` function.

Installation:
    To install ``mymetal`` in editable mode from the repository root:

    >>> python -m pip install -e .

Authors:
    - J. Pei

License:
    This package is licensed under MIT. See LICENSE for more information.

"""
# Importing submodules
# from . import build
# from . import calculate
# #from . import example
# from . import io
# #from . import ml
# from . import post
# from . import universal