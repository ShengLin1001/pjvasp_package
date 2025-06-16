"""
mymetal package

The `mymetal` package provides a comprehensive set of tools and utilities for materials science and computational modeling. 
It includes a range of subpackages that cover various aspects of material property calculations, structure building, post-processing, machine learning, and more.
This package is designed to support the development, analysis, and optimization of materials, particularly for thin film and other advanced material systems.

Subpackages:
    - build: Contains functions and classes for constructing material structures, including thin films and related materials.
    - calculate: Provides mathematical and computational utilities for material property calculations, simulations, and other scientific computations.
    - example: Includes example scripts and demonstrations to help users understand how to use the package and its features.
    - io: Handles input and output operations, including file reading and writing, data import/export, and interfacing with external tools.
    - ml: Implements machine learning algorithms and tools for predictive modeling, analysis, and optimization of materials.
    - post: Contains post-processing functions for analyzing and visualizing results from simulations, experiments, and other calculations.
    - universal: Provides general-purpose utilities and functions that are commonly used across the package, ensuring modularity and reusability.

Usage:
    To use this package, import the relevant subpackage and access the necessary functions and classes:
    
    Example:
    >>> from mymetal.build.film.stretch import generate_film
    >>> generate_film('Au', 'fcc', num_layers = 12, a_fcc = 4.08, slice_plane = (1, 1, 1))
    
    Each subpackage includes detailed documentation and example usage to guide users through its functionalities.
    You can see the examples and documentation for each subpackage by importing it and using the `help()` function:

Installation:
    To install `mymetal`, use the following:
    >>> pip install mymetal

Authors:
    - J. Pei

License:
    This package is licensed under [MIT]. See LICENSE for more information.

"""
# Importing submodules
# from . import build
# from . import calculate
# #from . import example
# from . import io
# #from . import ml
# from . import post
# from . import universal