"""
deformation module

This module provides functions for calculating deformation matrices. 
Deformation analysis is important for studying the mechanical behavior of materials under stress.

Functions:
    - cal_deform_matrix: Calculate the deformation matrix for a given structure.
"""

from ase import Atoms
import numpy as np

# simple version
def cal_deform_matrix(initial_cell: np.array = None, deformed_cell: np.array = None) -> np.array:
    """To calculate the deformation matrix.

    Args:
        initial_cell (np.array): initial configuration, column vector
        deformed_cell (np.array): deformed configuration, column vector

    Returns:
        np.array: a n*n matrix

    Raises:
        ValueError: Input matrix cannot be None.
    """
    if initial_cell is None or deformed_cell is None:
        raise ValueError("Input atoms cannot be None.")
    # F = B * A^(-1)
    deformation_matrix = np.dot(deformed_cell, np.linalg.inv(initial_cell))
    return deformation_matrix