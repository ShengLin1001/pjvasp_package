"""
matrix module

This module provides functions for performing matrix operations, particularly 
focused on the computation of Hermite normal forms and other linear algebra 
tools. These are useful for analyzing crystal structures and other material properties 
that can be represented through matrices.

Functions:
    - hermite_normal_form: Compute the Hermite normal form of a matrix.
"""

import numpy as np

def hermite_normal_form(matrix):
    """
    Convert a matrix to its Hermite Normal Form (HNF).

    Args:
        matrix (numpy.ndarray): Input matrix to be transformed.

    Returns:
        numpy.ndarray: Matrix in Hermite Normal Form.

    Notes:
        This function uses integer arithmetic to ensure the resulting
        matrix is in Hermite Normal Form, where each step maintains
        the integer nature of the matrix elements.

    Examples:
        >>> matrix = np.array([[1, 2], [0, 1]])
        >>> hermite_normal_form(matrix)
        array([[1, 0],
               [0, 1]])
    """
    m, n = matrix.shape
    hnf = matrix.copy().astype(np.int64)

    for i in range(min(m, n)):
        # Find the pivot element
        pivot_row = np.argmax(np.abs(hnf[i:, i])) + i
        if hnf[pivot_row, i] == 0:
            continue

        # Swap rows
        if pivot_row != i:
            hnf[[i, pivot_row]] = hnf[[pivot_row, i]]

        # Normalize the pivot row
        pivot = hnf[i, i]

        # Ensure the pivot is positive
        if pivot < 0:
            hnf[i] = -hnf[i]
            pivot = -pivot

        # Eliminate the current column in other rows
        for j in range(m):
            if j != i:
                q = hnf[j, i] // pivot
                hnf[j] -= q * hnf[i]

        # Process elements above the pivot
        for k in range(i):
            q = hnf[i, k] // pivot
            hnf[i] -= q * hnf[:, k]

    return hnf
