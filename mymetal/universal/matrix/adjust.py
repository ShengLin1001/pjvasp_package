"""
adjust submodule

This submodule contains functions for adjusting the size of a matrix.

Functions:
    - adjust_matrix: Adjusts the specified row and/or column of a matrix by setting all their elements to a specified value.
"""

from numpy import array
def adjust_matrix(matrix: array= None,row_num: int = None, column_num: int = None, value: float = 0):
    """Adjusts the specified row and/or column of a matrix by setting all their elements to a specified value.

    Args:
        matrix (np.ndarray): The matrix to adjust.
        row_num (int, optional): The row index to adjust. If None, no row will be adjusted.
        column_num (int, optional): The column index to adjust. If None, no column will be adjusted.
        value (float): The value to set the specified row and/or column elements to.

    Raises:
        ValueError: Input matrix cannot be None.
        IndexError: row_num is out of bounds.
        IndexError: column_num is out of bounds.

    Returns:
        array: a changed matrix
    """
    if matrix is None:
        raise ValueError("Input matrix cannot be None.")
    if row_num != None:
        if row_num >= matrix.shape[0] or row_num < 0:
            raise IndexError("row_num is out of bounds.")
        matrix[row_num, :] = value  # Set the entire row to 'value'

    if column_num != None:
        if column_num >= matrix.shape[1] or column_num < 0:
            raise IndexError("column_num is out of bounds.")
        matrix[:, column_num] = value  # Set the entire column to 'value'
    return matrix