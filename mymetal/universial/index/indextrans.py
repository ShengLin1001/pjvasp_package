"""
indextrans submodule

This submodule contains functions for converting between three-index and four-index notation for HCP crystals.

Functions:
    - three_index_to_four_index: Converts a three-index notation to a four-index notation for HCP crystals, or vice versa.
"""


import numpy as np

def three_index_to_four_index(index: list = None, reverse: bool = False) -> list:
    """
    Converts a three-index notation to a four-index notation for HCP crystals, or vice versa.

    Args:
        index (list): A list of three elements [U, V, W] for forward conversion or 
                      four elements [u, v, t, w] for reverse conversion.
        reverse (bool): If True, converts from four-index to three-index notation. 
                        If False, converts from three-index to four-index notation.

    Returns:
        list: A list representing the converted index.
        
    Raises:
        ValueError: If the index is None or does not have the correct number of elements 
                    for the specified conversion direction.
    """
    if index is None:
        raise ValueError("Index must be a list of three/four elements [u, v, w] / [u, v, t, w]")
    
    if reverse:
        if len(index) != 4:
            raise ValueError("For reverse conversion, index must be a list of four elements [u, v, t, w]")
        
        index = np.array([index[0], index[1], index[3]])
        # 四指数转换为三指数矩阵
        transformation_matrix = np.array([
            [2, 1, 0],
            [1, 2, 0],
            [0, 0, 1]
        ])
        four_index = np.array(index).reshape(3, 1)
        three_index = np.dot(transformation_matrix, four_index).flatten()
        return three_index.tolist()
    else:
        if len(index) != 3:
            raise ValueError("For forward conversion, index must be a list of three elements [U, V, W]")
        # 三指数转换为四指数矩阵
        transformation_matrix = np.array([
            [2/3, -1/3, 0],
            [-1/3, 2/3, 0],
            [0, 0, 1]
        ])
        three_index = np.array(index).reshape(3, 1)
        #print(three_index)
        four_index = np.dot(transformation_matrix, three_index).flatten()
        #print(four_index)
        index = four_index
        four_index = np.array([index[0], index[1], -(index[0] + index[1]), index[2]])
        #print(four_index)
        return four_index.tolist()