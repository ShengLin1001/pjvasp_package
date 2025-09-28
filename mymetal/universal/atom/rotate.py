import numpy as np

# taken from my matlab code written in 2023
def cal_rotate(ini_basis, theta, a10, a20, a30):
    """
    Rotate a vector around an arbitrary axis by a given angle using Rodrigues' rotation formula.

    Args:
        ini_basis (array-like): Initial basis vector (3-element array or list), column vector.
        theta (float): Rotation angle in radians.
        a10, a20, a30 (float): Components of the rotation axis vector.

    Returns:
        np.ndarray: Rotated vector (3-element array).
    """
    # Normalize the rotation axis
    axis = np.array([a10, a20, a30], dtype=float)
    axis = axis / np.linalg.norm(axis)
    a1, a2, a3 = axis

    cs = np.cos(theta)
    si = np.sin(theta)

    # Rodrigues rotation matrix
    trans_matrix = np.array([
        [cs + (1 - cs) * a1**2,      (1 - cs) * a1 * a2 - si * a3, (1 - cs) * a1 * a3 + si * a2],
        [(1 - cs) * a1 * a2 + si * a3, cs + (1 - cs) * a2**2,      (1 - cs) * a2 * a3 - si * a1],
        [(1 - cs) * a1 * a3 - si * a2, (1 - cs) * a2 * a3 + si * a1, cs + (1 - cs) * a3**2]
    ])

    ini_basis = np.array(ini_basis, dtype=float).reshape(3)  # Ensure shape (3,)
    return trans_matrix @ ini_basis  # Matrix multiplication
