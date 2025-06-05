"""
kpoints module

This module provides functions to generate k-points for electronic structure calculations,
calculate reciprocal lattice vectors, and display lattice information from VASP POSCAR files.
It includes functions for Monkhorst-Pack and Gamma-centered k-points, as well as methods
to compute reciprocal lattice vectors using both cross product and matrix inversion methods.

Functions:
    - get_kpoints_by_size(): Generates Monkhorst-Pack and Gamma-centered k-points.
    - get_size_by_distance(): Computes automatic k-point mesh based on RK or KSPACING.
    - get_distance_by_size(): Estimates RK value for a given KPOINTS grid size.
    - cal_reciprocal_matrix(): Calculates reciprocal lattice vectors using cross product.
    - cal_reciprocal_matrix2(): Calculates reciprocal lattice vectors using matrix inversion.
    - get_lattice_information(): Displays real-space and reciprocal lattice information from a VASP POSCAR file.
"""

from ase.io.vasp import read_vasp
import numpy as np
from ase.dft.kpoints import *


def get_kpoints_by_size(size: tuple=(1, 1, 1), offset: tuple=(0.5, 0.5, 0.5)):
    """
    Generate Monkhorst-Pack and Gamma-centered k-points.

    For each direction with an even number of k-points, an offset is added
    to shift the k-point grid appropriately.

    Args:
        size (tuple): Grid size in each reciprocal direction (nx, ny, nz).
        offset (tuple): Offset values for each direction.

    Returns:
        tuple: A tuple containing:
            - mpkpoints (ndarray): Original Monkhorst-Pack k-points in fractional coordinates.
            - gkpoints (ndarray): Offset-adjusted Monkhorst-Pack k-points.

    Example:
        >>> from mymetal.calculate.calqm.kpoints import get_kpoints_by_size
        >>> import matplotlib.pyplot as plt
        >>> mpk, gk = get_kpoints_by_size((4, 4, 4), (0.5, 0.5, 0.5))
        >>> plt.plot(mpk[:, 0], mpk[:, 1], 'o')
        >>> plt.plot(gk[:, 0], gk[:, 1], 'x')
        >>> plt.show()
    """   
    mpkpoints = np.array(monkhorst_pack(size))
    gkpoints = mpkpoints.copy() 
    # Apply the transformed offset only to columns corresponding to even size elements
    for i in range(len(size)):
        if size[i] % 2 == 0: 
            transformed_offset = offset[i] / size[i] 
            gkpoints[:, i] += transformed_offset
    
    return mpkpoints, gkpoints

def get_size_by_distance(file: str = None, rk: int = 100, kspacing: float =None) -> tuple:
    """
    Compute automatic k-point mesh based on RK or KSPACING scheme in VASP.

    This function reads the unit cell from a VASP CONTCAR file, calculates
    the reciprocal lattice lengths (in 2π/Å), and determines the number of 
    k-points along each direction using:

    - **Old VASP method** (RK): uses rounding (`round`) via `+0.5` and `int(...)`
    - **New VASP method** (KSPACING): uses the ceiling function (`ceil`)
        

    Args:
        file (str): Path to the VASP CONTCAR file.
        rk (int): RK product for the old method (unit: 2π/Å). R_k = 2π / KSPACING.
        kspacing (float, optional): KSPACING value (unit: Å⁻¹). If None, derived from rk.

    Returns:
        tuple: (old_auto_kpoints, new_kspacing_kpoints), each a (3,) np.ndarray.

    Example:
        >>> oldk, newk = get_kpoints_by_rk_or_kspacing(file='CONTCAR', rk=81)
        >>> print("Old method:", oldk)
        Old method: [32 32  3]
        >>> print("New method:", newk)
        New method: [33 33  3]
    """
    if kspacing is None:
        kspacing = 2 * np.pi / rk  # units: 1/Å, rk - 2*pi/Å
    # kspacing (R_{k}) in R_k = 2π / KSPACING
    np.set_printoptions(precision=8, suppress=True)
    atoms = read_vasp(file)
    reciprocal_cell_length = np.array(atoms.cell.reciprocal().cellpar())[:3] # units: 2*pi/Å
    
    # see here: https://www.vasp.at/wiki/index.php/KSPACING
    # +0.5 is to round up, 四舍五入
    old_auto_kpoints =  np.maximum(1, (rk * reciprocal_cell_length + 0.5).astype(int))
    # KSPACING tag 用的是上限函数，3.2 => 4
    new_kspacing_kpoints = np.maximum(1, np.ceil(reciprocal_cell_length * 2 * np.pi / kspacing).astype(int))
    return old_auto_kpoints, new_kspacing_kpoints

def get_distance_by_size(file: str='CONTCAR', size: tuple=None)-> list:
    """
    Estimate the matching RK value for a given KPOINTS grid size.

    This function calculates the reciprocal lattice vector lengths,
    computes the unit length per k-point, and searches for the highest
    RK value (radius cutoff in reciprocal space) that matches the target grid size.

    Args:
        file (str): Path to the VASP CONTCAR file. Default is 'CONTCAR'.
        size (tuple): A tuple of 3 integers representing the desired k-point grid size.

    Returns:
        list: The maximum matching RK value if found (int), otherwise the estimated RK values (numpy.ndarray).
    """
    #np.set_printoptions(precision=8, suppress=True)
    atoms = read_vasp(file)
    reciprocal_cell_length = np.array(atoms.cell.reciprocal().cellpar())[:3] # units: 2*pi/Å
    nsize = np.array(size)
    unit_lengths = reciprocal_cell_length / nsize              # units: 2*pi/Å

    rks = 1/ unit_lengths
    #print(r"$R_{k}$:", rks)

    matchrk = []
    for rk in range(1, 101):
        otemp, ntemp = get_size_by_distance(file=file, rk=rk)
        if (otemp == nsize).all() or (ntemp == nsize).all():
            matchrk.append(rk)
    if matchrk != []:
        return max(matchrk)
    else:
        return rks

def cal_reciprocal_matrix(cell_matrix: np.ndarray=None, scale: float=1.0):
    """
    Calculate the reciprocal lattice vectors using cross product.

    The result is in units of (2π/Å) if scale=2π.

    Args:
        cell_matrix (ndarray): Real-space lattice vectors (3x3 matrix).
        scale (float): Scaling factor, typically 2π for reciprocal lattice.

    Returns:
        tuple: Reciprocal lattice vectors (b1, b2, b3) as three 1D arrays.

    Example:
        >>> b1, b2, b3 = cal_reciprocal_matrix(cell_matrix, scale=2*np.pi)
    """
    # unit: 2*np.pi/N, for comparing with VASP
    V = np.dot(cell_matrix[0], np.cross(cell_matrix[1], cell_matrix[2]))
    b1 = scale * np.cross(cell_matrix[1], cell_matrix[2]) / V
    b2 = scale * np.cross(cell_matrix[2], cell_matrix[0]) / V
    b3 = scale * np.cross(cell_matrix[0], cell_matrix[1]) / V
    return b1, b2, b3

def cal_reciprocal_matrix2(cell_matrix: np.ndarray=None, scale: float=1.0):
    """
    Calculate reciprocal lattice vectors using matrix inversion.

    This method is mathematically equivalent to the cross-product method
    but uses matrix inversion and transpose.

    Args:
        cell_matrix (ndarray): Real-space lattice vectors (3x3 matrix).
        scale (float): Scaling factor, typically 2π for reciprocal lattice.

    Returns:
        ndarray: A 3x3 reciprocal lattice matrix.

    Example:
        >>> reciprocal_matrix = cal_reciprocal_matrix2(cell_matrix, scale=2*np.pi)
    """
    inv_cell_matrix = np.linalg.inv(cell_matrix)
    reciprocal_matrix =  scale * inv_cell_matrix.T
    # unit: 2*np.pi
    # [[b1];
    #  [b2];
    #  [b3]]
    return reciprocal_matrix

def get_lattice_information(file: str = None):
    """
    Display real-space and reciprocal lattice information from a VASP POSCAR file.

    Reads the atomic structure, prints the lattice vectors, reciprocal vectors,
    Bravais lattice type, and high-symmetry points.

    Args:
        file (str): Path to the VASP CONTCAR file.

    Example:
        >>> get_lattice_information('CONTCAR')
    """
    # Read the VASP CONTCAR file to get the atomic structure
    np.set_printoptions(precision=8, suppress=True)
    atoms = read_vasp(file)
    cell = np.array(atoms.cell)
    reciprocal_cell = np.array(atoms.cell.reciprocal())

    print('Cell vectors:')
    print(cell)
    print('Reciprocal cell vectors:')
    print(reciprocal_cell)

    lat = atoms.cell.get_bravais_lattice()
    print("Lattice type:",lat.name)
    print("High symmetry points:", list(lat.get_special_points()))
    lat.plot_bz()