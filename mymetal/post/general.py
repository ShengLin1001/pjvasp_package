"""
general submodule

This module provides general utility functions for post-processing VASP calculations,
including sorting data, polynomial fitting, reading CONTCAR files from directories,
and extracting structural information.

Functions:
    - my_sort: Sort convergence data based on input parameters.
    - my_ployfit: Perform polynomial fitting of given x-y data.
    - my_read_y_dir_contcar: Read CONTCAR files from y_dir subfolders.
    - get_structure_info: Extract lattice vectors, cell parameters, and c/a ratios.
"""

# Written by J. P.
# 2025.11.04

import numpy as np
import os
from myvasp import vasp_func as vf 
from ase.io.vasp import read_vasp

def my_sort(refx: list=[], Etot: list=None, Eent: list = None, pres: list = None, reverse: bool = False) -> tuple:
    """
    Sort convergence data based on input parameter (cutoff or k-point grid).

    Args:
        refx (list): List or array of input values used for sorting (1D or 2D).
        Etot (list): Total energy values.
        Eent (list): Electronic entropy values.
        pres (list): Pressure values.
        reverse (bool): If True, sort in descending order.

    Returns:
        tuple:
            x_sorted (np.ndarray): Sorted input values.
            Etot_sorted (np.ndarray): Sorted total energies.
            Eent_sorted (np.ndarray): Sorted entropy values.
            pres_sorted (np.ndarray): Sorted pressure values.
    """
    refx = np.array(refx)
    Etot = np.array(Etot)
    Eent = np.array(Eent)
    pres = np.array(pres)
    x = refx
    if x.ndim == 1:
        sort_idx = np.argsort(x)        # 1D 
    else:
        sort_idx = np.argsort(x[:, 0])  # 2D 

    if reverse:
        sort_idx = sort_idx[::-1]

    x_sorted = x[sort_idx]
    Etot_sorted = np.array(Etot)[sort_idx] 
    Eent_sorted = np.array(Eent)[sort_idx] 
    pres_sorted = np.array(pres)[sort_idx]
    return x_sorted, Etot_sorted, Eent_sorted, pres_sorted


def my_ployfit(x: np.array = None, y: np.array = None, deg: int = 2) -> tuple:
    """
    Perform polynomial fitting of given x–y data.

    Args:
        x (np.ndarray): Independent variable values.
        y (np.ndarray): Dependent variable values.
        deg (int): Degree of polynomial. Default 2.

    Returns:
        tuple: 
            - coeffs (np.ndarray): Polynomial coefficients, highest degree first.
            - y_fit (np.ndarray): Evaluated polynomial at input x values.
    """
    if x is None or y is None:
        raise ValueError("x and y must be provided for polynomial fitting.")
    if len(x) != len(y):
        raise ValueError("x and y must have the same length.")

    coeffs = np.polyfit(x, y, deg)
    p = np.poly1d(coeffs)
    y_fit = p(x)

    return coeffs, y_fit


def my_read_y_dir_contcar(dir: str = '.', post_data_file: str = './y_post_data.txt', file: str = 'CONTCAR') -> list:
    """
    Read CONTCAR files from y_dir subfolders.

    Args:
        dir (str): Directory containing subfolders with CONTCAR files. Default current directory.
        post_data_file (str): Path to post data file listing job names. Default './y_post_data.txt'.
        file (str): Name of the atoms file in each subfolder. Default 'CONTCAR'.

    Returns:
        list: List of ASE Atoms objects corresponding to each job.
    """
    jobn, Etot, Eent, pres = vf.vasp_read_post_data(post_data_file)
    latoms = []
    for job in jobn:
        contcar_path = os.path.join(dir, job, file)
        if os.path.isfile(contcar_path):
            atoms = read_vasp(contcar_path)
            latoms.append(atoms)
        else:
            raise FileNotFoundError(f"❌ File {contcar_path} does not exist.")
    
    return latoms

def get_structure_info(latoms: list = None) -> tuple:
    """
    Extract lattice vectors, row vectors, cell parameters, and c/a ratios from Atoms objects.

    Args:
        latoms (list): List of ASE Atoms objects.

    Returns:
        tuple:
            - lcvectors (np.ndarray): Column vector lengths for each structure.
            - lrvectors (np.ndarray): Row vector lengths for each structure.
            - lcellpars (np.ndarray): Cell parameters for each structure.
            - lca (np.ndarray): c/a ratios for each structure.
    """
    lcvectors = []
    lrvectors = []
    lcellpars = []
    lca       = [] # c/a ratio
    for atoms in latoms:
        cell = np.array(atoms.get_cell())
        cellpars = np.array(atoms.cell.cellpar())
        # Col vector [x, y, z]
        cvectors = np.linalg.norm(cell, axis=0)
        # Row vector [a1, a2, a3]
        rvectors = np.linalg.norm(cell, axis=1)
        ca = rvectors[2] / rvectors[0]

        # check the row vector [:] with cellpars[:3]
        for i in range(3):
            if not np.isclose(rvectors[i], cellpars[i], atol=1e-5):
                raise ValueError(f"Error: Row vector {rvectors[i]} and cell parameter {cellpars[i]} do not match.")
            
        lcvectors.append(cvectors)
        lrvectors.append(rvectors)
        lcellpars.append(cellpars)
        lca.append(ca)  
    return np.array(lcvectors), np.array(lrvectors), np.array(lcellpars), np.array(lca)

