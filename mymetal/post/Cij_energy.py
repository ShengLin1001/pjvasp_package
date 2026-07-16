"""
Cij_energy post-processing submodule.

This module provides functions to post-process LAMMPS deformation calculations to extract elastic constants (Cij) from energy data.


Functions:
    - post_lammps_Cij_energy: Main function to process deformation data and generate plots
    - check_ldata: Check consistency of reference states
    - fit_cij_energy: Fit energy-strain data and prepare the plotting quantities
    - write_cij_energy: Write calculated Cij values to a text file
    - read_cij_energy: Read calculated Cij values back from the text file
    - read_deform_data: Read deformation data from a specified directory

Change log:
    - Written by B. Y. (https://github.com/BinglunYin/vasp_utils/blob/master/vasp_workflow_bulk/yin_vasp_plot_cij_energy.py).
    - Revised by J. P. on 2025-11-04 to adapt to LAMMPS output format.
    - Added read_cij_energy by J. P. to mirror write_cij_energy.
"""
# This page is taken from https://github.com/BinglunYin/vasp_utils/blob/master/vasp_workflow_bulk/yin_vasp_plot_cij_energy.py
# Changed to fit the lammps output by J. P. 
# 2025.11.04

import numpy as np
from myvasp import vasp_func as vf
import os
from pathlib import Path
from ase import Atoms
from ase.io import read

from mymetal.universal.plot.workflow import my_plot_cij_energy

def post_lammps_Cij_energy(dir: str = None,
                           refcontcar: str = None,
                           save_fig_path: str = './y_post_cij_energy.pdf',
                           save_txt_path: str = './y_post_cij_energy.txt',):
    """
    Post-process LAMMPS deformation calculations to extract elastic constants (Cij) from energy.

    Args:
        dir (str): Root directory containing subdirectories for different deformation modes.
        refcontcar (str): Path to reference structure in LAMMPS data format.
        save_fig_path (str): File path to save the resulting plot (PDF).
        save_txt_path (str): File path to save Cij and derived properties (TXT).

    Returns:
        None
    """
    # zero-strain reference
    # style='atomic' 必须显式指定：LAMMPS write_data 用 atom_style atomic（id type x y z [ix iy iz]），
    # 而 ASE read_lammps_data 默认 style='full'，会按 7/10 字段解析 -> RuntimeError。只取几何（cell/natoms），元素无所谓。
    atoms_ref = read(refcontcar, format='lammps-data', style='atomic')
    cell_ref = np.array(atoms_ref.get_cell())
    natoms = len(atoms_ref)
    # Col vector [x, y, z]
    cvectors_ref = np.linalg.norm(cell_ref, axis=0)
    # row vector [a1, a2, a3]
    rvectors_ref = np.linalg.norm(cell_ref, axis=1)

    # subdir within dir
    dirlist=[
        'y_cij_energy_c11', 
        'y_cij_energy_c12', 
        'y_cij_energy_c13', 
        'y_cij_energy_c33', 
        'y_cij_energy_c44', 
    ]

    ldata = []   # list of data 

    for i in np.arange(len(dirlist)):
        dirn = dirlist[i]
        dirnpath = os.path.join(dir, dirn)
        print(dirn) 
        data = read_deform_data(dirnpath, atoms_ref)
        ldata.append(data) 

    check_ldata(ldata)

    dict_fit = fit_cij_energy(ldata)
    my_plot_cij_energy(
        ldata=ldata,
        dict_fit=dict_fit,
        if_save=True,
        savefile=save_fig_path,
    )
    write_cij_energy(dict_fit['cij'], save_txt_path)

    return None

def check_ldata(ldata: list = None):
    """
    Check consistency of reference states among different deformation datasets.

    Args:
        ldata (list): List of dictionaries containing deformation data.

    Raises:
        AssertionError: If energy, energy density, or stress differ at reference points.
    """
    iref0 = ldata[0]['iref']

    for i in np.arange(len(ldata)):
        iref1 = ldata[i]['iref']
        
        vf.confirm_0( ldata[0]['e'][iref0]  - ldata[i]['e'][iref1]  )
        vf.confirm_0( ldata[0]['ed'][iref0] - ldata[i]['ed'][iref1] )
        vf.confirm_0( ldata[0]['s'][iref0]  - ldata[i]['s'][iref1]  )

def fit_cij_energy(ldata: list = None) -> dict:
    """Fit the Cij energy-strain data and prepare quantities for plotting.

    Args:
        ldata (list): List of dictionaries containing deformation data.

    Returns:
        dict: Fit grid, fitted energy-density curves, equilibrium shifts,
            elastic constants, and stress slopes.
    """
    xi = np.linspace(min(ldata[0]['e']), max(ldata[0]['e']), 1000)
    lp0 = []
    le0 = []
    ly_fit = []
    for data in ldata:
        e = data['e']
        ed = data['ed']
        param = np.polyfit(e, ed, 2)
        lp0.append(param[0])
        le0.append(-param[1] / (2 * param[0]))
        ly_fit.append(np.polyval(param, xi))

    c11 = lp0[0] * 2
    c12 = lp0[1] - c11
    c33 = lp0[3] * 2
    c13 = lp0[2] - (c11 + c33) / 2
    c44 = lp0[4] * 2
    cij = np.array([c11, c12, c13, c33, c44])
    lslope = np.array([
        [c11, c12, c13, 0, 0, 0],
        [c11 + c12, c12 + c11, c13 * 2, 0, 0, 0],
        [c11 + c13, c12 + c13, c13 + c33, 0, 0, 0],
        [c13, c13, c33, 0, 0, 0],
        [0, 0, 0, c44, 0, 0],
    ])
    return {
        'xi': xi,
        'ly_fit': ly_fit,
        'le0': le0,
        'cij': cij,
        'lslope': lslope,
    }

def write_cij_energy( cij_hcp: np.array = None, save_txt_path: str = './y_post_cij_energy.txt',):
    """
    Write calculated Cij values and derived elastic properties to a text file.

    Args:
        cij_hcp (np.array): Array of elastic constants [C11, C12, C13, C33, C44].
        save_txt_path (str): Path to save the output TXT file.

    Returns:
        None
    """
    f = open(save_txt_path,'w+')
    f.write('# Cij from energy method: \n' )

    f.write('\n%12s %12s %12s %12s %12s \n'  
        %('C11', 'C12', 'C13', 'C33', 'C44') )

    f.write('%12.4f %12.4f %12.4f %12.4f %12.4f \n' \
        %(cij_hcp[0], cij_hcp[1], cij_hcp[2], cij_hcp[3], cij_hcp[4]) )



    from myalloy import calc_elastic_constant as ce 
    E_x, E_z, nu_xy, nu_xz, mu_xz = ce.calc_transverse_isotropy( cij_hcp ) 
  

    f.write('\n%12s %12s %12s %12s %12s \n'  
        %('E_x', 'E_z', 'nu_xy', 'nu_xz', 'mu_xz') )

    f.write('%12.4f %12.4f %12.4f %12.4f %12.4f \n' \
        %( E_x,   E_z,   nu_xy,   nu_xz,   mu_xz ) )



    f.write('\n%25s %25s \n'  
        %('E_x/(1-nu_xy)', 'nu_xz/(1-nu_xy)' ) )

    f.write('%25.4f %25.4f \n' \
        %( E_x/(1-nu_xy),   nu_xz/(1-nu_xy)  ) )


    f.write('\n%25s \n'  
        %('E_x/(1-nu_xy**2)'  ) )

    f.write('%25.4f \n' \
        %( E_x/(1-nu_xy**2)   ) )


    f.close()

    print(f"✅ Cij-energy written to '{save_txt_path}'")

def read_cij_energy(save_txt_path: str = './y_post_cij_energy.txt') -> dict:
    """
    Read elastic constants written by :func:`write_cij_energy`.

    Parses the ``C11 C12 C13 C33 C44`` block (always present) and, when the
    file contains them, the derived transverse-isotropy quantities
    ``E_x E_z nu_xy nu_xz mu_xz``. Each block is a header line of names
    followed by one line of numeric values.

    Args:
        save_txt_path (str): Path to the file written by ``write_cij_energy``.

    Returns:
        dict: Keys ``C11, C12, C13, C33, C44`` (GPa) are always returned; keys
        ``E_x, E_z, nu_xy, nu_xz, mu_xz`` are added when that block is present.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the C11..C44 block cannot be located.
    """
    p = Path(save_txt_path)
    if not p.is_file():
        raise FileNotFoundError(f"File not found: {p}")

    lines = p.read_text(encoding='utf-8', errors='ignore').splitlines()
    blocks = [['C11', 'C12', 'C13', 'C33', 'C44'],
              ['E_x', 'E_z', 'nu_xy', 'nu_xz', 'mu_xz']]

    out = {}
    for i, ln in enumerate(lines):
        names = ln.split()
        for block in blocks:
            if names[:5] == block and i + 1 < len(lines):
                vals = lines[i + 1].split()
                for k, v in zip(block, vals[:5]):
                    out[k] = float(v)

    if 'C11' not in out:
        raise ValueError(f"Could not find C11..C44 block in {p}.")
    return out

def read_deform_data(dirn: str = None, atoms_ref: Atoms = None) -> dict:
    """
    Read energy, stress, and strain data from a deformation calculation directory.

    Args:
        dirn (str): Directory containing deformation output (y_post_data.txt and movie.lammpstrj).
        atoms_ref (Atoms): Reference structure for calculating strain.

    Returns:
        dict: A dictionary containing:
            - 'e' (np.array): Applied strain values.
            - 'ed' (np.array): Energy density relative to reference.
            - 's' (np.array): Stress tensor components (GPa).
            - 'iref' (int): Index of reference (zero-strain) structure.
    """
    x  = np.array([])      # applied strain, with respect to the ref state 
    e_ss_V = np.zeros(6)   # strain components 
    
    #os.chdir(dirn)
    jobn, Etot, Eent, pres = vf.vasp_read_post_data(f'{dirn}/y_post_data.txt')

    iref = -1    # id of the ref state 
    latoms = read(f'{dirn}/movie.lammpstrj', format='lammps-dump-text', index=':')

    for i in np.arange(len(jobn)):
        temp = float(jobn[i])-1
        x = np.append(x, temp )
        
        if np.abs(temp)<1e-10:
            iref = i 

        temp = calc_strain(atoms_ref.get_cell()[:], latoms[i].get_cell()[:])
        e_ss_V = np.vstack([ e_ss_V, temp])
    

    if iref < 0:
        os.exit('ABORT: no ref state. ')    
    
    e_ss_V = np.delete(e_ss_V, 0, 0)
    np.savetxt(f'{dirn}/y_post_calc_strain_all.txt', e_ss_V )


    # energy density
    Etot = Etot - Etot[iref]
    ed   = Etot / latoms[iref].get_volume()  * vf.phy_const('qe')*1e21    
        
    # stress     
    s = pres *(-0.1)   # GPa

    temp = s.copy()
    s[:,3] = temp[:,4] # yz
    s[:,4] = temp[:,5] # xz
    s[:,5] = temp[:,3] # xy
    
      
    data = {}
    data['e']  = x
    data['ed'] = ed
    data['s']  = s
    data['iref']  = iref 
    return data 



def calc_strain(latt1, latt2):
    # latt1 - ref
    # latt2 - deformed 

    F = np.linalg.solve(latt1, latt2)
    F = F.T 
   
    # ================================
    e_ss      = 1/2*( F+F.T ) - np.eye(3)
    e_green   = 1/2*( F.T@F - np.eye(3) )
    e_almansi = 1/2*( np.eye(3) - np.linalg.inv(F@F.T) )

    e_ss_V = vf.calc_to_Voigt(e_ss)

    return e_ss_V
