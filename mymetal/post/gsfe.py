
# This page is taken from https://github.com/BinglunYin/vasp_utils/blob/master/vasp_workflow_planar_defects/yin_vasp_plot_gsfe.py
# Changed to fit the lammps output by J. P. 
# 2025.11.04

import numpy as np
from myvasp import vasp_func as vf 
import sys
from ase.io import read
import os
import matplotlib
import matplotlib.pyplot as plt
from mymetal.universal.plot.general import general_set_all_rcParams

def post_gsfe( post_data_path: str = None,
            latoms_path: str = None,
            latoms_format: str = 'lammps-dump-text',
            save_fig_path: str = './y_post_gsfe.pdf',
            save_u3_fig_path: str = './y_post_gsfe.u3.pdf',
            save_txt_path: str = './y_post_gsfe.txt',) -> None:
    """
    Post-process generalized stacking fault energy (GSFE) data from VASP or LAMMPS.

    This function reads post-calculation data and atomic configurations,
    checks for structural consistency, computes GSFE (γ-surface), and
    generates both plots and text reports.

    Args:
        post_data_path (str): Path to the VASP post data file (e.g., OUTCAR summary or parsed file).
        latoms_path (str): Path to atomic structure sequence (e.g., LAMMPS dump or trajectory).
        latoms_format (str, optional): ASE-supported format of atomic data. Defaults to 'lammps-dump-text'.
        save_fig_path (str, optional): Output path for GSFE plot. Defaults to './y_post_gsfe.pdf'.
        save_u3_fig_path (str, optional): Output path for displacement (u3) plot. Defaults to './y_post_gsfe.u3.pdf'.
        save_txt_path (str, optional): Output path for numerical results. Defaults to './y_post_gsfe.txt'.

    Returns:
        None
    """
    qe = vf.phy_const('qe')
    
    jobn, Etot, Eent, pres = vf.vasp_read_post_data(post_data_path)
    njobs = len(jobn)  # number of jobs

    if njobs < 1.5:
        sys.exit('==> ABORT! more structures needed. ')

    if abs(float(jobn[0])) > 1e-5 :
        sys.exit('==> ABORT! no reference state. ')
        
    latoms = read(latoms_path, format = latoms_format, index = ':')
    
    Asf = np.linalg.norm( \
        np.cross(latoms[0].cell[0, :], latoms[0].cell[1, :] ) )
    a11 = latoms[0].cell[0, 0]
    a22 = latoms[0].cell[1, 1]

    natoms = latoms[0].get_positions().shape[0]    
    E0bulk = Etot[0] / natoms


    dE, da3, dpos3 = check_constraints(Etot, latoms)
   
    gamma = dE/Asf *qe*1e23   
    sf, usf = find_sf_usf(gamma)

    from ase.formula import Formula
    str2 = latoms[0].get_chemical_formula()
    str2 = Formula(str2).format('latex')
    str_all = '%s\n$A =$%.4f $\\mathrm{\\AA}^2$' %(str2, Asf)  

    
    # #=========================
    plot_GSFE(jobn, gamma, da3, dpos3, latoms, str_all, save_fig_path, save_u3_fig_path)
    write_output(Asf, a11, a22, E0bulk, sf, usf, jobn, dE, gamma, da3, save_txt_path)

def check_constraints(Etot: list = None, latoms: list = None) -> tuple:
    """
    Check geometric and atomic constraints for GSFE calculations.

    Ensures that all structures in the sliding sequence:
      - Have identical atomic composition.
      - Maintain in-plane lattice constants (no relaxation in-plane).
      - Show displacement only along the slip direction.

    Args:
        Etot (list): Total energies of all configurations (in eV).
        latoms (list): Sequence of ASE `Atoms` objects representing configurations.

    Returns:
        tuple:
            - dE (np.ndarray): Energy differences relative to the reference configuration.
            - da3 (np.ndarray): Lattice vector changes along the third direction.
            - dpos3 (np.ndarray): Out-of-plane atomic displacements.
    """
    njobs = len(latoms)
    natoms = latoms[0].positions.shape[0]
    print('njobs, natoms:', njobs, natoms)

    dE = np.array([])
    da3 = np.zeros([1, 3])
    dpos3 = np.zeros([1, natoms])
    
    for i in np.arange(njobs):
        dE = np.append(dE, Etot[i]-Etot[0])

   
    # check elem
        if latoms[i].get_chemical_formula() \
            != latoms[0].get_chemical_formula():
            sys.exit('ABORT: wrong chemical formula. ')

        temp = latoms[i].get_atomic_numbers() \
            - latoms[0].get_atomic_numbers()
        vf.confirm_0( temp )

 
        # check latt 
        dlatt = latoms[i].cell[:] - latoms[0].cell[:]
        temp = np.linalg.norm(dlatt[0:2, :])
        if temp > 1e-10:
            print('dlatt:', dlatt)
            print('\n==> i, norm: {0}'.format([i, temp]) )
            sys.exit("==> ABORT: in-plane lattices relaxed. \n" )
    
        temp = dlatt[2,:].copy()
        da3 = np.vstack([ da3, temp[np.newaxis,:] ])
        

        # check pos
        dpos = latoms[i].positions - latoms[0].positions
        dposD = dpos @ np.linalg.inv(latoms[i].cell[:])
        dposD = dposD - np.around(dposD)  
        dpos = dposD @ latoms[i].cell[:]

        temp = np.linalg.norm(dpos[:,0:2])
        if temp > 1e-10:
            print('dpos', dpos)
            print('\n==> i, norm: {0}'.format([i, temp]) )
            #print(dpos)
            sys.exit("==> ABORT: atoms show in-plane relaxation. \n" )
    
        temp = dpos[:, 2].copy()
        temp = temp - temp.mean()
        dpos3 = np.vstack([ dpos3, temp[np.newaxis, :] ])
        
    da3 = np.delete(da3, 0, 0)
    dpos3 = np.delete(dpos3, 0, 0)
    
    if (dE.shape[0] != njobs)  \
        or (da3.shape[0] != njobs) \
        or (da3.shape[1] != 3) \
        or (dpos3.shape[0] != njobs) \
        or (dpos3.shape[1] != natoms):
        sys.exit("==> ABORT: wrong dimensions. \n" )

    for i in np.arange(njobs):
        temp = np.abs(da3[i, 1]*da3[1, 0] - da3[1, 1]*da3[i, 0])
        if temp > 1e-10:
            print('\n==> i, norm: {0}'.format([i, temp]) )
            sys.exit("==> ABORT: slip is not along a line. \n" )
    
    return dE, da3, dpos3

def find_sf_usf(gamma: np.ndarray = None) -> tuple:
    """
    Identify local minima and maxima in GSFE data.

    Local minima correspond to stable stacking faults (SF),
    and local maxima correspond to unstable stacking faults (USF).

    Args:
        gamma (np.ndarray): GSFE values (in mJ/m²).

    Returns:
        tuple:
            - sf (np.ndarray): Local minima (stable faults).
            - usf (np.ndarray): Local maxima (unstable faults).
    """
    njobs = gamma.shape[0]
    print(njobs)
    sf = np.array([])
    usf = np.array([])
    for i in np.arange(1, njobs-1, 1):
        if (gamma[i] < gamma[i-1]) and (gamma[i] < gamma[i+1]):
            sf = np.append(sf, gamma[i])
        if (gamma[i] > gamma[i-1]) and (gamma[i] > gamma[i+1]):
            usf = np.append(usf, gamma[i])
    return sf, usf

def write_output(Asf, a11, a22, E0bulk, sf, usf, jobn, dE, gamma, da3, save_txt_path = './y_post_gsfe.txt',):
    """
    Write GSFE results and derived quantities to a formatted text file.

    Includes surface area, bulk energy, stacking fault energies, and
    displacement-related quantities for each sliding configuration.

    Args:
        Asf (float): Surface area of the slip plane (Å²).
        a11 (float): Lattice constant along a1 direction (Å).
        a22 (float): Lattice constant along a2 direction (Å).
        E0bulk (float): Per-atom bulk energy (eV).
        sf (np.ndarray): Local minima of GSFE (stable faults).
        usf (np.ndarray): Local maxima of GSFE (unstable faults).
        jobn (list): Job or displacement indices.
        dE (np.ndarray): Energy differences (eV).
        gamma (np.ndarray): GSFE values (mJ/m²).
        da3 (np.ndarray): Lattice changes along slip direction.
        save_txt_path (str, optional): Output file path. Defaults to './y_post_gsfe.txt'.

    Returns:
        None
    """
    njobs = gamma.shape[0]
    print(njobs)
       
    f = open(save_txt_path,'w+')
    f.write('# VASP GSFE: \n' )
    f.write('# gamma = dE/Asf \n' )


    f.write('\n%16s %16s %16s %16s \n' \
        %('Asf (Ang^2)', 'a11 (Ang)', 'a22 (Ang)', 'E0_bulk (eV)' ) )
    f.write('%16.8f %16.8f %16.8f %16.8f \n\n' \
        %(Asf, a11, a22, E0bulk ))


    if sf.shape[0] > 0.5 :
        f.write('%20s: %16.8f \n' %('local min (mJ/m^2)', sf.min() ) )
    
    if usf.shape[0] > 0.5:
        f.write('%20s: %16.8f \n' %('local max (mJ/m^2)', usf.min() ) )
        
        
    f.write('\n%10s %10s %16s %10s %10s %10s %10s \n' \
        %('jobn', 'dE (eV)', 'gamma (mJ/m^2)', 
        'da31/a11', 'da32/a22', 'slip (Ang)',
        'da33 (Ang)'  
        ))
    
    for i in np.arange(njobs):
        f.write('%10s %10.4f %16.8f %10.4f %10.4f %10.4f %10.4f \n' \
            %(jobn[i], dE[i], gamma[i], 
            da3[i, 0]/a11,  da3[i, 1]/a22, np.linalg.norm(da3[i, 0:2]),
            da3[i, 2] 
            ))

    f.write(' \n')
    f.close() 
    print(f"✅ gsfe written to '{save_txt_path}'\n")

def plot_GSFE(jobn, gamma, da3, dpos3, latoms, str_all, save_fig_path = './y_post_gsfe.pdf', save_u3_fig_path = './y_post_gsfe.u3.pdf',):
    """
    Plot GSFE energy curve, displacement, and shear stress evolution.

    Generates two plots:
      1. GSFE vs. normalized slip vector (with stress and normal displacement).
      2. Atomic out-of-plane displacement (u3) profiles for all configurations.

    Args:
        jobn (list): Normalized displacement steps or job indices.
        gamma (np.ndarray): GSFE values (mJ/m²).
        da3 (np.ndarray): Lattice vector changes along slip direction.
        dpos3 (np.ndarray): Out-of-plane atomic displacements.
        latoms (list): Sequence of ASE `Atoms` objects.
        str_all (str): Summary label (chemical formula, area info).
        save_fig_path (str, optional): Output file for GSFE plot. Defaults to './y_post_gsfe.pdf'.
        save_u3_fig_path (str, optional): Output file for u3 plot. Defaults to './y_post_gsfe.u3.pdf'.

    Returns:
        None
    """
    njobs = gamma.shape[0]
    print(njobs)

    fig_subp = [3, 1]
    lg = general_set_all_rcParams(figure_subp=fig_subp, legend_framealpha=0.4,
                                  font_family=['Arial'])
    fig, ax1 = plt.subplots(fig_subp[0], fig_subp[1], sharex = True)

    xi = np.array([])
    disp = np.array([])
    for i in np.arange( njobs ):
        xi = np.append(xi, float(jobn[i]) )
        disp = np.append(disp, np.linalg.norm(da3[i, 0:2]))

    xi = xi / np.around( np.max([xi.max(), 10]) , -1)

    ax1[0].plot(xi, gamma, '-o')   
    ax1[1].plot(xi, da3[:,2], '-o')
    
    tau = np.diff(gamma) / np.diff(disp) *1e-2  #[GPa]
    x_tau = xi[0:-1].copy() + np.diff(xi)/2
    ax1[2].plot(x_tau, tau, '-o')
    ax1[2].plot([xi.min(), xi.max()], [0, 0], '--k')


    ax1[-1].set_xlabel('Normalized slip vector')
    ax1[0].set_ylabel('GSFE (mJ/m$^2$)')
    ax1[1].set_ylabel('Inelastic normal displacement ($\\mathrm{\\AA}$)')
    ax1[2].set_ylabel('Shear stress $\\tau$ (GPa)')

    dxi = np.around( xi.max()/6, 1)
    ax1[-1].set_xticks( np.arange(0, xi.max()+dxi, dxi ) )

    ax1[2].text(0, tau.min()/2, \
        '$\\tau_\\mathrm{max} =$ %.1f GPa' %( tau.max() ))

    ax1[0].text( 0.3, gamma.max()*0.2, str_all )


    plt.savefig(save_fig_path)


    #=====================
    fig_subp = [1, 1]
    lg = general_set_all_rcParams(figure_subp=fig_subp, legend_framealpha=0.4,
                                  font_family=['Arial'])
    fig, ax2 = plt.subplots(fig_subp[0], fig_subp[1], sharex = True)


    xi = latoms[0].positions[:, 2]
    for i in np.arange(njobs):
        temp = np.array([ xi, dpos3[i, :] ])
        ind = np.argsort(temp[0, :])
        temp2 = temp[:,ind]

        ax2.plot(temp2[0, :], temp2[1, :], '-o', label=f'{float(jobn[i]):.2f}')
    
    ax2.legend(loc='lower center', ncol=5,fontsize = 10, bbox_to_anchor=(0.5, 0.0))
    ax2.set_xlabel('Atom positions in $x_3$ ($\\mathrm{\\AA}$)')
    ax2.set_ylabel('Displacement $u_3$ ($\\mathrm{\\AA}$)')

    plt.savefig(save_u3_fig_path)
