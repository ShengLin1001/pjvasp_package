
# This page is taken from https://github.com/BinglunYin/vasp_utils/blob/master/vasp_workflow_bulk/yin_vasp_plot_cij_energy.py
# Changed to fit the lammps output by J. P. 
# 2025.11.04

import numpy as np
from myvasp import vasp_func as vf 
import os
import matplotlib.pyplot as plt
from mymetal.universal.plot.general import general_set_all_rcParams
from ase import Atoms
from ase.io import read

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
    atoms_ref = read(refcontcar, format='lammps-data')
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

    plot_cij_energy(ldata, save_fig_path, save_txt_path)

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

def plot_cij_energy(ldata: list = None, 
                    save_fig_path: str = './y_post_cij_energy.pdf',
                    save_txt_path: str = './y_post_cij_energy.txt',):
    """
    Fit energy-strain curves to quadratic functions and plot both energy and stress.

    Args:
        ldata (list): List of dictionaries containing deformation data.
        save_fig_path (str): Path to save the plot figure.
        save_txt_path (str): Path to save Cij and derived elastic properties.

    Returns:
        None
    """
    fig_subp = [5, 2]
    lg = general_set_all_rcParams(figure_subp=fig_subp, legend_framealpha=0.4,
                                  font_family=['Arial'])
    fig, axes = plt.subplots(fig_subp[0], fig_subp[1], sharex = True)

    # plot energy-strain 

    xi = np.linspace(min(ldata[0]['e']), max(ldata[0]['e']), 1000)

    p0_tot = np.array([])     # for Cij
    e0_tot = np.array([])     # to check shift 

    for i in np.arange(len(ldata)):
        e  = ldata[i]['e']
        ed = ldata[i]['ed'] 
        s  = ldata[i]['s'] 

        ax = axes[i, 0]
        ax.plot(e, ed, 'o')

        param = np.polyfit(e, ed, 2)
        fun = np.poly1d(param)
        yi = fun(xi)
        p0_tot = np.append(p0_tot, param[0]) 

        e0 = -param[1]/(2*param[0])
        e0_tot = np.append(e0_tot, e0 ) 

        ax.plot(xi, yi, '-C1')

        str1 = '$\\epsilon_0$ = %.6f' %(e0)
        ymin, ymax = ax.get_ylim()
        ax.text(0, ymax*0.85, str1)

    C11 = p0_tot[0]*2
    C12 = p0_tot[1]-C11
    C33 = p0_tot[3]*2
    C13 = p0_tot[2]-(C11+C33)/2
    C44 = p0_tot[4]*2

    # plot stress-strain 

    lcolor=['C0','C1','C2','C3','C4','C5']

    llegend=['$\\sigma_{xx}$','$\\sigma_{yy}$','$\\sigma_{zz}$','$\\sigma_{yz}$','$\\sigma_{xz}$','$\\sigma_{xy}$']

    lslope = np.array([
        [ C11    , C12    , C13    , 0, 0, 0 ],          # c11
        [ C11+C12, C12+C11, C13*2  , 0, 0, 0 ],          # c12
        [ C11+C13, C12+C13, C13+C33, 0, 0, 0 ],          # c13
        [ C13    , C13    , C33    , 0, 0, 0 ],          # c33
        [ 0      , 0      , 0      , C44,  0,  0 ],      # c44
    ])
 
    for i in np.arange(len(ldata)):
        e    = ldata[i]['e']
        ed   = ldata[i]['ed'] 
        s    = ldata[i]['s'] 
        iref = ldata[i]['iref'] 

        ax = axes[i, 1]
        for j in np.arange(6):
            ax.plot(e, s[:, j], 'o', color=lcolor[j], label = llegend[j], markersize=6, alpha = 0.4)
            # \sigma = C * \epsilon
            ax.plot(xi, xi * lslope[i,j] + s[iref,j], '-', color=lcolor[j], markersize=6, alpha = 0.4 )

    lg(axes[0,1].legend(loc='upper left', ncol=2))

    for i in np.arange(len(ldata)):
        for j in np.arange(2):
            #axes[i,j].set_xticks(np.arange(-0.002, 0.004, 0.002))
            axes[4,j].set_xlabel('Strain')

        axes[i,0].set_ylabel('Elastic energy density (GPa)')
        axes[i,1].set_ylabel('DFT stress (GPa)')


    str1 = '$C_{11}$, $C_{12}$, $C_{13}$, $C_{33}$, $C_{44}$ (GPa): \n%.2f,   %.2f,   %.2f,   %.2f,   %.2f ' \
        %( C11, C12, C13, C33, C44 )
    ymin, ymax = axes[0,0].get_ylim()
    axes[0,0].text(-0.003, ymax*1.05, str1 )

    plt.savefig(save_fig_path)

    write_cij_energy( np.array([ C11, C12, C13, C33, C44]), save_txt_path )

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