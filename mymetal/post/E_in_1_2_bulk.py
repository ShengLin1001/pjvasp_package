"""
E_in_1_2_bulk post-processing submodule.

This module provides functions to analyze and visualize the results of E_in_1_2 bulk calculations,
including generating contour and profile plots, and writing results to text files.

Functions:
    - post_E_in_1_2_bulk: Main function to post-process E_in_1_2 bulk calculations.
    - my_plot_E_in_1_2_bulk: Generate 2D contour and profile plots for E_in_1_2 bulk analysis.
    - my_write_E_in_1_2_bulk: Write E_in_1_2 bulk results to formatted text file.
    - my_read_E_in_1_2_bulk: Read E_in_1_2 bulk results from formatted text file.

Change log:
    - Writed by Pei on 2025-11-18
"""


import numpy as np
import pandas as pd
from ase import Atoms
from mymetal.post.general import get_structure_info
from mymetal.io.general import general_read, general_write
from mymetal.universal.plot.workflow import my_plot_E_in_1_2_bulk

def post_E_in_1_2_bulk(jobn: list = None, Etot: list = None, 
                       atoms_ref: Atoms = None,
                       latoms: list = None,
                       save_fig_path: str = 'p_post_E_in_1_2_bulk.pdf', save_fig_path2: str = 'p_post_E_in_1_2_bulk2.pdf',save_txt_path: str = 'p_post_E_in_1_2_bulk.txt'):
    """Post-process E_in_1_2 bulk calculations and generate analysis plots.
    
    Args:
        jobn: List of job names containing a1 and a2 lattice parameters
        Etot: List of total energies in eV
        atoms_ref: Reference ASE Atoms object
        latoms: List of atomic structures
        save_fig_path: Path for saving contour plot
        save_fig_path2: Path for saving profile plots  
        save_txt_path: Path for saving results data
    """
    df = pd.DataFrame(jobn, columns = ["jobn"])
    pattern = r"a1-([\d\.]+)-a2-([\d\.]+)"
    natoms = len(atoms_ref)
    df["a1"] = df["jobn"].str.extract(pattern)[0]
    df["a2"] = df["jobn"].str.extract(pattern)[1]
    df["Etot_eV"] = Etot # eV
    df["Etot_eV_per_atom"] = df["Etot_eV"] / natoms
    lcvectors, lrvectors, lcellpars, lca = get_structure_info(latoms)
    df["lca"] = lca
    #print(df)
    eq_info = my_plot_E_in_1_2_bulk(la1 = df["a1"].astype(float).tolist(),
                                la2 = df["a2"].astype(float).tolist(),
                                lenergy = df["Etot_eV_per_atom"].astype(float).tolist(),
                                lca = df["lca"].astype(float).tolist(),
                                save_fig_path= save_fig_path,
                                save_fig_path2= save_fig_path2
                                )
    
    my_write_E_in_1_2_bulk( df = df,
                            natoms = natoms, 
                           eq_info = eq_info,
                           filename = save_txt_path)


def my_write_E_in_1_2_bulk(df, 
                           natoms,
                           eq_info,
                            filename='p_post_E_in_1_2_bulk.txt'):
    """Write E_in_1_2 bulk results to formatted text file.
    
    Args:
        df: DataFrame containing calculation results
        natoms: Number of atoms in system
        eq_info: Equilibrium parameters (a1, a2, energy)
        filename: Output file path
    """

    (eq_a1, eq_a2, eq_energy) = eq_info

    with open(filename, 'w') as f:
        f.write('# E_in_1_2_bulk post-processing results\n\n')

        f.write(f'{"natoms":>16}\n')
        f.write(f'{natoms:16d}\n\n')

        f.write(f'{"eq infos":>16}\n')
        f.write(f'{"eq_a1":>16} {"eq_a2":>16} {"eq_Etot":>16}\n')
        f.write(f'{eq_a1:16.8f} {eq_a2:16.8f} {eq_energy:16.8f}\n\n')

        # f.write(f'{"jobn":>16} {"Etot (eV)":>16} {"lca":>16}\n')
        # for j, e, ca in zip(jobn, Etot, lca):
        #     f.write(f'{str(j):>16} {e:16.8f} {ca:16.8f}\n')
    general_write(filename=filename, if_append=True, dfc = df, if_write_col_num=True, if_write_row_num=False)
    
    print(f"✅ E_in_1_2_bulk data written to '{filename}'")
    my_read_E_in_1_2_bulk(filename=filename)

def my_read_E_in_1_2_bulk(filename: str = 'p_post_E_in_1_2_bulk.txt'):
    """Read E_in_1_2 bulk results from formatted text file.
    
    Args:
        filename: Input file path
        
    Returns:
        tuple: (natoms, eq_info, df) containing system info and data

    Note:
        df is a DataFrame with the following columns:
            - jobn: Job name
            - a1: Lattice parameter a1
            - a2: Lattice parameter a2
            - Etot_eV: Total energy in eV
            - Etot_eV_per_atom: Total energy per atom in eV
            - lca: Lattice constant a
        eq_info is a tuple (eq_a1, eq_a2, eq_energy)
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    i = 0
    while lines[i].startswith('#') or lines[i].strip() == '':
        # print(lines[i].strip())
        i += 1
    # natoms header line now
    i += 1
    natoms = int(lines[i].strip())
    i += 4  
    eq_info = np.array(list(map(float, lines[i].split())))
    # hearder is in the 9-line (idx from 0), but before has 1 comment line and 3 blank lines
    df = general_read(filepath=filename, header_row = 5)

    return natoms, eq_info, df


