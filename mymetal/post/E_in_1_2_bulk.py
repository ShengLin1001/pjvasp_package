"""
E_in_1_2_bulk post-processing submodule.

This module provides functions to analyze and visualize the results of E_in_1_2 bulk calculations,
including generating contour and profile plots, and writing results to text files.

Functions:
    - post_E_in_1_2_bulk: Main function to post-process E_in_1_2 bulk calculations.
    - my_plot_E_in_1_2_bulk: Generate 2D contour and profile plots for E_in_1_2 bulk analysis.
    - my_write_E_in_1_2_bulk: Write E_in_1_2 bulk results to formatted text file.
    - my_read_E_in_1_2_bulk: Read E_in_1_2 bulk results from formatted text file.
    - get_E_in_1_2_bulk_ref: Match this tree to a reference tree case by case, one energy zero.
    - write_E_in_1_2_bulk_ref: Write the case-by-case comparison with the reference tree.

Change log:
    - Writed by Pei on 2025-11-18
"""


import numpy as np
import pandas as pd
from ase import Atoms
from mymetal.post.general import get_structure_info
from mymetal.io.general import general_read, general_write
from mymetal.universal.plot.workflow import my_plot_E_in_1_2_bulk, my_plot_E_in_1_2_bulk_ref
from mymetal.universal.print.print import fail

def post_E_in_1_2_bulk(jobn: list = None, Etot: list = None, 
                       atoms_ref: Atoms = None,
                       latoms: list = None,
                       save_fig_path: str = 'p_post_E_in_1_2_bulk.pdf', save_fig_path2: str = 'p_post_E_in_1_2_bulk2.pdf',save_txt_path: str = 'p_post_E_in_1_2_bulk.txt',
                       ref_txt: str = None,
                       save_ref_fig_path: str = 'p_post_E_in_1_2_bulk_ref.pdf',
                       save_ref_txt_path: str = 'p_post_E_in_1_2_bulk_ref.txt'):
    """Post-process E_in_1_2 bulk calculations and generate analysis plots.
    
    Args:
        jobn: List of job names containing a1 and a2 lattice parameters
        Etot: List of total energies in eV
        atoms_ref: Reference ASE Atoms object
        latoms: List of atomic structures
        save_fig_path: Path for saving contour plot
        save_fig_path2: Path for saving profile plots  
        save_txt_path: Path for saving results data
        ref_txt: p_post_E_in_1_2_bulk.txt of a reference tree (e.g. the DFT one). Its points
            are overlaid on the profile plots, and the common cases go to save_ref_txt_path /
            save_ref_fig_path (parity). None: no comparison.
        save_ref_fig_path: Path for the parity plot against the reference
        save_ref_txt_path: Path for the case-by-case comparison table
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
    lenergy_ref, lca_ref = None, None
    if ref_txt is not None:
        df_cmp, lenergy_ref, lca_ref = get_E_in_1_2_bulk_ref(df, ref_txt)
        write_E_in_1_2_bulk_ref(df_cmp, ref_txt, save_ref_txt_path)
        my_plot_E_in_1_2_bulk_ref(df_cmp, if_save=True, savefile=save_ref_fig_path)
    eq_info = my_plot_E_in_1_2_bulk(la1 = df["a1"].astype(float).tolist(),
                                la2 = df["a2"].astype(float).tolist(),
                                lenergy = df["Etot_eV_per_atom"].astype(float).tolist(),
                                lca = df["lca"].astype(float).tolist(),
                                save_fig_path= save_fig_path,
                                save_fig_path2= save_fig_path2,
                                lenergy_ref = lenergy_ref, lca_ref = lca_ref
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




def get_E_in_1_2_bulk_ref(df: pd.DataFrame = None, ref_txt: str = None) -> tuple:
    """Match this tree to a reference tree case by case (same a1, a2) on one energy zero.

    The zero (anchor) is the common case of lowest reference energy, so the two trees line up
    even when their grids differ or their absolute energies do not compare (DFT vs a potential).

    Args:
        df: This tree, columns a1, a2, Etot_eV_per_atom, lca (as built in post_E_in_1_2_bulk).
        ref_txt: p_post_E_in_1_2_bulk.txt of the reference tree.

    Returns:
        tuple: ``(df_cmp, lenergy_ref, lca_ref)``. df_cmp holds the common cases: a1, a2,
        de, de_ref (E - E(anchor), meV/atom), dde = de - de_ref, ca, ca_ref (c/a),
        dca_pct = (ca/ca_ref - 1) * 100. lenergy_ref / lca_ref follow the rows of df: reference
        eV/atom shifted onto this tree's scale at the anchor, and c/a; NaN where the reference
        has no such case.
    """
    _, _, dfr = my_read_E_in_1_2_bulk(filename=ref_txt)
    # round: the case names carry 8 decimals, but a float read back may differ in the last bit
    dict_key = {'k1': lambda d: d['a1'].astype(float).round(6),
                'k2': lambda d: d['a2'].astype(float).round(6)}
    dfr = dfr.assign(**dict_key)
    if dfr.duplicated(['k1', 'k2']).any():
        fail('reference %s has repeated (a1, a2) cases' % ref_txt)
    dfm = df.assign(**dict_key).merge(dfr[['k1', 'k2', 'Etot_eV_per_atom', 'lca']],
                                      on=['k1', 'k2'], how='left', suffixes=('', '_ref'))
    energy = dfm['Etot_eV_per_atom'].astype(float)
    energy_ref = dfm['Etot_eV_per_atom_ref'].astype(float)
    mask = energy_ref.notna()
    if not mask.any():
        fail('no case of %s has the (a1, a2) of a case of this tree' % ref_txt)
    i0 = energy_ref[mask].idxmin()
    dfm['de'] = (energy - energy[i0]) * 1e3
    dfm['de_ref'] = (energy_ref - energy_ref[i0]) * 1e3
    dfm['dde'] = dfm['de'] - dfm['de_ref']
    dfm['ca'] = dfm['lca'].astype(float)
    dfm['ca_ref'] = dfm['lca_ref'].astype(float)
    dfm['dca_pct'] = (dfm['ca'] / dfm['ca_ref'] - 1) * 100
    df_cmp = (dfm.loc[mask, ['k1', 'k2', 'de', 'de_ref', 'dde', 'ca', 'ca_ref', 'dca_pct']]
              .rename(columns={'k1': 'a1', 'k2': 'a2'}).reset_index(drop=True))
    lenergy_ref = (energy_ref - energy_ref[i0] + energy[i0]).tolist()
    return df_cmp, lenergy_ref, dfm['ca_ref'].tolist()


def write_E_in_1_2_bulk_ref(df_cmp: pd.DataFrame = None, ref_txt: str = None,
                            filename: str = 'p_post_E_in_1_2_bulk_ref.txt') -> None:
    """Write the case-by-case comparison of get_E_in_1_2_bulk_ref (readable by
    ``pd.read_csv(filename, sep=r'\\s+', comment='#')``).

    Args:
        df_cmp: Common cases from get_E_in_1_2_bulk_ref.
        ref_txt: The reference file, recorded in the header.
        filename: Output file path.
    """
    anchor = df_cmp.loc[df_cmp['de_ref'].abs().idxmin()]
    row_e = df_cmp.loc[df_cmp['dde'].abs().idxmax()]
    row_c = df_cmp.loc[df_cmp['dca_pct'].abs().idxmax()]
    summary = ('npt %d  dde MAE %.3f max %.3f at (%.4f, %.4f)  dca_pct MAE %.3f max %.3f at '
               '(%.4f, %.4f)' % (len(df_cmp), df_cmp['dde'].abs().mean(), row_e['dde'],
                                 row_e['a1'], row_e['a2'], df_cmp['dca_pct'].abs().mean(),
                                 row_c['dca_pct'], row_c['a1'], row_c['a2']))
    with open(filename, 'w', encoding='utf-8') as f:
        f.write('# this tree vs reference %s, matched by (a1, a2)\n' % ref_txt)
        f.write('# de = E - E(anchor) (meV/atom); anchor = common case of lowest reference '
                'energy: a1 %.6f a2 %.6f\n' % (anchor['a1'], anchor['a2']))
        f.write('# dde = de - de_ref (meV/atom); ca = c/a; dca_pct = (ca/ca_ref - 1) * 100\n')
        f.write('# %s\n' % summary)
        f.write(df_cmp.to_string(index=False, float_format=lambda x: '%.6f' % x) + '\n')
    print("✅ comparison with the reference written to '%s'" % filename)
    print('📊 %s' % summary)
