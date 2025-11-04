# Written by J. P.
# 2025.11.04

from ase.io.vasp import read_vasp
import numpy as np
import os
from myvasp import vasp_func as vf 
from mymetal.universal.plot.workflow import my_plot_stretch
from mymetal.post.general import my_read_y_dir_contcar, get_structure_info
from ase.io import read
from mymetal.io.general import general_read

def post_stretch(dirsurf: str = 'y_stretch', refcontcar: str='./y_full_relax/CONTCAR'):
    """
    Post-process VASP stretch calculations.

    Args:
        dirsurf (str): Directory containing stretch calculation outputs. Default 'y_stretch'.
        refcontcar (str): Path to reference CONTCAR file. Default './y_full_relax/CONTCAR'.

    Returns:
        None. Generates plots and writes summary file 'p_post_stretch.txt'.
    """
    myroot     = os.getcwd()
    refcontcar = os.path.join(myroot, refcontcar)
    dirsurf    = os.path.join(myroot, dirsurf)

    if os.path.isdir(dirsurf):
        print(f"✅ Directory {dirsurf} exists.")
    else:
        raise FileNotFoundError(f"❌ Directory {dirsurf} does not exist. Please run 'yin_vasp_run_stretch' first.")


    atoms_ref = read_vasp(refcontcar)
    cell_ref = np.array(atoms_ref.get_cell())
    # Col vector [x, y, z]
    cvectors_ref = np.linalg.norm(cell_ref, axis=0)
    # row vector [a1, a2, a3]
    rvectors_ref = np.linalg.norm(cell_ref, axis=1)
    natoms    = atoms_ref.get_positions().shape[0]

    os.chdir(dirsurf)
    # general post
    os.system("yin_vasp_univ_post")

    # read post data
    jobn, Etot, Eent, pres = vf.vasp_read_post_data() # list, array, array, array | str, eV, eV, kB
    latoms = my_read_y_dir_contcar(os.path.join(dirsurf, "y_dir"))
    lcvectors, lrvectors, lcellpars, lca = get_structure_info(latoms)
    _, stretch_type = get_stretch_type(lcvectors, jobn, cvectors_ref)

    _, (coeffs, extr_x, extr_y, extr_rvectors)  = my_plot_stretch(
        jobn = jobn,
        Etot = Etot,
        natoms = natoms,
        stretch_type = stretch_type,
        rvectors_ref = rvectors_ref,
        lca = lca,
        if_save = True,
        savefile = 'p_post_stretch.pdf',
        dpi = 300,
    )
    my_write_stretch(jobn, Etot, natoms, stretch_type,  rvectors_ref, lca,
                     cell_ref,
                     (coeffs, extr_x, extr_y, extr_rvectors) )
    # jobn, Etot, natoms, stretch_type,  rvectors_ref, lca, cell_ref,\
    # (coeffs, extr_x, extr_y, extr_rvectors) =
    # my_read_stretch('p_post_stretch.txt')


# Post-process LAMMPS stretch calculation
def post_lammps_stretch(post_data_file: str = './y_post_data.txt', refcontcar: str = './CONTCAR', latoms_lammpstrj: str = None,
                        save_fig_path: str = 'p_post_stretch.pdf', save_txt_path: str = 'p_post_stretch.txt'):
    """
    Post-process LAMMPS stretch calculations and extract strain information.

    This function reads stretch data, compares it with a zero-strain reference
    structure, calculates stretch types, and generates related plots and output files.

    Args:
        file (str): Path to the stretch data file (default: './stretch.txt').
        refcontcar (str): Path to the zero-strain reference CONTCAR file (default: './CONTCAR').
        latoms_lammpstrj (str): Path to the LAMMPS trajectory or dump file containing
            strained atomic structures.
        save_fig_path (str): Path to save the generated stretch plot (default: 'p_post_stretch.pdf').
        save_txt_path (str): Path to save the processed stretch data and coefficients
            (default: 'p_post_stretch.txt').

    Returns:
        None

    Side Effects:
        - Generates a plot of stretch vs energy using `my_plot_stretch` (saved as PDF if enabled).
        - Writes processed stretch data and coefficients to an output file via `my_write_stretch`.
    """

    # read stretch data
    jobn, Etot, Eent, pres = vf.vasp_read_post_data(post_data_file)

    # zero-strain reference
    atoms_ref = read(refcontcar, format='lammps-data')
    cell_ref = np.array(atoms_ref.get_cell())
    natoms = len(atoms_ref)
    # Col vector [x, y, z]
    cvectors_ref = np.linalg.norm(cell_ref, axis=0)
    # row vector [a1, a2, a3]
    rvectors_ref = np.linalg.norm(cell_ref, axis=1)

    # list of atoms at different strains
    latoms    = read(latoms_lammpstrj, format='lammps-dump-text', index=':')
    lcvectors, lrvectors, lcellpars, lca = get_structure_info(latoms)
    _, stretch_type = get_stretch_type(lcvectors, jobn, cvectors_ref)

    _, (coeffs, extr_x, extr_y, extr_rvectors)  = my_plot_stretch(
        jobn = jobn,
        Etot = Etot,
        natoms = natoms,
        stretch_type = stretch_type,
        rvectors_ref = rvectors_ref,
        lca = lca,
        if_save = True,
        savefile = save_fig_path,
        dpi = 300,
    )
    my_write_stretch(jobn, Etot, natoms, stretch_type,  rvectors_ref, lca,
                     cell_ref,
                     (coeffs, extr_x, extr_y, extr_rvectors),
                     save_txt_path)


def get_stretch_type(lcvectors: list = None, jobn: list = None, cvectors_ref: np.ndarray = None) -> str:
    """
    Determine stretched lattice axes compared to reference.

    Args:
        lcvectors (np.ndarray): Column vector lengths for each structure.
        jobn (list): Job identifiers corresponding to each structure.
        cvectors_ref (np.ndarray): Reference column vector lengths.

    Returns:
        tuple:
            - if_stretch (np.ndarray): Boolean array indicating which axes are stretched.
            - stretch_type (str): String of stretched axes ('x', 'y', 'z', combinations).
    """
    # More advanced method using strain matrix
    if lcvectors is None or cvectors_ref is None:
        raise ValueError("lcvectors and lengths_ref must be provided.")

    lif_stretch = []
    for cvectors, job in zip(lcvectors, jobn):
        stretch = np.array(cvectors) / np.array(cvectors_ref)
        if_stretch = np.abs(stretch - float(job)) < 1e-5
        lif_stretch.append(if_stretch)
        #print(if_stretch)

    lif_stretch_arr = np.array(lif_stretch)
    # [[True, False, False],
    #  [True, False, False]]
    # => [True, False, False] (if_stretch )
    if_stretch = np.all(lif_stretch_arr, axis=0)
    
    if not np.array_equal(if_stretch, np.array([False, False, False])):
        stretch_axes = ['x', 'y', 'z']
        stretch_type = ''.join([axis for axis, stretched in zip(stretch_axes, if_stretch) if stretched])
        # only those return value: 'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz'
    else:
        raise ValueError("Inconsistent stretch patterns found in the structures.")
    
    return if_stretch, stretch_type

def my_write_stretch(jobn, Etot, natoms, stretch_type,  rvectors_ref, lca, cell_ref,
                     plot_return,
                      file_name='p_post_stretch.txt'):
    """
    Write stretch analysis results and structural data to a text file.

    Args:
        jobn (list): Job identifiers or stretch factors.
        Etot (list): Total energies (eV) corresponding to jobn.
        natoms (int): Number of atoms.
        stretch_type (str): Type of stretch ('x', 'y', 'z', etc.).
        rvectors_ref (array): Reference row vectors.
        lca (list): c/a ratios for each structure.
        cell_ref (array): Reference lattice cell vectors.
        plot_return (tuple): (coeffs, extr_x, extr_y, extr_rvectors) from quadratic fitting.
        file_name (str): Output filename. Default 'p_post_stretch.txt'.

    Returns:
        None.
    """
    (coeffs, extr_x, extr_y, extr_rvectors) = plot_return

    with open(file_name, 'w') as f:
        f.write('# Stretch post-processing results\n\n')

        f.write(f'{"stretch_type":>16}\n')
        f.write(f'{stretch_type:>16s}\n\n')

        f.write(f'{"natoms":>16}\n')
        f.write(f'{natoms:16d}\n\n')

        f.write(f'{"Poly coeffs":>16}\n')
        f.write(f'{"a":>16} {"b":>16} {"c":>16}\n')
        f.write(f'{coeffs[0]:16.8f} {coeffs[1]:16.8f} {coeffs[2]:16.8f}\n\n')

        f.write(f'{"Extr infos":>16}\n')
        f.write(f'{"factor":>16} {"E (eV/atom)":>16}\n')
        f.write(f'{extr_x:16.8f} {extr_y:16.8f}\n\n')

        f.write(f'{"rvectors_ref (A)":>16}\n')
        f.write(f'{rvectors_ref[0]:16.8f} {rvectors_ref[1]:16.8f} {rvectors_ref[2]:16.8f}\n\n')
        
        f.write(f'{"Extr rvector (A)":>16}\n')
        f.write(f'{extr_rvectors[0]:16.8f} {extr_rvectors[1]:16.8f} {extr_rvectors[2]:16.8f}\n\n')
        
        f.write(f'{"cell_ref (A)":>16}\n')
        for vec in cell_ref:
            f.write(f'{vec[0]:16.8f} {vec[1]:16.8f} {vec[2]:16.8f}\n')
        f.write('\n')

        f.write(f'{"jobn":>16} {"Etot (eV)":>16} {"lca":>16}\n')
        for j, e, ca in zip(jobn, Etot, lca):
            f.write(f'{str(j):>16} {e:16.8f} {ca:16.8f}\n')
    
    print(f"✅ Stretch data written to '{file_name}'")

def my_read_stretch(file_name='p_post_stretch.txt'):
    """
    Read stretch analysis results from a formatted text file.

    Args:
        file_name (str): Path to file written by my_write_stretch. Default 'p_post_stretch.txt'.

    Returns:
        tuple:
            - jobn (list): Job identifiers.
            - Etot (np.ndarray): Total energies (eV).
            - natoms (int): Number of atoms.
            - stretch_type (str): Type of stretch.
            - rvectors_ref (np.ndarray): Reference row vectors.
            - lca (np.ndarray): c/a ratios.
            - cell_ref (np.ndarray): Reference cell vectors.
            - plot_return (tuple): (coeffs, extr_x, extr_y, extr_rvectors) from quadratic fitting.
    """
    with open(file_name, 'r') as f:
        lines = f.readlines()

    i = 0
    while lines[i].startswith('#') or lines[i].strip() == '':
        print(lines[i].strip())
        i += 1
    # stretch_type header line now
    i += 1
    stretch_type = lines[i].strip()
    i += 3  
    natoms = int(lines[i].strip())
    i += 4  
    coeffs = np.array(list(map(float, lines[i].split())))
    i += 4
    extr_x, extr_y = map(float, lines[i].split())
    i += 3  
    rvectors_ref = np.array(list(map(float, lines[i].split())))
    i += 3
    extr_rvectors = np.array(list(map(float, lines[i].split())))
    i += 3
    cell_ref = []
    for _ in range(3):
        cell_ref.append(list(map(float, lines[i].split())))
        i += 1
    cell_ref = np.array(cell_ref)
    i += 2
    jobn, Etot, lca = [], [], []
    while i < len(lines):
        line = lines[i].strip()
        if line:
            parts = line.split()
            jobn.append(parts[0])
            Etot.append(float(parts[1]))
            lca.append(float(parts[2]))
        i += 1
    Etot = np.array(Etot)
    lca = np.array(lca)
    
    return jobn, Etot, natoms, stretch_type, rvectors_ref, lca, cell_ref, (coeffs, extr_x, extr_y, extr_rvectors)

