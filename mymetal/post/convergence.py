from ase.io.vasp import read_vasp
import numpy as np
import os
from myvasp import vasp_func as vf 
from mymetal.universial.plot.plot import my_plot_convergence
from mymetal.post.general import my_sort

def post_convergence(dirsurf: str = 'y_convergence', dirlists = ['y_convergence_encuts', 'y_convergence_kpoints'], refcontcar: str='./y_full_relax/CONTCAR'):
    """
    Post-process VASP convergence calculations and generate plots and summaries.

    Args:
        dirsurf (str): Root directory containing convergence subdirectories.
        dirlists (list): List of subdirectories to analyze (e.g., encuts or kpoints).
        refcontcar (str): Path to reference CONTCAR file for determining atom count.

    Raises:
        FileNotFoundError: If `dirsurf` or any required subdirectory is missing.

    Side Effects:
        - Executes `yin_vasp_univ_post` to extract VASP outputs.
        - Saves convergence results to 'p_post_convergence.txt'.
        - Generates annotated convergence plots.
    """
    myroot = os.getcwd()
    refcontcar = os.path.join(myroot, refcontcar)
    dirsurf    = os.path.join(myroot, dirsurf)

    if os.path.isdir(dirsurf):
        print(f"✅ Directory {dirsurf} exists.")
    else:
        raise FileNotFoundError(f"❌ Directory {dirsurf} does not exist. Please run the convergence calculations first.")


    atoms_ref = read_vasp(refcontcar)
    natoms    = atoms_ref.get_positions().shape[0]


    for dir in dirlists:
        os.chdir(dirsurf)
        dir = os.path.join(dirsurf, dir)
        if os.path.isdir(dir):
            print(f"✅ Directory {dir} exists.")
        else:
            continue
        os.chdir(dir)
        # general post
        # os.system("yin_vasp_univ_post > /dev/null 2>&1")

        # read post data
        jobn, Etot, Eent, pres = vf.vasp_read_post_data() # list, array, array, array | str, eV, eV, kB

        encuts  = False
        kpoints = False
        reverse = False
        # re-read and sort jobn
        x  = np.array([])
        # change encuts, kpoints, reverse tag
        if dir.endswith('y_convergence_encuts'):
            x = np.array([float(i) for i in jobn])
            encuts  = True
        elif dir.endswith('y_convergence_kpoints'):
            x = np.array([[float(num) for num in i.split('-')] for i in jobn])
            kpoints = True

        x_sorted, Etot_sorted, Eent_sorted, pres_sorted = my_sort(x, Etot, Eent, pres, reverse)

        # my_plot_convergence(x_sorted , Etot_sorted, encuts=encuts, kpoints=kpoints)
        # my_plot_convergence(x_sorted , Etot_sorted, encuts=encuts, kpoints=kpoints, if_difference=True, if_text=True,)
        my_plot_convergence(x_sorted , Etot_sorted / natoms, encuts=encuts, kpoints=kpoints, if_difference=True, if_text=True, if_mask=True, 
                            mask_condition=lambda y: np.abs(y) > 1, if_save=True)
        my_write_convergence(natoms, jobn, Etot_sorted)

def my_write_convergence(natoms, jobn, Etot):
    """
    Write convergence results to a formatted text file.

    Args:
        natoms (int): Number of atoms in the system.
        jobn (list or array): List of job identifiers (e.g., cutoff or k-points).
        Etot (list or array): Total energy values (eV).
    
    Output:
        - Writes to 'p_post_convergence.txt' in current working directory.
    """
    f = open('p_post_convergence.txt','w')
    f.write('# Convergence test for encuts and kpoints: \n' )


    f.write('\n%16s\n' \
        %('natoms') )

    f.write('%16d\n' \
        %( natoms) )

    f.write('\n%16s %16s \n' \
        %('jobn', 'Etot (eV)') )

    for i in np.arange(len(jobn)):
        f.write('%16s %16.8f\n' \
            %(jobn[i], Etot[i]) )

    f.close()  

def my_read_convergence(file: str = 'p_post_convergence.txt') -> tuple:
    """
    Read convergence data from the result file.

    Args:
        file (str): Path to the result file (default is 'p_post_convergence.txt').

    Returns:
        tuple:
            natoms (int): Number of atoms in the system.
            jobn (np.ndarray): Array of job identifiers (1D or 2D).
            Etot (np.ndarray): Array of total energy values (eV).
    """

    with open(file, 'r') as f:
        lines = f.readlines()

    # Extract number of atoms
    for i, line in enumerate(lines):
        if 'natoms' in line:
            natoms = int(lines[i + 1].strip().split()[0])
            break

    # Initialize data reading flags and containers
    data_started = False
    jobn_raw = []
    Etot = []

    # Parse the data block following 'Etot'
    for i, line in enumerate(lines):
        if 'Etot' in line:
            data_started = True
            continue
        if data_started:
            if line.strip() == "":
                break  # Stop at the first blank line after the data block
            parts = line.strip().split()
            jobn_raw.append(parts[0])
            Etot.append(float(parts[1]))

    # Convert job identifiers to numeric format
    if any('-' in s for s in jobn_raw):
        jobn = np.array([[float(x) for x in s.split('-')] for s in jobn_raw])
    else:
        jobn = np.array([float(s) for s in jobn_raw])

    return natoms, jobn, np.array(Etot)

