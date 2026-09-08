"""
convergence post-processing submodule

This module provides functions to post-process VASP convergence calculations,
generate plots, and summarize results in formatted text files.

Functions:
    - get_axis_of_dir: Map a y_convergence_* subtree to its scan axis.
    - parse_case_names: Turn case directory names into a numeric x axis.
    - post_convergence: Main function to process convergence data and generate outputs.
    - my_write_convergence: Write convergence results to a formatted text file.
    - my_read_convergence: Read convergence data from the result file.

Change log:
    - Written by J. P. on 2025.11.04.
    - 2026.09.01: auto-detect which y_convergence_* subtrees exist instead of taking a
      fixed list, add the klength / ediff / kpar_ncore axes, dispatch kpar_ncore to
      mymetal.post.kpar_ncore, and read the self-describing case names
      (encut_400, kpoints_3-3-3, klen_81, ediff_1E-6) with a fallback to the older
      bare-value names (400, 3-3-3) so finished trees stay readable.
"""

# Written by J. P.
# 2025.11.04

from ase.io.vasp import read_vasp
import numpy as np
import os
from myvasp import vasp_func as vf 
from mymetal.universal.plot.workflow import my_plot_convergence
from mymetal.post.general import my_sort
from mymetal.post.kpar_ncore import post_kpar_ncore

# One row per scan axis: the y_convergence_<axis> subtree name is the key the
# post-processor dispatches on, and the label pair is what the figure gets. Adding an
# axis means adding a row here and one in mymetal.build.workflow.convergence.
DICT_AXIS_PREFIX = {'encuts': 'encut_', 'kpoints': 'kpoints_',
                    'klength': 'klen_', 'ediff': 'ediff_'}
DICT_AXIS_LABEL = {
    'encuts':  ['Energy cutoff (eV)', r'$E-E^{\text{Ref}}$ (meV per atom)'],
    'kpoints': ['K-point grid', r'$E-E^{\text{Ref}}$ (meV per atom)'],
    'klength': [r'Automatic $k$-mesh length $A$ (-)', r'$E-E^{\text{Ref}}$ (meV per atom)'],
    'ediff':   ['EDIFF (eV)', r'$E-E^{\text{Ref}}$ (meV per atom)'],
}
# The parallel axis is a timing benchmark, not a convergence curve, so it is dispatched
# to its own reader/plotter instead of going through the curve path below.
AXIS_KPAR_NCORE = 'kpar_ncore'
# Reference end of the sorted axis. Every curve is drawn as E - E_ref, and my_plot_
# convergence takes the LAST point as the reference -- so each axis has to be sorted so
# that its most accurate setting lands last. Denser k mesh and higher ENCUT are more
# accurate ascending; a TIGHTER EDIFF is a SMALLER number, so that one sorts descending.
DICT_AXIS_REVERSE = {'encuts': False, 'kpoints': False, 'klength': False, 'ediff': True}


def get_axis_of_dir(dirname: str = None) -> str:
    """Map a ``y_convergence_*`` subtree name to its scan axis.

    Args:
        dirname (str): Directory name directly under ``y_convergence``.

    Returns:
        str: Axis name, or ``None`` if the directory is not a scan subtree.
    """
    for axis in list(DICT_AXIS_PREFIX) + [AXIS_KPAR_NCORE]:
        if dirname == 'y_convergence_%s' % axis:
            return axis
    return None


def parse_case_names(ljobn: list = None, axis: str = None):
    """Turn case directory names into the numeric x axis of one scan.

    Accepts both the self-describing names written since 2026-09 (``encut_400``,
    ``kpoints_3-3-3``, ``klen_81``, ``ediff_1E-6``) and the bare-value names the
    older bash workflow wrote (``400``, ``3-3-3``), so finished trees stay readable.

    Args:
        ljobn (list): Case directory names, as returned by vasp_read_post_data.
        axis (str): Scan axis these cases belong to.

    Returns:
        np.ndarray: 1D values, or one ``[n1, n2, n3]`` row per case for ``kpoints``.
    """
    prefix = DICT_AXIS_PREFIX[axis]
    lraw = [jobn[len(prefix):] if jobn.startswith(prefix) else jobn
            for jobn in ljobn]
    if axis == 'kpoints':
        return np.array([[float(num) for num in raw.split('-')] for raw in lraw])
    return np.array([float(raw) for raw in lraw])


def post_convergence(dirsurf: str = 'y_convergence', dirlists: list = None,
                     refcontcar: str = './y_full_relax/CONTCAR'):
    """
    Post-process a VASP parameter scan and generate plots and summaries.

    Which axes exist is read off the directory tree rather than passed in: a run that
    only scanned ENCUT has only y_convergence_encuts, and asking for the others would
    be an error rather than a no-op. Pass ``dirlists`` to restrict the set by hand.

    Args:
        dirsurf (str): Root directory containing the y_convergence_* subtrees.
        dirlists (list): Subtrees to analyze; ``None`` auto-detects every known one.
        refcontcar (str): Reference CONTCAR/POSCAR used only for the atom count.

    Raises:
        FileNotFoundError: If `dirsurf` is missing.

    Side Effects:
        - Executes `yin_vasp_univ_post` in each curve subtree.
        - Writes p_post_convergence.{txt,pdf} into each curve subtree.
        - Writes p_post_kpar_ncore.{txt,pdf} into the parallel subtree.
    """
    ### check ------------------------------------------------------------------
    myroot = os.getcwd()
    refcontcar = os.path.join(myroot, refcontcar)
    dirsurf = os.path.join(myroot, dirsurf)

    if os.path.isdir(dirsurf):
        print(f"✅ Directory {dirsurf} exists.")
    else:
        raise FileNotFoundError(f"❌ Directory {dirsurf} does not exist. Please run the convergence calculations first.")

    # A y_full_relax built from a snapshot rather than an actual relaxation has only a
    # POSCAR; the atom count is all that is read here, so either file will do.
    if not os.path.isfile(refcontcar):
        refposcar = os.path.join(os.path.dirname(refcontcar), 'POSCAR')
        if not os.path.isfile(refposcar):
            raise FileNotFoundError(
                f"❌ neither {refcontcar} nor {refposcar} exists; the atom count "
                f"cannot be determined.")
        print(f"⚠️  {refcontcar} not found, using POSCAR")
        refcontcar = refposcar
    atoms_ref = read_vasp(refcontcar)
    natoms    = atoms_ref.get_positions().shape[0]
    print(f"📄 Reference   : {refcontcar} ({natoms} atoms)")

    if dirlists is None:
        dirlists = sorted(name for name in os.listdir(dirsurf)
                          if get_axis_of_dir(name) is not None
                          and os.path.isdir(os.path.join(dirsurf, name)))
    if not dirlists:
        print("⚠️  no y_convergence_* subtree found; nothing to post-process")
        return
    print(f"🔭 Scan axes   : {len(dirlists)} ({', '.join(dirlists)})")
    # to here: the tree exists and we know which axes it actually carries

    ### main -------------------------------------------------------------------
    ndone = 0
    lskip = []
    for dirn in dirlists:
        os.chdir(dirsurf)
        axis = get_axis_of_dir(dirn)
        dir = os.path.join(dirsurf, dirn)
        if axis is None or not os.path.isdir(dir):
            print(f"⚠️  skipping {dirn}: not a known scan subtree")
            lskip.append(dirn)
            continue
        print(f"\n================ 📁 {dirn}  (axis: {axis})")
        os.chdir(dir)

        if axis == AXIS_KPAR_NCORE:
            # Timing benchmark: its own reader, table and two-panel figure, and it
            # runs pei_vasp_univ_post itself.
            post_kpar_ncore(path_workflow=dir)
            ndone += 1
            continue

        # general post
        os.system("yin_vasp_univ_post")

        # read post data
        jobn, Etot, Eent, pres = vf.vasp_read_post_data() # list, array, array, array | str, eV, eV, kB
        x = parse_case_names(jobn, axis)

        # Sorted here rather than through my_sort because jobn has to travel with the
        # same permutation: the old code wrote unsorted jobn against sorted energies,
        # which mislabels every row as soon as the directory order is not numeric order.
        sort_idx = np.argsort(x) if x.ndim == 1 else np.argsort(x[:, 0])
        if DICT_AXIS_REVERSE[axis]:
            sort_idx = sort_idx[::-1]
        x_sorted = x[sort_idx]
        Etot_sorted = np.array(Etot)[sort_idx]
        ljobn_sorted = [jobn[index] for index in sort_idx]

        my_plot_convergence(x_sorted, Etot_sorted / natoms,
                            kpoints=(axis == 'kpoints'),
                            laxis_label=DICT_AXIS_LABEL[axis],
                            if_logx=(axis == 'ediff'),
                            if_difference=True, if_text=True, if_mask=True,
                            mask_condition=lambda y: np.abs(y) > 1, if_save=True)
        my_write_convergence(natoms, ljobn_sorted, Etot_sorted)
        ndone += 1
        print(f"✅ {dirn}: p_post_convergence.txt / p_post_convergence.pdf written")
    os.chdir(myroot)
    # to here: every detected axis has its own table and figure beside its cases

    print("\n================ 📊 Summary")
    print(f"axes processed={ndone}    skipped={len(lskip)}")
    for dirn in lskip:
        print(f"❌ {dirn}")
    print("🎉 Done.")


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

