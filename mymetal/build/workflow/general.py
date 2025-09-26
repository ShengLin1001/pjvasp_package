import os
import numpy as np
import shutil
import glob

def cp_vaspfiles(cpfiles = ["CONTCAR", "KPOINTS", "INCAR", "POTCAR", "Y_CONSTR_LATT"], globfiles = ["sub*"], oldir = "y_full_relax/", newdir = "y_full_relax_temp/"):

    for file in cpfiles:
        if os.path.isfile(os.path.join(oldir, file)):
            shutil.copy(os.path.join(oldir, file), newdir)
        else:
            print(f"Warning: {file} not found in {oldir}")
    
    for globfile in globfiles:
        for file in glob.glob(os.path.join(oldir, globfile)):
            shutil.copy(file, newdir)

def rm_i(workdir = None):
    if os.path.exists(workdir):
        ans = input(f"The directory {workdir} already exists. Do you want to delete it? (y/n): ").strip().lower()
        if ans == "y":
            shutil.rmtree(workdir)
            print(f"Deleted the old directory {workdir}")
        else:
            print("Keeping the old directory. No changes made.")

def compare_three_lattices(lattice_before, lattice_after, lattice_after_adjusted):

    print("\nLattices")
    print("Original:\n", lattice_before)
    print("Stretched:\n", lattice_after)
    print("Adjusted:\n", lattice_after_adjusted)

    print("\nVolumes")
    volume_before = np.abs(np.linalg.det(lattice_before))
    volume_after  = np.abs(np.linalg.det(lattice_after))
    volume_after_adjusted  = np.abs(np.linalg.det(lattice_after_adjusted))
    print(f"Original: {volume_before:.6f}")
    print(f"Stretched: {volume_after:.6f}  (Change: {(volume_after - volume_before)/volume_before*100:.4f} %)")
    print(f"Adjusted: {volume_after_adjusted:.6f}  (Change: {(volume_after_adjusted - volume_before)/volume_before*100:.4f} %)")

    print("\nLattice constants (lengths)")
    lengths_before = np.linalg.norm(lattice_before, axis=1)
    lengths_after  = np.linalg.norm(lattice_after, axis=1)
    lengths_after_adjusted  = np.linalg.norm(lattice_after_adjusted, axis=1)
    for i, (l_before, l_after, l_after_adj) in enumerate(zip(lengths_before, lengths_after, lengths_after_adjusted)):
        print(f"a{i+1}: {l_before:.6f} -> {l_after:.6f} -> {l_after_adj:.6f}  (Change: {(l_after - l_before)/l_before*100:.4f} % -> {(l_after_adj - l_before)/l_before*100:.4f} %)")

    print("\nCol vectors (lengths)")
    lengths_before = np.linalg.norm(lattice_before, axis=0)
    lengths_after  = np.linalg.norm(lattice_after, axis=0)
    lengths_after_adjusted  = np.linalg.norm(lattice_after_adjusted, axis=0)
    for i, (l_before, l_after, l_after_adj) in enumerate(zip(lengths_before, lengths_after, lengths_after_adjusted)):
        print(f"c{i+1}: {l_before:.6f} -> {l_after:.6f} -> {l_after_adj:.6f}  (Change: {(l_after - l_before)/l_before*100:.4f} % -> {(l_after_adj - l_before)/l_before*100:.4f} %)")
