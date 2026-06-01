"""
convergence.py

This script provides a function to check the convergence status of VASP calculations in subdirectories and optionally resubmit jobs if they are not converged or missing output files.

Functions:
    check_and_submit_jobs(path_ydir: Path, if_sbatch: bool = False): 
        Checks each subdirectory in the given parent directory for VASP job convergence and optionally submits jobs.

Usage:
    Import the function and call it with the target directory and submission option.
"""

import os
from pathlib import Path

def check_and_submit_jobs_in_ydir(path_ydir: Path = None,
                          if_sbatch: bool = False):
    """Check VASP job convergence and optionally resubmit jobs.

    Args:
        path_ydir (str or Path): The parent directory containing job subdirectories.
        if_sbatch (bool, optional): If True, submit jobs when not converged or OUTCAR missing. Default is False.
    """
    path_ydir = Path(path_ydir) if path_ydir else None
    path_cwd  = os.getcwd()

    for subdir in sorted(path_ydir.iterdir()) if path_ydir else []:
        os.chdir(path_cwd)  # For safety, always return to the original directory before processing each subdir
        check_and_submit_jobs_in_subdir(subdir, if_sbatch)


def check_and_submit_jobs_in_subdir(path_subdir: Path = None,
                                   if_sbatch: bool = False):

    if path_subdir.is_dir():
        path_outcar = path_subdir / "OUTCAR"
        os.chdir(path_subdir)
        print(f"===================={path_subdir}")
        if path_outcar.exists():
            ret = os.system("grep -q 'reached required accuracy' OUTCAR")
            if ret == 0:
                print(f"✅ {path_subdir.name}")
            else:
                print(f"❌ {path_subdir.name}")

                ret2 = os.system("grep -q 'please rerun with smaller EDIFF' OUTCAR")
                if ret2 == 0:
                    os.system("pei_vasp_univ_find_and_change -ediff 1e-10")
                    os.system("pei_vasp_univ_find_and_change -algo Normal")
                    
                if if_sbatch:
                    os.system("cp CONTCAR POSCAR && pei_vasp_univ_clean_up_full && sbatch sub.*")
        else:
            print("❌OUTCAR does not exist")
            if if_sbatch:
                os.system("sbatch sub.*")
