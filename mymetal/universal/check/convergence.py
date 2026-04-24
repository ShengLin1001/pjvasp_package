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

def check_and_submit_jobs(path_ydir: Path = None, if_sbatch: bool = False):
    """Check VASP job convergence and optionally resubmit jobs.

    Args:
        path_ydir (str or Path): The parent directory containing job subdirectories.
        if_sbatch (bool, optional): If True, submit jobs when not converged or OUTCAR missing. Default is False.
    """
    path_ydir = Path(path_ydir)
    path_cwd  = os.getcwd()
    for subdir in sorted(path_ydir.iterdir()):
        if subdir.is_dir():
            path_outcar = subdir / "OUTCAR"
            os.chdir(subdir)
            print(f"===================={subdir}")
            if path_outcar.exists():
                ret = os.system("grep -q 'reached required accuracy' OUTCAR")
                if ret == 0:
                    print(f"✅ {subdir.name}")
                else:
                    print(f"❌ {subdir.name}")
                    if if_sbatch:
                        os.system("cp CONTCAR POSCAR && pei_vasp_univ_clean_up_full && sbatch sub.*")
            else:
                print("❌OUTCAR does not exist")
                if if_sbatch:
                    os.system("sbatch sub.*")
            os.chdir(path_cwd)