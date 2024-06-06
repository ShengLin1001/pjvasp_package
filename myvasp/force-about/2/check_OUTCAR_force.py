#！/usr/bin/env python3



import argparse as ap
import numpy as np
import re
from pathlib import Path


parser = ap.ArgumentParser(add_help=True,
                    formatter_class=ap.ArgumentDefaultsHelpFormatter,
                    description="""
                    Author:  Dr. Huan Wang, 
                    Version: v1.0,
                    Date:    May 11, 2022""")
parser.add_argument("-o",
                    metavar="<VASP OUTCAR file>",
                    type=Path,
                    help="OUTCAR file",
                    default=(Path.cwd() / "OUTCAR"),
                    )
args = parser.parse_args()


def check_OUTCAR(outfile):
    """
    This function checks whether the OUTCAR file exists.

    Parameters
    ----------
    outfile : Path
        The path to the OUTCAR file.

    Raises
    ------
    SystemExit
        If the OUTCAR file does Not exist, print the warnning.

    Returns
    -------
    None.

    """
    if not outfile.is_file():
        drawline = "-" * 79
        NoOutcar = "OUTCAR file does NOT exist! Please check your directory."
        warn_no_file = "\n".join((drawline, NoOutcar, drawline))
        raise SystemExit(warn_no_file)


def grab_info(outfile):
    """ This function grabs two data, e.g. the number of atoms and 
    the convergence criteria of force, from the OUTCAR file. Then,
    collects data from the position and total-force table.

    Parameters
    ----------
    outfile : Path
        The path to the OUTCAR file.

    Returns
    -------
    num_atom    : int
        The number of atoms in the system.
    ediffg      : float
        The convergence criteria of force.
    position    : numpy array
        The 3D data array of the x, y, z posistion in unit of angstrom
    total force : numpy array
        The 3D data array of the Total Force components in unit of eV/Angst
    """
    data = []
    table_list = []
    patten_NumAtom = r"\s+\w+\s+NIONS =\s+(\d+)"
    patten_EdiffG  = r"\s+EDIFFG =(\s[-+]?\.\w+[-+]?\d+)"
    patten_drift = r"\s+total drift:\s+"
    
    with open(outfile, "r") as fo:
        for ind, line in enumerate(fo):
            table_list.append(line.strip().split())
            
            if re.search(patten_NumAtom, line):
                num_atom = int(re.search(patten_NumAtom, line).group(1))
                
            elif  re.search(patten_EdiffG, line):
                ediffg = float(re.search(patten_EdiffG, line).group(1))
                
            elif re.search(patten_drift, line):
                data += table_list[-num_atom-2:ind-1]

    return num_atom, ediffg, \
           np.asfarray(data).reshape(-1,num_atom,6)[:,:,:3],\
           np.asfarray(data).reshape(-1,num_atom,6)[:,:,3:]


def calcuate_force(force, criteria, NumAtom):
    """ This function calculates the force for each atom.
    
    Parameters
    ----------
    force : numpy 3D array
        The force components of each atoms at each step.
    criteria : float
        The convergence criteria of force.
    NumAtom : int
        number of atoms.

    Returns
    -------
    Boolean array
        The Boolean array that True represents the force of the atom 
    does NOT converged.

    """
    F = np.sqrt(np.power(force, 2).sum(axis=-1))
    return F > np.abs(criteria)
    

def main():
    """ work flow
    1) check whether the OUTCAR file exists or not;
    2) grab the number of atoms, the convergence criteria, and the force 
       components of each atom from the OUTCAR file;
    3) calculate the force for each atom by using numpy;
    4) print the steps and atoms with the force that larger than 
       the convergence criteria.
    """

    outcar = args.o
    
    fmt = """There are {:} atoms in this system and {:} steps for calculatons.
          In the step {:}, the force of atom {:} does NOT converge.\n"""
    

    check_OUTCAR(outcar)
    num_atom, ediffg, position, total_force = grab_info(outcar)
    check = calcuate_force(total_force, ediffg, num_atom)
    for step, col in enumerate(check):
        for atom, convergence in enumerate(col):
            if convergence == True:
                print(fmt.format(len(col), len(check), step+1, atom+1))



if __name__ == "__main__":
    main()