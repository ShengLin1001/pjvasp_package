"""
density submodule

This submodule provides functions for calculating the density of an atomic structure. It includes a function to compute
the density based on the number of atoms and the volume of the structure.

Functions:
    - cal_density: Calculate the density of an atomic structure.

"""

from ase import Atoms

def cal_density(atoms: Atoms) -> float:
    volume = atoms.get_volume()  # in Å^3
    natoms = len(atoms)
    density = natoms / volume  # in atoms/Å^3
    return density