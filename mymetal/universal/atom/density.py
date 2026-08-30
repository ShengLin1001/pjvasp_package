"""
density submodule

This submodule provides functions for calculating the density of an atomic structure. It includes a function to compute
the density based on the number of atoms and the volume of the structure.

Functions:
    - cal_density: Calculate the density of an atomic structure.

"""

from ase import Atoms

def cal_density(atoms: Atoms) -> float:
    """Calculate the atomic number density of a structure.

    Args:
        atoms (Atoms): An ``ase.Atoms`` object with a non-degenerate cell.

    Returns:
        float: Number density in atoms/Å³.

    Raises:
        Exception: If the cell is degenerate (volume ≈ 0), ASE raises.
    """
    volume = atoms.get_volume()  # in Å^3
    natoms = len(atoms)
    density = natoms / volume  # in atoms/Å^3
    return density