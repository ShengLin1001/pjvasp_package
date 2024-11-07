"""
fixatom submodule

This submodule provides functions for fixing atoms in a structure. It includes functions for fixing atoms based on
their indices or a boolean mask. These functions are designed to streamline common tasks in materials science simulations
and data handling.

Functions:
    - fixatoms: Fix atoms in a structure based on their indices or a boolean mask.
"""


from ase.constraints import FixAtoms
from ase import Atoms

# take from ase FixAtoms
def fixatoms(atoms: Atoms = None, mask: list = None, indices: list = None) -> Atoms:
    """
    Constrain chosen atoms.

    Parameters
    ----------
    indices : list of int
        Indices for those atoms that should be constrained.
    mask : list of bool
        One boolean per atom indicating if the atom should be
        constrained or not.

    To fix some atoms for mask list.
    ----------
    like: mask=[atom.symbol == 'Cu' for atom in atoms]  ===>>>   mask = [True, ..., False, ...].
          indices=[atom.index for atom in atomcopy if atom.symbol == 'Cu']  ===>>>  [3, 6, 10, ...].
          mask = [ atom.position[2] -  min(hetero.get_positions()[:,2]) < 5 for atom in hetero].
    """
    atomcopy = atoms.copy()
    c = FixAtoms(indices = indices, mask = mask)
    atomcopy.set_constraint(c)
    return atomcopy