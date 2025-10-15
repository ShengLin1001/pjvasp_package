"""
neighbor submodule

This submodule provides functions for calculating neighbor distances in atomic structures. It includes a function to compute
all neighbor distances within a specified cutoff distance, with options for self-interaction and periodic boundary conditions
handling.

Functions:
    - get_neighbor_distances: Calculate all neighbor distances within a specified cutoff distance.
"""
from ase import Atoms
from ase.neighborlist import NeighborList
import numpy as np


def get_neighbor_distances(atoms: Atoms = None, cutoff: float = 10,
                           self_interaction: bool = False, bothways: bool = True, skin: float = 0.0) -> tuple:
    cutoffs = [cutoff / 2] * len(atoms)  # NeighborList needs each atom's radius
    # skin = 0.0 means no extra distance for neighbor searching
    # It's important!
    nl = NeighborList(cutoffs, self_interaction=self_interaction, bothways=bothways, skin = skin)
    nl.update(atoms)

    distances = []

    for i in range(len(atoms)):
        indices, offsets = nl.get_neighbors(i)   # offset - shift the periodic image
        for j, offset in zip(indices, offsets):
            rij = atoms.positions[j] + np.dot(offset, atoms.get_cell()) - atoms.positions[i]
            distances.append(np.linalg.norm(rij))

    return np.array(distances), np.unique(np.array(distances))