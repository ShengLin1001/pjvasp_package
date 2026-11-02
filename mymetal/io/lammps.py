"""
lammps io submodule

LAMMPS-side structure files for workflows shared with VASP. A case directory always keeps a
VASP-format POSCAR / CONTCAR (the structure contract every ``y_dir`` post-processor reads);
these helpers write the LAMMPS input structure next to it and turn LAMMPS output back into a
CONTCAR, so the VASP post-processing runs unchanged on LAMMPS results.

Frame convention: LAMMPS needs a lower-triangular cell (a along x, b in the xy plane). ASE's
:class:`~ase.calculators.lammps.Prism` rotates a general cell into that frame on write; a
CONTCAR rebuilt from a dump is left in the LAMMPS frame. Every ``y_dir`` workflow here starts
from cells that already are lower-triangular, so the rotation is the identity in practice.

Functions:
    - write_lmp_data: Write ``data.lmp`` (atomic style) for an Atoms object.
    - write_lmp_neb_coords: Write one replica coordinate file for ``neb ... each``.
    - read_lmp_dump_atoms: Rebuild an Atoms object from a final LAMMPS dump.
    - read_lmp_dump_forces: Forces of the last frame of a dump, if it has fx fy fz.
"""

from pathlib import Path

import numpy as np
from ase import Atoms
from ase.calculators.lammps import Prism
from ase.io import read, write


def write_lmp_data(path_data: Path = None, atoms: Atoms = None) -> None:
    """Write a LAMMPS data file in ``atom_style atomic``.

    Atom types follow ASE's default species order (alphabetical), the same as the existing
    ``lmp_utils`` model scripts, so ``mass <type> ...`` lines stay interchangeable.

    Args:
        path_data: Output path, e.g. ``<case>/data.lmp``.
        atoms: Structure to write.
    """
    write(Path(path_data), atoms, format='lammps-data', atom_style='atomic')


def write_lmp_neb_coords(path_coords: Path = None, atoms: Atoms = None) -> None:
    """Write a replica coordinate file read by ``neb ... each <file>``.

    Format (LAMMPS ``neb`` doc): first line the atom count, then ``id x y z`` per atom, in the
    LAMMPS frame. Positions go through the same Prism as :func:`write_lmp_data` so a replica
    lines up with the data file of replica 0.

    Args:
        path_coords: Output path, e.g. ``y_neb/y_dir/coords.3``.
        atoms: Image structure; must share cell and atom order with replica 0.
    """
    prism = Prism(np.array(atoms.get_cell()))
    lpos = prism.vector_to_lammps(atoms.get_positions(), wrap=True)
    llines = ['%d' % len(atoms)]
    for i, pos in enumerate(lpos, start=1):
        llines.append('%d %.12f %.12f %.12f' % (i, pos[0], pos[1], pos[2]))
    Path(path_coords).write_text('\n'.join(llines) + '\n')


def read_lmp_dump_atoms(path_dump: Path = None, atoms_ref: Atoms = None) -> Atoms:
    """Rebuild the final structure from a LAMMPS ``custom`` dump (``id type x y z``).

    The dump carries no element names, so symbols and atom order come from ``atoms_ref``
    (the case POSCAR the data file was written from); cell and positions come from the
    last dump frame (ASE sorts the atoms by id).

    Args:
        path_dump: Dump path, e.g. ``<case>/y_final.lammpstrj``.
        atoms_ref: Input structure of the same case.

    Returns:
        Atoms: Final structure in the LAMMPS frame, periodic in all directions.
    """
    atoms_dump = read(Path(path_dump), format='lammps-dump-text', index=-1)
    if len(atoms_dump) != len(atoms_ref):
        raise ValueError('dump %s has %d atoms, reference has %d'
                         % (path_dump, len(atoms_dump), len(atoms_ref)))
    # A box/relax that dilates about the box centre moves the origin off 0; ASE keeps it as
    # celldisp. Subtract it, or every fractional coordinate carries a rigid shift (harmless
    # for the energy, but a fake translation between NEB endpoints).
    lpos = atoms_dump.get_positions() - np.array(atoms_dump.get_celldisp()).reshape(3)
    atoms = Atoms(symbols=atoms_ref.get_chemical_symbols(), positions=lpos,
                  cell=atoms_dump.get_cell(), pbc=True)
    return atoms


def read_lmp_dump_forces(path_dump: Path = None) -> np.ndarray:
    """Forces (eV/A, sorted by atom id) of the last dump frame, or None without fx fy fz.

    Args:
        path_dump: Dump path, e.g. ``y_neb/y_dir/y_final.3.lammpstrj``.

    Returns:
        np.ndarray | None: (natom, 3) forces.
    """
    atoms_dump = read(Path(path_dump), format='lammps-dump-text', index=-1)
    if atoms_dump.calc is None:
        return None
    return atoms_dump.get_forces()
