"""
lmp_univ post submodule

LAMMPS counterpart of ``pei_vasp_univ_post`` and ``pei_vasp_univ_neb_nebbarrier.py``: turn
finished LAMMPS cases into the files the VASP post-processing already reads, so a
``y_full_relax``-style workflow run with ``-calc lammps`` reuses its plotting unchanged.

Per relaxation case (inputs from :mod:`mymetal.build.workflow.lammps_input`)::

    y_final.lammpstrj -> CONTCAR        (scale factor copied from the case POSCAR)
    y_lmp_energy.txt  -> one line of <workdir>/y_post_data.txt, same columns as VASP:
                         jobn  E  EENTRO  in kB  xx yy zz xy yz zx

Per NEB (``y_neb/y_dir`` from ``pei_vasp_run_neb -calc lammps``)::

    y_final.<i>.lammpstrj -> <ii>/CONTCAR
    log.lammps            -> p_neb_py.dat (id dist_cum energy_rel force_b image), the
                             nebbarrier table; force_b is nan (LAMMPS does not report the
                             tangential force per replica)

Functions:
    - read_lmp_case_energy: The energy / stress line of one case.
    - write_lmp_case_contcar: Final dump -> CONTCAR of one case.
    - post_lmp_cases: y_post_data.txt + CONTCAR for a list of cases.
    - read_lmp_neb_log: Last line of the universe NEB log.
    - get_lmp_neb_ftol: The force tolerance of the ``neb`` command in lmp.in.
    - post_lmp_neb: p_neb_py.dat + image CONTCARs of one NEB.
"""

import math
import re
from pathlib import Path

from ase import Atoms

from mymetal.io.lammps import read_lmp_dump_atoms
from mymetal.io.vasp import my_read_vasp, my_write_vasp
from mymetal.universal.print.print import fail, warn

ENERGY_FILE = 'y_lmp_energy.txt'
DUMP_FILE = 'y_final.lammpstrj'
DONE_MARKER = 'y_lmp_done.txt'
POST_DATA_HEADER = 'i1   energy(sigma->0)(eV)   EENTRO(eV)   -stress(kB)'
NEB_HEADER_KEY = 'Step MaxReplicaForce'
NEB_NCOL_FIXED = 9          # Step MaxReplicaForce MaxAtomForce GradV0 GradV1 GradVc EBF EBR RDT


def read_lmp_case_energy(path_case: Path = None) -> str:
    """Return ``"E EENTRO in kB xx yy zz xy yz zx"`` of a finished case, or None.

    Args:
        path_case: Case directory.

    Returns:
        str | None: The line, or None when the case is unfinished.
    """
    path_energy = path_case / ENERGY_FILE
    if not (path_case / DONE_MARKER).is_file() or not path_energy.is_file():
        return None
    line = path_energy.read_text().strip()
    if len(line.split()) != 10:
        warn('malformed %s: %r' % (path_energy, line))
        return None
    return line


def write_lmp_case_contcar(path_case: Path = None, name_dump: str = DUMP_FILE,
                           name_ref: str = 'POSCAR') -> Atoms:
    """Write ``CONTCAR`` from the final dump, with the scale factor of the case POSCAR.

    Args:
        path_case: Case directory.
        name_dump: Final dump file name.
        name_ref: Input structure the case started from.

    Returns:
        Atoms: The final structure.
    """
    atoms_ref, a0 = my_read_vasp(str(path_case / name_ref))
    atoms = read_lmp_dump_atoms(path_case / name_dump, atoms_ref)
    my_write_vasp(str(path_case / 'CONTCAR'), atoms, lattice_scale_factor=a0,
                  label='LAMMPS final structure of %s' % path_case.name)
    return atoms


def post_lmp_cases(path_work: Path = None, lpath_case: list = None) -> tuple:
    """Collect finished cases into ``<path_work>/y_post_data.txt`` and write their CONTCAR.

    Unfinished cases are left out of the table, as ``pei_vasp_univ_post`` leaves out
    directories without an OUTCAR.

    Args:
        path_work: Directory receiving y_post_data.txt (the parent of ``y_dir``).
        lpath_case: Case directories, in table order.

    Returns:
        tuple: ``(ldone, lunfinished)`` case names.
    """
    ldone = []
    lunfinished = []
    llines = [POST_DATA_HEADER]
    for path_case in lpath_case:
        line = read_lmp_case_energy(path_case)
        if line is None:
            lunfinished.append(path_case.name)
            continue
        write_lmp_case_contcar(path_case)
        lfield = line.split()
        llines.append(path_case.name + '\t' + lfield[0] + '\t' + lfield[1] + '\t'
                      + ' '.join(lfield[2:]))
        ldone.append(path_case.name)
    (Path(path_work) / 'y_post_data.txt').write_text('\n'.join(llines) + '\n')
    return ldone, lunfinished


def read_lmp_neb_log(path_log: Path = None) -> dict:
    """Parse the last data line of a universe NEB log (``log.lammps``).

    Columns: ``Step MaxReplicaForce MaxAtomForce GradV0 GradV1 GradVc EBF EBR RDT`` then
    ``RD<i> PE<i>`` per replica (LAMMPS ``neb`` doc). The last header wins, so the climbing
    stage is read when there is one.

    Args:
        path_log: Universe log path.

    Returns:
        dict: ``step``, ``max_replica_force``, ``max_atom_force``, ``rdt``, ``lrd``, ``lpe``.
    """
    llines = Path(path_log).read_text().splitlines()
    lidx = [i for i, line in enumerate(llines) if line.split()[:2] == NEB_HEADER_KEY.split()]
    if not lidx:
        fail('no NEB table (%r) in %s' % (NEB_HEADER_KEY, path_log))
    ncol = len(llines[lidx[-1]].split())
    lrow = []
    for line in llines[lidx[-1] + 1:]:
        lfield = line.split()
        if len(lfield) != ncol:
            break
        try:
            lrow = [float(x) for x in lfield]
        except ValueError:
            break
    if not lrow:
        fail('NEB table in %s has no data line' % path_log)
    nrep = (ncol - NEB_NCOL_FIXED) // 2
    return {
        'step': int(lrow[0]),
        'max_replica_force': lrow[1],
        'max_atom_force': lrow[2],
        'rdt': lrow[8],
        'lrd': [lrow[NEB_NCOL_FIXED + 2 * i] for i in range(nrep)],
        'lpe': [lrow[NEB_NCOL_FIXED + 2 * i + 1] for i in range(nrep)],
    }


def get_lmp_neb_ftol(path_in: Path = None) -> float:
    """Force tolerance (2nd argument) of the ``neb`` command in ``lmp.in``, or None.

    Args:
        path_in: The NEB input.

    Returns:
        float | None: ftol.
    """
    for line in Path(path_in).read_text().splitlines():
        match = re.match(r'\s*neb\s+(\S+)\s+(\S+)', line)
        if match:
            return float(match.group(2))
    return None


def post_lmp_neb(path_neb: Path = None) -> dict:
    """Write ``p_neb_py.dat`` and the image CONTCARs of a finished LAMMPS NEB.

    Args:
        path_neb: ``y_neb/y_dir`` (image dirs 00 .. N+1, log.lammps, y_final.<i>.lammpstrj).

    Returns:
        dict: The parsed log plus ``converged``, ``ftol``, ``barrier_eV`` (max PE - PE_0).
    """
    path_neb = Path(path_neb)
    if not (path_neb / DONE_MARKER).is_file():
        fail('NEB not finished (no %s) in %s' % (DONE_MARKER, path_neb))
    dict_neb = read_lmp_neb_log(path_neb / 'log.lammps')
    nrep = len(dict_neb['lpe'])
    for i in range(nrep):
        write_lmp_case_contcar(path_neb / ('%02d' % i), name_dump='../y_final.%d.lammpstrj' % i)
    pe0 = dict_neb['lpe'][0]
    llines = ['%3s %16s %16s %16s %8s' % ('id', 'dist_cum(Å)', 'energy_rel(eV)',
                                          'force_b(eV/Å)', 'image')]
    for i, (rd, pe) in enumerate(zip(dict_neb['lrd'], dict_neb['lpe'])):
        llines.append('%3d %16.8f %16.8f %16s %8d'
                      % (i, rd * dict_neb['rdt'], pe - pe0, 'nan', i))
    (path_neb / 'p_neb_py.dat').write_text('\n'.join(llines) + '\n')
    ftol = get_lmp_neb_ftol(path_neb / 'lmp.in')
    dict_neb['ftol'] = ftol
    dict_neb['converged'] = (ftol is not None and dict_neb['max_replica_force'] <= ftol)
    dict_neb['barrier_eV'] = max(dict_neb['lpe']) - pe0
    if not math.isfinite(dict_neb['barrier_eV']):
        warn('non-finite NEB barrier in %s' % path_neb)
    return dict_neb
