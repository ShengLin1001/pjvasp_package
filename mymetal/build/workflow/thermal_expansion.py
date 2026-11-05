"""
thermal_expansion workflow module

Lattice parameters vs temperature from NPT MD: the generation core of
``pei_vasp_run_thermal_expansion`` (LAMMPS only). Driven from a directory that contains
``y_full_relax`` it builds::

    y_thermal_expansion/                 (-outdir)
        y_dir/
            T-<T>/                       POSCAR (supercell) + data.lmp + lmp.in

``y_full_relax`` holds the unit cell (CONTCAR, else POSCAR) and ``lmp.in.npt``, the NPT input
rendered by :func:`mymetal.build.workflow.lammps_input.generate_lmp_npt_input`. Like INCAR /
POTCAR for VASP, that file carries the potential and the MD settings; the workflow itself
knows nothing about either. Each case is an independent run from the same supercell: its
``lmp.in`` is ``lmp.in.npt`` with ``variable temp`` / ``variable seed`` put on top, and the
POSCAR comment line records ``T:<T> size:<nx>x<ny>x<nz> seed:<seed>`` for the post-processing.
The seed is drawn per case and never fixed. It only prepares directories -- no Slurm
script, no sbatch.

Cells must be orthogonal (FCC cubic, HCP orthorhombic a x sqrt(3)a x c): the input couples
x and y, which keeps such a hexagonal cell hexagonal, while a 120-degree cell sits exactly on
the LAMMPS tilt limit |xy| = lx/2 and flips under NPT.

Functions:
    - check_temperatures: Positive, unique temperatures in K.
    - check_size: Three positive supercell multipliers.
    - check_orthogonal: Fail for a tilted reference cell.
    - get_case_input: ``lmp.in.npt`` with this case's temperature and seed on top.
    - generate_thermal_expansion_dirs: Build the whole tree.
"""

import glob
import shutil
from pathlib import Path

import numpy as np

from mymetal.build.workflow.E_in_1_2_bulk import get_reference, write_case_structure
from mymetal.universal.print.print import confirm_prepare_outdir, fail

NPT_INPUT = 'lmp.in.npt'
CASE_INPUT = 'lmp.in'
CASE_FORMAT = 'T-%04d'
LGLOB_EXTRA = ['sub.*']          # copied along, as for the other -calc lammps workflows
SEED_MAX = 900000000             # velocity create takes a positive 32-bit integer
TILT_TOL = 1e-8                  # A, off-diagonal cell components treated as zero


def check_temperatures(ltemp: list = None) -> list:
    """Return ``ltemp`` sorted if all are positive and unique, else fail.

    Args:
        ltemp: Temperatures, K (whole numbers; they name the case dirs).

    Returns:
        list: Sorted temperatures.
    """
    if not ltemp:
        fail('give at least one temperature')
    lbad = [temp for temp in ltemp if temp <= 0 or int(temp) != temp]
    if lbad:
        fail('temperatures must be whole K > 0 (velocity create, dir names), got %s' % lbad)
    if len(set(ltemp)) != len(ltemp):
        fail('duplicate temperatures in %s' % ltemp)
    return sorted(int(temp) for temp in ltemp)


def check_size(lsize: list = None) -> list:
    """Return ``lsize`` if it is three positive integers, else fail.

    Args:
        lsize: Supercell multipliers ``[nx, ny, nz]`` of the unit cell.

    Returns:
        list: ``lsize``.
    """
    if lsize is None or len(lsize) != 3 or any(int(n) != n or n <= 0 for n in lsize):
        fail('size must be three positive integers, got %s' % lsize)
    return [int(n) for n in lsize]


def check_orthogonal(atoms=None) -> None:
    """Fail when the reference cell has off-diagonal components (see module docstring).

    Args:
        atoms: Reference unit cell.
    """
    cell = np.array(atoms.get_cell())
    if np.abs(cell - np.diag(np.diag(cell))).max() > TILT_TOL:
        fail('the reference cell must be orthogonal (FCC cubic, HCP orthorhombic), got\n%s'
             % cell)


def get_case_input(text_npt: str = None, temp: int = None, seed: int = None) -> str:
    """``lmp.in.npt`` with the variables it leaves open defined on top.

    Args:
        text_npt: Content of ``y_full_relax/lmp.in.npt``.
        temp: Temperature, K.
        seed: Velocity seed.

    Returns:
        str: The case ``lmp.in``.
    """
    return ('# written by pei_vasp_run_thermal_expansion: this case\'s temperature (K) and '
            'velocity seed\n'
            'variable        temp equal ' + str(temp) + '\n'
            'variable        seed equal ' + str(seed) + '\n\n' + text_npt)


def generate_thermal_expansion_dirs(path_root: Path = None, ltemp: list = None,
                                    lsize: list = None, srcdir: str = 'y_full_relax',
                                    outdir: str = 'y_thermal_expansion',
                                    force: bool = False) -> list:
    """Build ``<outdir>/y_dir/T-<T>/`` for every temperature.

    Args:
        path_root: Directory that contains ``srcdir``; must be absolute.
        ltemp: Temperatures, K.
        lsize: Supercell multipliers ``[nx, ny, nz]``.
        srcdir: Reference directory name.
        outdir: Output directory name.
        force: Delete an existing ``outdir`` without prompting.

    Returns:
        list: Created case directories.
    """
    ### check
    path_root = Path(path_root)
    if not path_root.is_absolute():
        fail('path_root must be absolute, got %s' % path_root)
    ltemp = check_temperatures(ltemp)
    lsize = check_size(lsize)
    path_src = path_root / srcdir
    if not path_src.is_dir():
        fail('no %s found in %s' % (srcdir, path_root))
    if not (path_src / NPT_INPUT).is_file():
        fail('no %s in %s (render it with generate_lmp_npt_input)' % (NPT_INPUT, path_src))
    ### to here

    ### prepare
    atoms_ref, a0 = get_reference(path_src)
    check_orthogonal(atoms_ref)
    atoms_super = atoms_ref.repeat(lsize)
    text_npt = (path_src / NPT_INPUT).read_text()
    path_out = path_root / outdir
    confirm_prepare_outdir(path_out, force=force)
    (path_out / 'y_dir').mkdir(parents=True)
    rng = np.random.default_rng()
    print('📍 supercell %s -> %d atoms' % ('x'.join(str(n) for n in lsize), len(atoms_super)))
    ### to here

    ### main
    lpath_case = []
    for temp in ltemp:
        name = CASE_FORMAT % temp
        path_case = path_out / 'y_dir' / name
        print('================ 📁 %s' % name)
        path_case.mkdir()
        seed = int(rng.integers(1, SEED_MAX))
        label = 'T:%d size:%s seed:%d' % (temp, 'x'.join(str(n) for n in lsize), seed)
        (path_case / CASE_INPUT).write_text(get_case_input(text_npt, temp, seed))
        for pattern in LGLOB_EXTRA:
            for name_extra in sorted(glob.glob(str(path_src / pattern))):
                shutil.copy(name_extra, path_case / Path(name_extra).name)
        write_case_structure(path_case, atoms_super, a0, label, 'lammps')
        lpath_case.append(path_case)
    ### to here

    print('================ 📊 summary')
    print('✅ %d case dirs (%d atoms each) in %s'
          % (len(lpath_case), len(atoms_super), path_out / 'y_dir'))
    print('▶️  it does NOT sbatch; submit %s/y_dir/* yourself' % outdir)
    return lpath_case
