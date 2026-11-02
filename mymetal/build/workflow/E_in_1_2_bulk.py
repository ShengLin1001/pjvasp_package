"""
E_in_1_2_bulk workflow module

Uniaxial / biaxial in-plane stretch of a bulk cell: the generation core of
``pei_vasp_run_E_in_1_2_bulk`` (formerly a bash script with inline Python; same CLI,
same output). Driven from a directory that contains ``y_full_relax`` it builds::

    y_E_in_1_2_bulk/                     (-outdir)
        y_perturb_seed.txt               only with -perturb > 0
        y_dir/
            a1-<a1>-a2-<a2>/             POSCAR + calc files (+ data.lmp for -calc lammps)

Deformation (unchanged from the bash version): a0 is the scale factor (line 2) of
``y_full_relax/CONTCAR`` (POSCAR if there is no CONTCAR). The x components of all cell vectors
are multiplied by a1/a0, the y components by a2/a0, and the z components divided by their
product, so the volume is conserved; fractional coordinates are kept. a1 / a2 are absolute
lengths in Angstrom, not strains. With a hexagonal reference a = (a0, 0, 0),
b = (-a0/2, a0*sqrt(3)/2, 0) this gives a = (a1, 0, 0), b = (-a1/2, a2*sqrt(3)/2, 0).

Only the calc files differ between engines (:data:`DICT_CALC_FILES`): VASP copies
INCAR / KPOINTS / POTCAR / sub.* / Y_* / CHGCAR, LAMMPS copies ``lmp.in`` (see
:mod:`mymetal.build.workflow.lammps_input`) and also writes ``data.lmp``. The POSCAR is written
either way, so ``pei_vasp_plot_E_in_1_2_bulk.py`` post-processes both. It only prepares
directories -- no Slurm script, no sbatch.

Functions:
    - check_calc: Validate the engine name.
    - get_reference: Read the reference cell and its a0.
    - get_a1_a2_pairs: Grid (every a1 x a2) or paired (a1[i], a2[i]) in-plane targets.
    - get_strained_atoms: Apply the in-plane stretch with volume conservation.
    - get_perturbed_atoms: Random uniform displacement of every coordinate.
    - copy_calc_files: Copy one engine's calc files into a case directory.
    - write_case_structure: POSCAR (+ data.lmp) of one case.
    - generate_E_in_1_2_bulk_dirs: Build the whole tree.
"""

import glob
import shutil
from pathlib import Path

import numpy as np
from ase import Atoms

from mymetal.io.lammps import write_lmp_data
from mymetal.io.vasp import my_read_vasp, my_write_vasp
from mymetal.universal.print.print import confirm_prepare_outdir, fail, warn

# engine -> (required calc files, optional glob patterns); the VASP entry is the cp list of
# the old bash version, where CHGCAR only speeds up the SCF and may be absent
DICT_CALC_FILES = {
    'vasp':   (['INCAR', 'KPOINTS', 'POTCAR'], ['sub.*', 'Y_*', 'CHGCAR']),
    'lammps': (['lmp.in'], ['sub.*']),
}

VOLUME_TOL = 1e-6           # A^3, volume conservation is exact up to round-off
SEED_FILE = 'y_perturb_seed.txt'


def check_calc(calc: str = None) -> str:
    """Return ``calc`` if it is a known engine, else fail.

    Args:
        calc: Engine name, a key of :data:`DICT_CALC_FILES`.

    Returns:
        str: ``calc``.
    """
    if calc not in DICT_CALC_FILES:
        fail('calc must be one of %s, got %r' % (sorted(DICT_CALC_FILES), calc))
    return calc


def get_reference(path_src: Path = None) -> tuple:
    """Read the reference cell from ``CONTCAR`` (fallback ``POSCAR``) of ``path_src``.

    Args:
        path_src: The ``y_full_relax`` directory.

    Returns:
        tuple: ``(atoms, a0)``; ``a0`` is the POSCAR scale factor.
    """
    path_contcar = path_src / 'CONTCAR'
    if not path_contcar.is_file():
        path_contcar = path_src / 'POSCAR'
        if not path_contcar.is_file():
            fail('neither CONTCAR nor POSCAR found in %s' % path_src)
        warn('no CONTCAR in %s, using POSCAR as the reference' % path_src)
    atoms, a0 = my_read_vasp(str(path_contcar))
    print('📍 reference: %s  (a0 = %.8f, %d atoms)' % (path_contcar, a0, len(atoms)))
    return atoms, a0


def get_a1_a2_pairs(la1: list = None, la2: list = None, if_pair: bool = False) -> list:
    """In-plane targets: every a1 x a2 combination, or a1[i] with a2[i].

    Args:
        la1: Target x lengths, Angstrom.
        la2: Target y lengths, Angstrom.
        if_pair: Pair element-wise (``-pair``) instead of the full grid.

    Returns:
        list: ``[(a1, a2), ...]`` in generation order.
    """
    if not la1 or not la2:
        fail('both -a1 and -a2 need at least one value')
    if if_pair:
        if len(la1) != len(la2):
            fail('-pair needs -a1 and -a2 of the same length, got %d and %d'
                 % (len(la1), len(la2)))
        return list(zip(la1, la2))
    return [(a1, a2) for a1 in la1 for a2 in la2]


def get_strained_atoms(atoms_ref: Atoms = None, a0: float = None, a1: float = None,
                       a2: float = None) -> tuple:
    """Stretch x by a1/a0 and y by a2/a0, shrink z by their product (volume conserved).

    Args:
        atoms_ref: Reference cell, its a1 / a2 vectors in units of ``a0``.
        a0: Reference scale factor.
        a1: Target x length.
        a2: Target y length.

    Returns:
        tuple: ``(atoms, strainx, strainy)``.
    """
    strainx = a1 / a0
    strainy = a2 / a0
    atoms = atoms_ref.copy()
    cell = np.array(atoms.get_cell())
    cell[:, 0] *= strainx
    cell[:, 1] *= strainy
    cell[:, 2] /= (strainx * strainy)
    atoms.set_cell(cell, scale_atoms=True)
    if abs(atoms.get_volume() - atoms_ref.get_volume()) > VOLUME_TOL:
        warn('volume not conserved: %.8f -> %.8f'
             % (atoms_ref.get_volume(), atoms.get_volume()))
    return atoms, strainx, strainy


def get_perturbed_atoms(atoms: Atoms = None, amp: float = None,
                        rng: np.random.Generator = None) -> Atoms:
    """Displace every Cartesian coordinate by a uniform random number in [-amp, amp].

    Args:
        atoms: Structure to perturb (not modified).
        amp: Amplitude, Angstrom.
        rng: Random generator; one generator walks through all cases in order, so the
            whole set is reproducible from the seed in :data:`SEED_FILE`.

    Returns:
        Atoms: Perturbed copy.
    """
    atoms = atoms.copy()
    atoms.set_positions(atoms.get_positions()
                        + rng.uniform(-amp, amp, size=(len(atoms), 3)))
    return atoms


def copy_calc_files(path_src: Path = None, path_case: Path = None, calc: str = None) -> None:
    """Copy one engine's calc files from ``path_src`` into ``path_case``.

    Args:
        path_src: ``y_full_relax``.
        path_case: Target case directory.
        calc: Engine name.
    """
    lrequired, lglob = DICT_CALC_FILES[calc]
    for name in lrequired:
        shutil.copy(path_src / name, path_case / name)
    for pattern in lglob:
        for name in sorted(glob.glob(str(path_src / pattern))):
            shutil.copy(name, path_case / Path(name).name)


def write_case_structure(path_case: Path = None, atoms: Atoms = None, a0: float = None,
                         label: str = None, calc: str = None) -> None:
    """POSCAR with scale factor a0 (the contract the post-processing reads), plus data.lmp.

    Args:
        path_case: Case directory.
        atoms: Case structure.
        a0: Scale factor written on POSCAR line 2.
        label: POSCAR comment line.
        calc: Engine name.
    """
    my_write_vasp(str(path_case / 'POSCAR'), atoms, lattice_scale_factor=a0, label=label)
    if calc == 'lammps':
        write_lmp_data(path_case / 'data.lmp', atoms)


def generate_E_in_1_2_bulk_dirs(path_root: Path = None, la1: list = None, la2: list = None,
                                if_pair: bool = False, calc: str = 'vasp',
                                perturb: float = 0.0, srcdir: str = 'y_full_relax',
                                outdir: str = 'y_E_in_1_2_bulk',
                                force: bool = False) -> list:
    """Build ``<outdir>/y_dir/a1-<a1>-a2-<a2>/`` for every in-plane target.

    Args:
        path_root: Directory that contains ``srcdir``; must be absolute.
        la1: Target x lengths, Angstrom.
        la2: Target y lengths, Angstrom.
        if_pair: Pair a1[i] with a2[i] instead of the full grid.
        calc: ``'vasp'`` or ``'lammps'``.
        perturb: Random displacement amplitude, Angstrom; 0 disables it.
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
    calc = check_calc(calc)
    if perturb < 0:
        fail('perturb must be >= 0, got %g' % perturb)
    path_src = path_root / srcdir
    if not path_src.is_dir():
        fail('no %s found in %s' % (srcdir, path_root))
    lmissing = [name for name in DICT_CALC_FILES[calc][0] if not (path_src / name).is_file()]
    if lmissing:
        fail('-calc %s needs %s in %s' % (calc, ', '.join(lmissing), path_src))
    lpair = get_a1_a2_pairs(la1, la2, if_pair)
    ### to here

    ### prepare
    atoms_ref, a0 = get_reference(path_src)
    path_out = path_root / outdir
    confirm_prepare_outdir(path_out, force=force)
    (path_out / 'y_dir').mkdir(parents=True)
    rng = None
    if perturb > 0:
        seed = np.random.SeedSequence().entropy      # drawn per run, recorded for replay
        rng = np.random.default_rng(seed)
        (path_out / SEED_FILE).write_text('seed %d\namp %.8f\n' % (seed, perturb))
        print('🎲 perturb ±%.4f Å, seed %d -> %s' % (perturb, seed, path_out / SEED_FILE))
    ### to here

    ### main
    lpath_case = []
    for a1, a2 in lpair:
        name = 'a1-%.8f-a2-%.8f' % (a1, a2)
        path_case = path_out / 'y_dir' / name
        print('================ 📁 %s' % name)
        path_case.mkdir()
        atoms, strainx, strainy = get_strained_atoms(atoms_ref, a0, a1, a2)
        label = 'strainx:%s strainy:%s val1:%.8f  val2:%.8f' % (strainx, strainy, a1, a2)
        if rng is not None:
            atoms = get_perturbed_atoms(atoms, perturb, rng)
            label += ' perturb:%g' % perturb
        copy_calc_files(path_src, path_case, calc)
        write_case_structure(path_case, atoms, a0, label, calc)
        lpath_case.append(path_case)
    ### to here

    print('================ 📊 summary')
    print('✅ %d case dirs (%s, %s) in %s'
          % (len(lpath_case), 'paired' if if_pair else 'grid', calc, path_out / 'y_dir'))
    print('▶️  it does NOT sbatch; submit %s/y_dir/* yourself' % outdir)
    return lpath_case
