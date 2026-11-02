"""
neb workflow module

Prepare a NEB from relaxed endpoints: the generation core of ``pei_vasp_run_neb``
(formerly a bash script; same CLI and output). Driven from a directory that contains::

    y_full_relax_neb/
        ini/CONTCAR  fin/CONTCAR         relaxed endpoints (VASP, or LAMMPS via pei_lmp_univ_post)
        ini/<calc templates>.neb         INCAR.neb KPOINTS.neb POTCAR sub.*.neb Y_*.neb (VASP)
                                         lmp.in.neb sub.*.neb                          (LAMMPS)

it builds ``y_neb/y_dir/00 .. N+1`` with VTST ``nebmake.pl`` (linear interpolation) and
``nebavoid.pl 1`` (no pair closer than 1 A), copies the templates with ``.neb`` stripped, and
for ``-calc lammps`` also writes ``data.lmp`` (image 00) and ``coords.<i>`` (image i >= 1) for
``neb ... each coords.${i}`` -- so both engines start from the very same images. It only
prepares the directory -- no sbatch.

Functions:
    - check_endpoints: Both endpoint CONTCARs exist.
    - run_tool: Run an external helper (dist.pl, nebmake.pl ...) and echo its output.
    - get_image_dirs: The numbered image directories, in order.
    - copy_neb_templates: Copy ``*.neb`` templates (and POTCAR) and strip ``.neb``.
    - write_lmp_neb_structures: data.lmp + coords.<i> from the image POSCARs.
    - print_neb_inputs: Echo the inputs to check before submitting.
    - generate_neb_dirs: Build the whole ``y_neb`` tree.
"""

import glob
import shutil
import subprocess
from pathlib import Path

from mymetal.io.lammps import write_lmp_data, write_lmp_neb_coords
from mymetal.io.vasp import my_read_vasp
from mymetal.universal.print.print import confirm_prepare_outdir, fail, warn

# engine -> template patterns looked up in ini/ (the cp list of the old bash version)
DICT_NEB_TEMPLATES = {
    'vasp':   ['INCAR.neb', 'KPOINTS.neb', 'POTCAR', 'sub.*.neb', 'Y_*.neb'],
    'lammps': ['lmp.in.neb', 'sub.*.neb'],
}
# engine -> files echoed for a last look before sbatch
DICT_NEB_ECHO = {
    'vasp':   ['INCAR', 'KPOINTS', 'sub.*', 'Y_*'],
    'lammps': ['lmp.in', 'sub.*'],
}
NEBAVOID_MIN_DIST = '1'         # Angstrom, as in the bash version


def check_endpoints(path_src: Path = None) -> tuple:
    """Return the ini / fin CONTCAR paths, failing if one is missing.

    Args:
        path_src: ``y_full_relax_neb``.

    Returns:
        tuple: ``(path_ini, path_fin)``.
    """
    lpath = [path_src / 'ini' / 'CONTCAR', path_src / 'fin' / 'CONTCAR']
    for path in lpath:
        if path.is_file():
            print('✅ %s exists.' % path)
        else:
            fail('%s not found. Not fully relaxed.' % path)
    return lpath[0], lpath[1]


def run_tool(lcmd: list = None, path_cwd: Path = None) -> str:
    """Run an external helper in ``path_cwd``, echo stdout / stderr, fail on non-zero exit.

    Args:
        lcmd: Command and arguments.
        path_cwd: Working directory.

    Returns:
        str: Captured stdout.
    """
    print('▶️  %s' % ' '.join(str(x) for x in lcmd))
    result = subprocess.run([str(x) for x in lcmd], cwd=str(path_cwd),
                            capture_output=True, text=True)
    if result.stdout:
        print(result.stdout.rstrip())
    if result.stderr:
        print(result.stderr.rstrip())
    if result.returncode != 0:
        fail('%s exited with %d' % (lcmd[0], result.returncode))
    return result.stdout


def get_image_dirs(path_ydir: Path = None, n_image: int = None) -> list:
    """The image directories 00 .. n_image+1 that nebmake.pl created.

    Args:
        path_ydir: ``y_neb/y_dir``.
        n_image: Number of intermediate images.

    Returns:
        list: Image directory paths, in order.
    """
    lpath_image = [path_ydir / ('%02d' % i) for i in range(n_image + 2)]
    lmissing = [p.name for p in lpath_image if not (p / 'POSCAR').is_file()]
    if lmissing:
        fail('nebmake.pl did not produce POSCAR in %s' % ', '.join(lmissing))
    return lpath_image


def copy_neb_templates(path_ini: Path = None, path_ydir: Path = None, calc: str = None) -> list:
    """Copy the engine's templates from ``ini/`` and strip the ``.neb`` suffix.

    Args:
        path_ini: ``y_full_relax_neb/ini``.
        path_ydir: ``y_neb/y_dir``.
        calc: Engine name.

    Returns:
        list: Names written into ``path_ydir``.
    """
    lname = []
    for pattern in DICT_NEB_TEMPLATES[calc]:
        lsrc = sorted(glob.glob(str(path_ini / pattern)))
        if not lsrc:
            warn('no %s in %s' % (pattern, path_ini))
        for src in lsrc:
            name = Path(src).name
            name_new = name[:-len('.neb')] if name.endswith('.neb') else name
            shutil.copy(src, path_ydir / name_new)
            print('🔁 %s => %s' % (name, name_new))
            lname.append(name_new)
    return lname


def write_lmp_neb_structures(path_ydir: Path = None, lpath_image: list = None) -> None:
    """data.lmp from image 00 and coords.<i> from image i (i >= 1).

    Args:
        path_ydir: ``y_neb/y_dir``.
        lpath_image: Image directories, in order.
    """
    atoms0, _ = my_read_vasp(str(lpath_image[0] / 'POSCAR'))
    write_lmp_data(path_ydir / 'data.lmp', atoms0)
    print('📄 data.lmp <- %s/POSCAR' % lpath_image[0].name)
    for i, path_image in enumerate(lpath_image[1:], start=1):
        atoms, _ = my_read_vasp(str(path_image / 'POSCAR'))
        if len(atoms) != len(atoms0):
            fail('image %s has %d atoms, image 00 has %d'
                 % (path_image.name, len(atoms), len(atoms0)))
        write_lmp_neb_coords(path_ydir / ('coords.%d' % i), atoms)
    print('📄 coords.1 .. coords.%d <- image POSCARs' % (len(lpath_image) - 1))


def print_neb_inputs(path_ydir: Path = None, lpath_image: list = None, calc: str = None) -> None:
    """Echo the image headers and the copied inputs, as the bash version did.

    Args:
        path_ydir: ``y_neb/y_dir``.
        lpath_image: Image directories.
        calc: Engine name.
    """
    for path_image in lpath_image:
        print('-' * 59)
        print('📄 POSCAR in image %s' % path_image.name)
        print('\n'.join((path_image / 'POSCAR').read_text().splitlines()[:10]))
    for pattern in DICT_NEB_ECHO[calc]:
        for name in sorted(glob.glob(str(path_ydir / pattern))):
            print('================ 🧾 %s' % Path(name).name)
            for k, line in enumerate(Path(name).read_text().splitlines(), start=1):
                print('%6d  %s' % (k, line))
    if calc == 'vasp' and (path_ydir / 'POTCAR').is_file():
        print('================ 🧪 POTCAR head (should contain COPYR)')
        print('\n'.join((path_ydir / 'POTCAR').read_text().splitlines()[:10]))


def generate_neb_dirs(path_root: Path = None, n_image: int = None, calc: str = 'vasp',
                      srcdir: str = 'y_full_relax_neb', outdir: str = 'y_neb',
                      force: bool = False) -> Path:
    """Build ``<outdir>/y_dir`` with ``n_image`` interpolated images between ini and fin.

    Args:
        path_root: Directory containing ``srcdir``; must be absolute.
        n_image: Number of intermediate images (>= 0).
        calc: ``'vasp'`` or ``'lammps'``.
        srcdir: Endpoint directory name.
        outdir: Output directory name.
        force: Delete an existing ``outdir`` without prompting.

    Returns:
        Path: ``<outdir>/y_dir``.
    """
    ### check
    path_root = Path(path_root)
    if not path_root.is_absolute():
        fail('path_root must be absolute, got %s' % path_root)
    if calc not in DICT_NEB_TEMPLATES:
        fail('calc must be one of %s, got %r' % (sorted(DICT_NEB_TEMPLATES), calc))
    if not isinstance(n_image, int) or n_image < 0:
        fail('N_IMAGE must be a non-negative integer, got %r' % n_image)
    path_src = path_root / srcdir
    if not path_src.is_dir():
        fail('No %s found in %s! Abort!' % (srcdir, path_root))
    path_ini, path_fin = check_endpoints(path_src)
    if calc == 'lammps' and not (path_src / 'ini' / 'lmp.in.neb').is_file():
        fail('-calc lammps needs ini/lmp.in.neb in %s' % path_src)
    ### to here

    ### prepare
    path_out = path_root / outdir
    print('🧹 Preparing clean NEB workspace at: %s' % path_out)
    confirm_prepare_outdir(path_out, force=force)
    path_ydir = path_out / 'y_dir'
    path_ydir.mkdir(parents=True)
    ### to here

    ### main
    print('📏 Running dist.pl to check distance...')
    run_tool(['dist.pl', path_ini, path_fin], path_ydir)
    print('🧬 Running nebmake.pl to generate %d images...' % n_image)
    run_tool(['nebmake.pl', path_ini, path_fin, n_image], path_ydir)
    print('🛡️  Running nebavoid.pl to keep every interatomic distance >= %s A...'
          % NEBAVOID_MIN_DIST)
    run_tool(['nebavoid.pl', NEBAVOID_MIN_DIST], path_ydir)
    lpath_image = get_image_dirs(path_ydir, n_image)

    # endpoint results go into 00 / N+1 so the image chain carries its own ends
    for path_end_src, path_image in [(path_src / 'ini', lpath_image[0]),
                                     (path_src / 'fin', lpath_image[-1])]:
        if calc == 'vasp':
            if (path_end_src / 'OUTCAR').is_file():
                shutil.copy(path_end_src / 'OUTCAR', path_image / 'OUTCAR')
            else:
                warn('no OUTCAR in %s' % path_end_src)
        shutil.copy(path_image / 'POSCAR', path_image / 'CONTCAR')

    print('🎞️  Generating movie with pei_vasp_univ_neb_nebmovie...')
    run_tool(['pei_vasp_univ_neb_nebmovie'], path_ydir)
    for name in sorted(glob.glob(str(path_ydir / 'movie*.xyz'))):
        shutil.move(name, path_out / Path(name).name)

    print('📥 Copying input templates from ini/...')
    copy_neb_templates(path_src / 'ini', path_ydir, calc)
    if calc == 'lammps':
        write_lmp_neb_structures(path_ydir, lpath_image)
    print_neb_inputs(path_ydir, lpath_image, calc)
    ### to here

    print('================ 📊 summary')
    print('✅ %d images (00 .. %02d, %s) in %s' % (len(lpath_image), n_image + 1, calc, path_ydir))
    if calc == 'lammps':
        print('▶️  run in %s: lmp -partition %dx1 -in lmp.in (e.g. via pei_lmp_univ_sbatch)'
              % (path_ydir, n_image + 2))
    print('🚀 Ready to submit job in directory: %s' % path_ydir)
    return path_ydir
