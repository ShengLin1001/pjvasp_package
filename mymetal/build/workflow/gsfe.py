"""
gsfe workflow module

Extend a generalized-stacking-fault-energy (GSFE) gamma-line produced by the bash
workflow ``pei_vasp_run_gsfe`` from the short partial-slip segment it lays out by
default (``N1=10`` images spanning ``(bp1*a11, bp2*b22)``) to ONE FULL SLIP PERIOD --
the fault vector at which the sheared cell is the perfect crystal again.

Driven from the directory that holds ``y_full_relax`` and ``y_gsfe_<mode>``, it appends::

    y_gsfe_<mode>/
        y_dir/
            00 .. 10           # already there, never touched
            11 .. <period-1>   # created here; same inputs, only POSCAR line 5 differs
        y_sbatch_extend.sh     # auto-generated submit list, NEW dirs only

Two design points that matter:

1. The shear step is taken from the EXISTING directories (a least-squares fit through
   the in-plane part of ``a3`` over index), not recomputed from ``bp1``/``bp2``. New
   images then land exactly on the line the old ones sit on --
   ``mymetal.post.gsfe.check_constraints`` aborts on a 1e-10 collinearity deviation,
   and ``bc``-vs-numpy rounding is not worth risking there.

2. The period, and the symmetry-equivalent images inside it, are detected from the
   structures themselves. A mode whose second half is only the mirror of its first
   half is then reported as ALREADY COMPLETE instead of being recomputed: for HCP
   prism1w image ``i`` and image ``20-i`` are the same structure up to symmetry, so
   the missing half carries no information a NNP could learn from, and FCC_100 walks
   a whole period in its default 10 steps to begin with. Only genuinely inequivalent
   images are laid out; ``include_equivalent=True`` overrides that.

Functions:
    - get_gsfe_line: read ``y_dir/NN``, return reference atoms, shear step, existing indices.
    - get_structure_fingerprint: symmetry-invariant fingerprint of one sheared image.
    - shear_reference_cell: build the image at a given index by tilting ``a3`` in plane.
    - get_gsfe_period: full-period index plus the groups of equivalent images.
    - generate_gsfe_extend_dirs: create the missing ``y_dir/NN`` directories.

Change log:
    - Written by J. P. on 2026-08-20, to turn the partial gamma-lines of the Au FCC/HCP
      dataset (``construct_dataset/calculate/A21-2``, ``A22-2``) into closed full-slip
      curves without recomputing anything already converged.
"""

import glob
import shutil
import numpy as np
from pathlib import Path

from ase import Atoms
from ase.io import read
from ase.neighborlist import neighbor_list

from mymetal.slurm.submit import stamp_comment_tag
from mymetal.universal.print.print import fail, warn

# Slip modes of pei_vasp_run_gsfe, as (bp1, bp2) in units of a11 / b22. Only used to
# cross-check the step reconstructed from the existing directories -- a mismatch means
# the y_dir on disk was not made by the mode the caller named.
DICT_SLIP_MODE = {
    "FCC_111":     (-0.5, 1.0 / 3.0),
    "FCC_100":     (0.5, 0.5),
    "HCP_basal":   (-0.5, 1.0 / 3.0),
    "HCP_prism1w": (0.5, 0.0),
    "HCP_pyr1n":   (0.5, -0.11),
    "HCP_pyr1w":   (0.0, 0.5),
    "HCP_pyr2":    (0.0, 0.5),
    "a1":          (1.0, 0.0),
}

# Files pei_vasp_run_gsfe seeds every image with. Kept verbatim so an extended image is
# byte-identical to an old one apart from POSCAR line 5; CHGCAR is the perfect-crystal
# charge density and only serves as a starting guess (INCAR has ICHARG=1).
LFILE_SEED = ["INCAR", "KPOINTS", "POSCAR", "POTCAR", "Y_CONSTR_LATT", "CHGCAR"]
LGLOB_SEED = ["sub.*"]

# printf format of the bash workflow's line 5, so old and new POSCARs look the same.
FORMAT_A3 = "    %.16f    %.16f    %.16f"

# Fingerprint = sorted list of periodic pair distances below the cutoff. Symmetry-
# equivalent images give the same multiset; the converse can fail for homometric pairs,
# which is why the equivalence groups are always printed for the caller to eyeball.
FINGERPRINT_RCUT = 6.0
FINGERPRINT_TOL = 2.0e-4

# How far past the existing images to look for the closing index before giving up. A
# mode with an irrational (bp1, bp2) -- pyr1n's -0.11 -- never closes, and should fail
# loudly rather than probe forever.
INDEX_PROBE_MAX = 120

# Deviation of the existing images from a straight, evenly spaced line, in Angstrom.
# mymetal.post.gsfe rejects the whole run above 1e-10, so anything looser is a bug here.
TOL_LINE = 1.0e-10


def get_structure_fingerprint(atoms: Atoms = None, rcut: float = FINGERPRINT_RCUT) -> np.ndarray:
    """Sorted periodic pair distances, rounded, as a symmetry-invariant fingerprint.

    Args:
        atoms (Atoms): One sheared image.
        rcut (float): Neighbour cutoff in Angstrom.

    Returns:
        np.ndarray: Sorted distances rounded to 4 decimals.
    """
    return np.round(np.sort(neighbor_list("d", atoms, rcut)), 4)


def shear_reference_cell(atoms_ref: Atoms = None, step_xy: np.ndarray = None,
                         index: int = 0) -> Atoms:
    """Tilt ``a3`` in plane by ``index`` steps, leaving the Cartesian positions alone.

    This is exactly what pei_vasp_run_gsfe does by rewriting POSCAR line 5: the atoms do
    not move, so the fault sits at the cell boundary and the slab stays a rigid shift of
    the perfect crystal rather than a homogeneous shear.

    Args:
        atoms_ref (Atoms): Image 00 (the unsheared reference).
        step_xy (np.ndarray): In-plane shear per index, shape (2,).
        index (int): Image index.

    Returns:
        Atoms: The sheared image.
    """
    atoms = atoms_ref.copy()
    cell = np.array(atoms.cell)
    cell[2, 0:2] = np.array(atoms_ref.cell)[2, 0:2] + step_xy * index
    atoms.set_cell(cell)
    return atoms


def get_gsfe_line(path_ydir: Path = None) -> tuple:
    """Read an existing ``y_dir`` and reconstruct the gamma-line it walks.

    The step is least-squares fitted through all existing images instead of being taken
    from a single pair, so a truncated ``bc`` digit in any one POSCAR cannot tilt the
    extension off the line.

    Args:
        path_ydir (Path): The ``y_gsfe_<mode>/y_dir`` directory.

    Returns:
        tuple:
            - lindex (list): Existing image indices, ascending.
            - atoms_ref (Atoms): Image 00.
            - step_xy (np.ndarray): In-plane shear per index, shape (2,).
    """
    lindex = sorted(int(p.name) for p in path_ydir.iterdir()
                    if p.is_dir() and p.name.isdigit())
    if len(lindex) < 2:
        fail("need at least two existing images in %s, found %d" % (path_ydir, len(lindex)))
    if lindex[0] != 0:
        fail("no reference image 00 in %s" % path_ydir)
    if lindex != list(range(len(lindex))):
        fail("images in %s are not a gap-free 00..%02d run: %s"
             % (path_ydir, lindex[-1], lindex))

    atoms_ref = read(str(path_ydir / "00" / "POSCAR"), format="vasp")
    a3_xy = np.array([np.array(read(str(path_ydir / ("%02d" % i) / "POSCAR"),
                                    format="vasp").cell)[2, 0:2] for i in lindex])

    # a3_xy(i) = a3_xy(0) + step * i, fitted over every image; the intercept is refitted
    # too so a shifted reference shows up in the residual instead of being absorbed.
    coef = np.polyfit(np.array(lindex, dtype=float), a3_xy, 1)
    step_xy, origin_xy = coef[0], coef[1]
    resid = np.abs(a3_xy - (origin_xy + np.outer(lindex, step_xy))).max()
    if resid > TOL_LINE:
        fail("existing images in %s are not on one evenly spaced line "
             "(max deviation %.3e Ang > %.1e); mymetal.post.gsfe would abort on them"
             % (path_ydir, resid, TOL_LINE))
    print("   line fit residual: %.3e Ang  ✅" % resid)
    return lindex, atoms_ref, step_xy


def get_gsfe_period(atoms_ref: Atoms = None, step_xy: np.ndarray = None,
                    index_probe_max: int = INDEX_PROBE_MAX) -> tuple:
    """Find the closing index and the symmetry-equivalent images inside one period.

    Args:
        atoms_ref (Atoms): Image 00.
        step_xy (np.ndarray): In-plane shear per index, shape (2,).
        index_probe_max (int): Largest index to probe before giving up.

    Returns:
        tuple:
            - index_period (int): Smallest index > 0 whose structure equals image 00.
            - lgroup (list): Groups of equivalent indices in ``[0, index_period)``.
    """
    dict_fp = {i: get_structure_fingerprint(shear_reference_cell(atoms_ref, step_xy, i))
               for i in range(index_probe_max + 1)}

    def same(i, j):
        fi, fj = dict_fp[i], dict_fp[j]
        return fi.shape == fj.shape and np.allclose(fi, fj, atol=FINGERPRINT_TOL)

    index_period = next((i for i in range(1, index_probe_max + 1) if same(i, 0)), None)
    if index_period is None:
        fail("no closing image within %d steps -- this slip vector never returns to the "
             "perfect crystal; give an explicit -index_end" % index_probe_max)

    lgroup, sindex_seen = [], set()
    for i in range(index_period):
        if i in sindex_seen:
            continue
        lgroup_now = [j for j in range(index_period)
                      if j not in sindex_seen and same(i, j)]
        sindex_seen.update(lgroup_now)
        lgroup.append(lgroup_now)
    return index_period, lgroup


def check_slip_mode(mode: str = None, atoms_ref: Atoms = None,
                    step_xy: np.ndarray = None, index_max: int = None) -> None:
    """Warn if the step on disk disagrees with the named mode's (bp1, bp2).

    Args:
        mode (str): Slip mode name, a key of :data:`DICT_SLIP_MODE`.
        atoms_ref (Atoms): Image 00.
        step_xy (np.ndarray): In-plane shear per index, shape (2,).
        index_max (int): Highest existing index, i.e. the ``N1`` the bash run used.
    """
    if mode not in DICT_SLIP_MODE:
        warn("slip mode %s is not in DICT_SLIP_MODE; skipping the (bp1, bp2) cross-check"
             % mode)
        return
    bp1, bp2 = DICT_SLIP_MODE[mode]
    cell = np.array(atoms_ref.cell)
    step_expect = np.array([cell[0, 0] * bp1, cell[1, 1] * bp2]) / index_max
    dev = np.abs(step_xy - step_expect).max()
    if dev > 1.0e-6:
        warn("step on disk %s does not match mode %s (bp1=%.4f, bp2=%.4f -> %s); "
             "is this y_dir really that mode?"
             % (np.round(step_xy, 6), mode, bp1, bp2, np.round(step_expect, 6)))
    else:
        print("   (bp1, bp2) cross-check against mode %s: ✅" % mode)


def get_index_todo(lindex: list = None, index_period: int = None, lgroup: list = None,
                   index_end: int = None, include_equivalent: bool = False) -> list:
    """Decide which image indices still have to be computed.

    Args:
        lindex (list): Existing indices.
        index_period (int): Closing index.
        lgroup (list): Equivalence groups inside one period.
        index_end (int): Last index to lay out; ``None`` means ``index_period - 1``.
        include_equivalent (bool): Keep images that only mirror an existing one.

    Returns:
        list: Indices to create, ascending.
    """
    if index_end is None:
        # index_period itself is the perfect crystal again, i.e. a duplicate of 00; the
        # gamma-line closes on 00's energy without spending a job on it.
        index_end = index_period - 1
    lindex_todo = [i for i in range(max(lindex) + 1, index_end + 1)]
    if include_equivalent:
        return lindex_todo

    sindex_covered = set(lindex)
    lindex_keep = []
    for i in lindex_todo:
        lgroup_i = next((g for g in lgroup if (i % index_period) in g), None)
        if lgroup_i is not None and sindex_covered & set(lgroup_i):
            continue
        lindex_keep.append(i)
        if lgroup_i is not None:
            sindex_covered.update(lgroup_i)
    return lindex_keep


def write_image(path_image: Path = None, path_src: Path = None,
                a3_new: np.ndarray = None) -> None:
    """Seed one image directory and point its ``a3`` at the new fault vector.

    The copied ``sub.*`` scripts are each stamped with their own
    ``#SBATCH --comment`` de-duplication tag, so ``pei_slurm_univ_sbatch_retry`` can tell
    "this image already has a job" from "this image still needs one". Without it every
    image would carry the byte-identical script copied out of ``y_full_relax`` and the
    wrapper would have nothing but WorkDir and a clock to go on.

    Args:
        path_image (Path): Directory to create.
        path_src (Path): ``y_full_relax``, the source of the VASP inputs.
        a3_new (np.ndarray): The new third lattice vector, shape (3,).
    """
    path_image.mkdir(parents=True)
    for name in LFILE_SEED:
        path_file = path_src / name
        if path_file.is_file():
            shutil.copy(str(path_file), str(path_image))
        else:
            warn("%s not found in %s" % (name, path_src))
    for pattern in LGLOB_SEED:
        for path_file in glob.glob(str(path_src / pattern)):
            shutil.copy(path_file, str(path_image))
            stamp_comment_tag(path_image / Path(path_file).name)

    path_poscar = path_image / "POSCAR"
    lline = path_poscar.read_text().splitlines(keepends=True)
    lline[4] = (FORMAT_A3 % (a3_new[0], a3_new[1], a3_new[2])) + "\n"
    path_poscar.write_text("".join(lline))


def write_sbatch_script(path_script: Path = None, path_ydir: Path = None,
                        lindex: list = None) -> None:
    """Write the submit list for the newly created images.

    Plain concatenation, never an f-string: the body is full of ``${}`` and ``$(())``
    that would have to be escaped everywhere.

    Args:
        path_script (Path): File to write.
        path_ydir (Path): The ``y_dir`` holding the images.
        lindex (list): Indices to submit.
    """
    text = ('#!/bin/bash\n'
            '# ================ ⚠️  AUTO-GENERATED by mymetal.build.workflow.gsfe\n'
            '# Do not edit -- pei_vasp_run_gsfe_extend overwrites this file on every run.\n'
            '# Submits ONLY the images this extension created; the old ones stay put.\n'
            '\n'
            'set -u\n'
            'path_ydir="' + str(path_ydir) + '"\n'
            '\n'
            'for dirn1 in ' + " ".join("%02d" % i for i in lindex) + ' ; do\n'
            '    cd "$path_ydir/$dirn1" || exit 1\n'
            '    echo "▶️  sbatch $path_ydir/$dirn1"\n'
            '    pei_slurm_univ_sbatch_retry sub.*\n'
            '    cd "$path_ydir" || exit 1\n'
            'done\n'
            'echo "🎉 submitted ' + str(len(lindex)) + ' image(s)"\n')
    path_script.write_text(text)
    path_script.chmod(0o755)


def generate_gsfe_extend_dirs(path_root: str = None, mode: str = None,
                              index_end: int = None, include_equivalent: bool = False,
                              srcdir: str = "y_full_relax", workdir: str = None,
                              index_probe_max: int = INDEX_PROBE_MAX) -> None:
    """Extend an existing GSFE gamma-line to a full slip period.

    Args:
        path_root (str): Directory holding ``y_full_relax`` and ``y_gsfe_<mode>``.
        mode (str): Slip mode name, e.g. ``FCC_111``.
        index_end (int): Last image index to lay out; ``None`` = one full period.
        include_equivalent (bool): Also lay out images that only mirror an existing one.
        srcdir (str): Relaxed reference directory name.
        workdir (str): GSFE directory name; ``None`` = ``y_gsfe_<mode>``.
        index_probe_max (int): Largest index probed when looking for the period.
    """
    # ================ 🔎 check -- structural preconditions only
    path_root = Path(path_root).resolve()
    if not path_root.is_dir():
        fail("path_root %s is not a directory" % path_root)
    path_src = path_root / srcdir
    if not path_src.is_dir():
        fail("reference directory %s not found" % path_src)
    path_work = path_root / (workdir if workdir else "y_gsfe_" + mode)
    path_ydir = path_work / "y_dir"
    if not path_ydir.is_dir():
        fail("%s not found -- run pei_vasp_run_gsfe first" % path_ydir)
    if index_end is not None and index_end < 1:
        fail("index_end must be a positive integer, got %s" % index_end)

    # ================ 🧱 prepare -- reconstruct the line, then decide what is missing
    print("📁 work dir  : %s" % path_work)
    print("📁 reference : %s" % path_src)
    lindex, atoms_ref, step_xy = get_gsfe_line(path_ydir)
    index_max = max(lindex)
    print("   existing images: 00..%02d  (step %s Ang, |step| = %.4f Ang)"
          % (index_max, np.round(step_xy, 6), np.linalg.norm(step_xy)))
    check_slip_mode(mode, atoms_ref, step_xy, index_max)

    index_period, lgroup = get_gsfe_period(atoms_ref, step_xy, index_probe_max)
    print("   full slip period: index %d  (%.4f Ang, %.2f x the existing segment)"
          % (index_period, index_period * np.linalg.norm(step_xy), index_period / index_max))
    lgroup_dup = [g for g in lgroup if len(g) > 1]
    if lgroup_dup:
        print("   symmetry-equivalent images inside one period: %s" % lgroup_dup)
    else:
        print("   symmetry-equivalent images inside one period: none, all %d differ"
              % index_period)

    lindex_todo = get_index_todo(lindex, index_period, lgroup, index_end, include_equivalent)
    if not lindex_todo:
        print("📊 nothing to do: 00..%02d already covers every inequivalent image of the "
              "full period." % index_max)
        print("🎉 done")
        return

    # ================ 🚀 main -- lay the new images out
    cell_ref = np.array(atoms_ref.cell)
    print("▶️  creating %d image(s): %s"
          % (len(lindex_todo), " ".join("%02d" % i for i in lindex_todo)))
    lindex_made, lindex_skip = [], []
    for index in lindex_todo:
        path_image = path_ydir / ("%02d" % index)
        if path_image.exists():
            warn("%s already exists, skipped" % path_image)
            lindex_skip.append(index)
            continue
        a3_new = cell_ref[2].copy()
        a3_new[0:2] = cell_ref[2, 0:2] + step_xy * index
        write_image(path_image, path_src, a3_new)
        print("   ✅ %02d  a3 in-plane = (%12.6f, %12.6f)  slip = %8.4f Ang"
              % (index, a3_new[0], a3_new[1], np.linalg.norm(step_xy) * index))
        lindex_made.append(index)

    path_script = path_work / "y_sbatch_extend.sh"
    if lindex_made:
        write_sbatch_script(path_script, path_ydir, lindex_made)

    # ================ 📊 summary
    print("📊 summary")
    print("   created : %d" % len(lindex_made))
    print("   skipped : %d%s" % (len(lindex_skip),
                                 ("  " + str(lindex_skip)) if lindex_skip else ""))
    print("   gamma-line now runs 00..%02d of a %d-step period"
          % (max(lindex_made + lindex_skip + lindex), index_period))
    if lindex_made:
        print("   submit with: bash %s" % path_script)
    print("🎉 done")
