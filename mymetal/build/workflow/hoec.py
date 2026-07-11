"""
hoec workflow module

Generate the deformed VASP input directories for HIGHER-ORDER elastic constants
(2nd + 3rd + 4th order) by the energy-strain method of Wang & Li, Phys. Rev. B 79,
224102 (2009), generalized to cubic (FCC) and hexagonal (HCP) symmetry.

Driven from a directory that contains ``y_full_relax`` (a fully relaxed CONVENTIONAL
cell: cubic 4-atom for FCC, hexagonal 2-atom for HCP), it builds::

    y_hoec_energy/
        y_hoec_modes.json                 # manifest read by mymetal.post.hoec_energy
        y_hoec_energy_<mode>/
            y_dir/
                s_+0.0000/  (POSCAR, INCAR, KPOINTS, POTCAR, Y_CONSTR_LATT)
                s_+0.0100/  ...

For each mode the Green-Lagrange strain is eta = xi * d (engineering Voigt), the
symmetric deformation gradient F = (I + 2 eta)^(1/2) is applied to the reference cell
(fractional coordinates fixed), and the ions are relaxed at fixed cell shape (ISIF=2).
:func:`mymetal.post.hoec_energy.post_hoec_energy` then fits U(xi)/V0 and solves for
the elastic constants.

This is the higher-order sibling of the ``pei_vasp_run_cij_energy`` bash workflow: same
``y_full_relax`` entry point and same ``y_dir/<strain>`` layout. It only prepares input
directories -- it writes no Slurm script and submits no job.

Functions:
    - prepare_hoec_reference: copy y_full_relax to a temp dir and make it an
      elastic-constant reference (CONTCAR, ISIF=2, ISYM=0, Y_CONSTR_LATT).
    - deform_atoms: apply a deformation gradient at fixed fractional coordinates.
    - generate_hoec_dirs: the full generation, writing the manifest.

Change log:
    - Written by J. P., 2026, as the standalone ``pei_vasp_run_hoec_energy`` script.
    - Merged into mymetal by J. P. on 2026-07-10; the group-theory core moved to
      ``mymetal.calculate.calmechanics.hoec``, the CLI stayed in
      ``vasp_utils/vasp_workflow_bulk/pei_vasp_run_hoec_energy``.
"""

import json
import os
import shutil
import subprocess
import numpy as np
from pathlib import Path
from ase import Atoms
from ase.io import read, write

from myvasp import vasp_func as vf

from mymetal.calculate.calmechanics.hoec import (
    MODES, check_symmetry, get_deformation_gradient, get_strain_list)
from mymetal.universal.print.print import fail, warn

# INCAR tags forced on every deformed job. ISIF=2 fixes the box and relaxes the ions
# (relaxed constants, as in the 2nd-order workflow); ISYM=0 keeps every mode on the
# same footing, since the modes lower the symmetry by different amounts.
INCAR_TAGS = {
    'isif': '2',
    'isym': '1',
    'lcharg': 'F',
    'algo': 'Normal',
    'nbands': 'comment',
}

# fully fix the box for the cluster VASP build (constrained-lattice file)
Y_CONSTR_LATT_TEXT = "0\n0.0  0.0  0.0  0.0  0.0  0.0\n"

# copied into every strain directory
COPY_FILES = ["INCAR", "KPOINTS", "POTCAR", "Y_CONSTR_LATT"]

FIND_AND_CHANGE = "pei_vasp_univ_find_and_change"


# ---------------------------------------------------------------------------
# reference preparation
# ---------------------------------------------------------------------------
def prepare_hoec_reference(path_full_relax: Path = None, path_tmp: Path = None) -> Path:
    """Copy ``y_full_relax`` to a temp dir and turn it into an elastic-constant reference.

    Ensures a CONTCAR exists, rewrites the INCAR tags in :data:`INCAR_TAGS`, and writes
    a fully-constrained ``Y_CONSTR_LATT``. Mirrors the preamble of
    ``pei_vasp_run_cij_energy``.

    ``pei_vasp_univ_find_and_change`` edits ``INCAR`` relative to the current directory,
    so this chdir's into the temp dir and always restores the original cwd.

    Args:
        path_full_relax (Path): The relaxed reference directory.
        path_tmp (Path): Temp directory to create (removed first if present).

    Returns:
        Path: Path to the reference CONTCAR inside ``path_tmp``.
    """
    if shutil.which(FIND_AND_CHANGE) is None:
        fail("%s not on PATH; source the vasp_utils environment first" % FIND_AND_CHANGE)

    if path_tmp.exists():
        shutil.rmtree(path_tmp)
    shutil.copytree(path_full_relax, path_tmp)

    path_contcar = path_tmp / "CONTCAR"
    if not path_contcar.is_file() or path_contcar.stat().st_size == 0:
        warn("CONTCAR missing/empty in reference; copying POSCAR -> CONTCAR")
        shutil.copy(path_tmp / "POSCAR", path_contcar)

    path_back = Path.cwd()
    os.chdir(path_tmp)
    try:
        for tag, value in INCAR_TAGS.items():
            subprocess.run([FIND_AND_CHANGE, "-" + tag, value], check=True)
    finally:
        os.chdir(path_back)      # every caller continues from a known cwd

    (path_tmp / "Y_CONSTR_LATT").write_text(Y_CONSTR_LATT_TEXT, encoding='utf-8')
    for f in COPY_FILES:
        if not (path_tmp / f).is_file():
            fail("required input file missing in reference: %s" % (path_tmp / f))
    return path_contcar


def deform_atoms(atoms_ref: Atoms = None, F: np.ndarray = None) -> Atoms:
    """Apply a deformation gradient, keeping fractional coordinates fixed.

    ASE stores lattice vectors as rows, so a column-vector map a' = F a becomes
    ``cell @ F.T`` on the row-vector cell.

    Args:
        atoms_ref (Atoms): Reference structure (not modified).
        F (np.ndarray): The 3x3 deformation gradient.

    Returns:
        Atoms: A deformed copy.
    """
    atoms = atoms_ref.copy()
    new_cell = np.array(atoms_ref.get_cell()) @ F.T
    atoms.set_cell(new_cell, scale_atoms=True)
    return atoms


# ---------------------------------------------------------------------------
# main generation
# ---------------------------------------------------------------------------
def generate_hoec_dirs(path_root: str = None, symmetry: str = 'auto',
                       emax: float = 0.12, de: float = 0.01,
                       srcdir: str = 'y_full_relax',
                       outdir: str = 'y_hoec_energy') -> Path:
    """Generate the deformed VASP input directories and the mode manifest.

    Args:
        path_root (str): Directory containing ``srcdir``. Defaults to the cwd.
        symmetry (str): 'auto', 'cubic' or 'hex'.
        emax (float): Maximum |xi|.
        de (float): Step in xi.
        srcdir (str): Reference directory name.
        outdir (str): Output directory name (removed first if present).

    Returns:
        Path: The output directory containing ``y_hoec_modes.json``.
    """
    ### check ------------------------------------------------------------------
    path_root = Path(path_root or Path.cwd()).resolve()
    if not path_root.is_dir():
        fail("root directory not found: %s" % path_root)
    path_full_relax = path_root / srcdir
    if not path_full_relax.is_dir():
        fail("%s not found under %s" % (srcdir, path_root))
    try:
        lxi = get_strain_list(emax, de)
    except ValueError as e:
        fail(str(e))
    path_tmp = path_root / (srcdir + "_hoec_temp")
    path_out = path_root / outdir

    ### prepare ----------------------------------------------------------------
    print("🧹 Preparing reference from %s" % path_full_relax.name)
    path_contcar = prepare_hoec_reference(path_full_relax, path_tmp)
    atoms_ref = read(path_contcar, format='vasp')
    cell_ref = np.array(atoms_ref.get_cell())
    v0 = float(abs(np.linalg.det(cell_ref)))
    natoms = len(atoms_ref)

    try:
        symmetry = check_symmetry(symmetry, cell_ref)
    except ValueError as e:
        fail(str(e))
    dict_modes = MODES[symmetry]

    print("📐 symmetry     : %s" % symmetry)
    print("📦 reference    : %d atoms, V0 = %.4f Å³" % (natoms, v0))
    print("🔧 modes        : %d (%s)" % (len(dict_modes), ", ".join(dict_modes)))
    print("📏 strain xi    : %d points in [%.4f, %.4f] step %.4f"
          % (len(lxi), -emax, emax, de))
    print("📊 total calcs  : %d modes × %d strains = %d single-point relaxations"
          % (len(dict_modes), len(lxi), len(dict_modes) * len(lxi)))
    warn("higher-order constants need CONVERGED dense k-points and ENCUT in "
         "%s (see Wang-Li Fig. 1,2). This does NOT change them." % srcdir)

    if path_out.exists():
        shutil.rmtree(path_out)
    path_out.mkdir()

    ### main: loop over modes and strains --------------------------------------
    dict_manifest_modes = {}
    for name, d_dir in dict_modes.items():
        path_ydir = path_out / ("y_hoec_energy_" + name) / "y_dir"
        path_ydir.mkdir(parents=True)
        print("\n================ 📁 mode %s  d=%s" % (name, d_dir))

        lstrain_dir = []
        for xi in lxi:
            label = "s_%+.4f" % xi
            path_strain = path_ydir / label
            path_strain.mkdir()
            try:
                F = get_deformation_gradient(d_dir, xi)
            except ValueError as e:
                fail(str(e))
            atoms = deform_atoms(atoms_ref, F)
            write(path_strain / "POSCAR", atoms, format='vasp',
                  direct=True, vasp5=True, sort=False)
            for f in COPY_FILES:
                shutil.copy(path_tmp / f, path_strain / f)
            lstrain_dir.append(label)
        print("   ▶️  wrote %d strained POSCARs" % len(lstrain_dir))
        print("   📁 prepared dir: %s" % path_ydir.parent)

        dict_manifest_modes[name] = dict(direction=list(d_dir), strain_dirs=lstrain_dir)
    # to here: every mode has a full y_dir of strained inputs

    ### write the manifest the post script reads -------------------------------
    # record the same eV/Å³ -> GPa constant mymetal.post uses, so the two never drift
    dict_manifest = dict(symmetry=symmetry, natoms=natoms, v0_ang3=v0,
                         ev_per_a3_to_gpa=vf.phy_const('qe') * 1e21,
                         emax=emax, de=de, xi=lxi,
                         ref_dir=str(path_full_relax), modes=dict_manifest_modes)
    (path_out / "y_hoec_modes.json").write_text(
        json.dumps(dict_manifest, indent=2), encoding='utf-8')

    shutil.rmtree(path_tmp)
    print("\n📊 SUMMARY: %s | %d modes | %d strains | %d calcs | %s"
          % (symmetry, len(dict_modes), len(lxi), len(dict_modes) * len(lxi),
             "generated only; no submit script, no sbatch"))
    print("🎉 done — manifest at %s" % (path_out / "y_hoec_modes.json"))
    return path_out
