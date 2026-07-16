"""Generate VASP KPAR/NCORE benchmark inputs from ``y_full_relax``.

The workflow keeps the VASP workload and the 128-rank Slurm launcher fixed while
varying only ``KPAR`` and ``NCORE``.  Every generated INCAR is changed through
``pei_vasp_univ_find_and_change``; ``NPAR`` is explicitly kept commented so it
cannot compete with the NCORE strategy.

Functions:
    - get_default_pairs: Return the default 21 KPAR/NCORE combinations.
    - check_pairs: Validate benchmark pairs against the fixed MPI rank count.
    - write_case_control_files: Write the fixed KPOINTS and lattice constraint.
    - write_workflow_manifest: Record the exact benchmark input contract.
    - generate_kpar_ncore_dirs: Build ``y_kpar_ncore/y_dir`` benchmark inputs.

The command-line front end lives in ``pei_vasp_run_kpar_ncore``; this module stays
import-only, exactly like ``mymetal.build.workflow.hoec``.
"""

import json
import shutil
import subprocess
from pathlib import Path

from mymetal.universal.print.print import fail, confirm_prepare_outdir


DEFAULT_NTASKS = 128
DEFAULT_KPAR_NCORE = {
    128: [1],
    64: [2, 1],
    32: [1, 2, 4],
    16: [1, 2, 4, 8],
    8: [1, 2, 4, 8, 16],
    4: [1, 2, 4, 8, 16, 32],
}
FIND_AND_CHANGE = "pei_vasp_univ_find_and_change"
REQUIRED_SOURCE_FILES = ("INCAR", "POTCAR")
FIXED_INCAR_TAGS = (
    ("nsw", "0"),
    ("ibrion", "-1"),
    ("isif", "2"),
    ("algo", "Normal"),
    ("encut", "250"),
    ("lwave", "F"),
    ("lcharg", "F"),
    ("lelf", "F"),
)
KPOINTS_TEXT = "Automatic mesh\n0\nA\n40\n"
Y_CONSTR_LATT_TEXT = "0\n0 0 0 0 0 0\n"


def get_default_pairs() -> list[tuple[int, int]]:
    """Return the default KPAR/NCORE combinations in plotting order.

    Returns:
        list[tuple[int, int]]: The 21 ``(KPAR, NCORE)`` pairs requested for a
        128-rank calculation.
    """
    lpairs = []
    for kpar, lncore in DEFAULT_KPAR_NCORE.items():
        for ncore in lncore:
            lpairs.append((kpar, ncore))
    return lpairs


def check_pairs(lpairs: list[tuple[int, int]], ntasks: int) -> list[tuple[int, int]]:
    """Validate benchmark pairs against the fixed MPI rank count.

    Args:
        lpairs (list[tuple[int, int]]): Requested KPAR/NCORE pairs.
        ntasks (int): MPI rank count used by every benchmark job.

    Returns:
        list[tuple[int, int]]: Validated pairs.
    """
    if ntasks <= 0:
        fail("ntasks must be a positive integer")
    if not lpairs:
        fail("at least one KPAR:NCORE pair is required")
    if len(set(lpairs)) != len(lpairs):
        fail("duplicate KPAR:NCORE pairs are not allowed")

    lkpar_supported = list(DEFAULT_KPAR_NCORE)
    for kpar, ncore in lpairs:
        if kpar not in lkpar_supported:
            fail("KPAR=%d is not one of the six plotted values: %s"
                 % (kpar, lkpar_supported))
        if ntasks % (kpar * ncore) != 0:
            fail("KPAR=%d, NCORE=%d is incompatible with %d MPI ranks; "
                 "KPAR*NCORE must divide ntasks" % (kpar, ncore, ntasks))
    return lpairs


def check_source_inputs(path_source: Path) -> tuple[Path, list[Path]]:
    """Check the reusable VASP inputs in ``y_full_relax``.

    Args:
        path_source (Path): Source calculation directory.

    Returns:
        tuple[Path, list[Path]]: Structure file to copy as POSCAR and Slurm
        scripts matching ``sub.*``.
    """
    if not path_source.is_dir():
        fail("source directory not found: %s" % path_source)

    for name in REQUIRED_SOURCE_FILES:
        path_input = path_source / name
        if not path_input.is_file() or path_input.stat().st_size == 0:
            fail("required source input missing or empty: %s" % path_input)

    path_structure = path_source / "CONTCAR"
    if not path_structure.is_file() or path_structure.stat().st_size == 0:
        path_structure = path_source / "POSCAR"
    if not path_structure.is_file() or path_structure.stat().st_size == 0:
        fail("neither a usable CONTCAR nor POSCAR exists in %s" % path_source)

    lsubmit = sorted(path_source.glob("sub.*"))
    lsubmit = [path_submit for path_submit in lsubmit if path_submit.is_file()]
    if not lsubmit:
        fail("no Slurm input matching sub.* found in %s" % path_source)
    return path_structure, lsubmit


def copy_case_inputs(path_source: Path, path_structure: Path,
                     lsubmit: list[Path], path_case: Path) -> None:
    """Copy one minimal, restart-free benchmark input set.

    Args:
        path_source (Path): Validated ``y_full_relax`` directory.
        path_structure (Path): CONTCAR or POSCAR used as the static POSCAR.
        lsubmit (list[Path]): Slurm scripts matching ``sub.*``.
        path_case (Path): New benchmark case directory.
    """
    path_case.mkdir()
    for name in REQUIRED_SOURCE_FILES:
        shutil.copy2(path_source / name, path_case / name)
    shutil.copy2(path_structure, path_case / "POSCAR")
    for path_submit in lsubmit:
        shutil.copy2(path_submit, path_case / path_submit.name)


def write_case_control_files(path_case: Path) -> None:
    """Write the fixed automatic KPOINTS and fully constrained lattice file.

    Args:
        path_case (Path): Benchmark case directory.
    """
    # These files define the common benchmark workload, so generate them here
    # instead of inheriting values that may vary between y_full_relax sources.
    (path_case / "KPOINTS").write_text(KPOINTS_TEXT, encoding="utf-8")
    (path_case / "Y_CONSTR_LATT").write_text(
        Y_CONSTR_LATT_TEXT, encoding="utf-8")


def change_case_incar(path_case: Path, kpar: int, ncore: int) -> None:
    """Apply the fixed benchmark and parallel INCAR settings.

    Args:
        path_case (Path): Benchmark case containing INCAR.
        kpar (int): K-point parallelization groups.
        ncore (int): Cores working on one orbital.
    """
    for tag, value in FIXED_INCAR_TAGS:
        subprocess.run(
            [FIND_AND_CHANGE, "-" + tag, value], cwd=path_case, check=True)

    # Keep these three explicit calls together: they are the variable under test,
    # and NPAR must remain disabled whenever the NCORE strategy is active.
    subprocess.run(
        [FIND_AND_CHANGE, "-kpar", str(kpar)], cwd=path_case, check=True)
    subprocess.run(
        [FIND_AND_CHANGE, "-ncore", str(ncore)], cwd=path_case, check=True)
    subprocess.run(
        [FIND_AND_CHANGE, "-npar", "comment"], cwd=path_case, check=True)


def get_active_incar_values(path_incar: Path, tag: str) -> list[str]:
    """Return active values for one INCAR tag, ignoring commented lines.

    Args:
        path_incar (Path): INCAR to inspect.
        tag (str): Case-insensitive tag name.

    Returns:
        list[str]: Active values found before inline comments.
    """
    lvalue = []
    for line in path_incar.read_text(encoding="utf-8").splitlines():
        stripped = line.lstrip()
        if not stripped or stripped.startswith("#") or "=" not in stripped:
            continue
        name, value = stripped.split("=", maxsplit=1)
        if name.strip().upper() == tag.upper():
            lvalue.append(value.split("#", maxsplit=1)[0].strip())
    return lvalue


def check_case_incar(path_case: Path, kpar: int, ncore: int) -> None:
    """Verify all fixed and variable INCAR settings.

    Args:
        path_case (Path): Generated benchmark case.
        kpar (int): Expected KPAR.
        ncore (int): Expected NCORE.
    """
    path_incar = path_case / "INCAR"
    for tag, expected in FIXED_INCAR_TAGS:
        lvalue = get_active_incar_values(path_incar, tag)
        if lvalue != [expected]:
            fail("%s has active %s values %s, expected [%s]"
                 % (path_incar, tag.upper(), lvalue, expected))

    lkpar = get_active_incar_values(path_incar, "KPAR")
    lncore = get_active_incar_values(path_incar, "NCORE")
    lnpar = get_active_incar_values(path_incar, "NPAR")
    if lkpar != [str(kpar)]:
        fail("%s has active KPAR values %s, expected [%d]" % (path_incar, lkpar, kpar))
    if lncore != [str(ncore)]:
        fail("%s has active NCORE values %s, expected [%d]"
             % (path_incar, lncore, ncore))
    if lnpar:
        fail("%s still has active NPAR values: %s" % (path_incar, lnpar))


def check_case_control_files(path_case: Path) -> None:
    """Verify the exact KPOINTS and Y_CONSTR_LATT benchmark content.

    Args:
        path_case (Path): Generated benchmark case.
    """
    path_kpoints = path_case / "KPOINTS"
    path_constraint = path_case / "Y_CONSTR_LATT"
    if path_kpoints.read_text(encoding="utf-8") != KPOINTS_TEXT:
        fail("unexpected automatic-mesh KPOINTS content: %s" % path_kpoints)
    if path_constraint.read_text(encoding="utf-8") != Y_CONSTR_LATT_TEXT:
        fail("unexpected lattice constraint content: %s" % path_constraint)


def write_workflow_manifest(path_manifest: Path, path_source: Path,
                            ntasks: int,
                            lpairs: list[tuple[int, int]]) -> None:
    """Record the complete benchmark contract for later auditing.

    Args:
        path_manifest (Path): JSON output path.
        path_source (Path): Source ``y_full_relax`` directory.
        ntasks (int): Fixed MPI rank count.
        lpairs (list[tuple[int, int]]): Generated KPAR/NCORE pairs.
    """
    dict_manifest = {
        "source_dir": str(path_source),
        "ntasks": ntasks,
        "pairs": [{"kpar": kpar, "ncore": ncore} for kpar, ncore in lpairs],
        "static": True,
        "npar_commented": True,
        "incar": {tag.upper(): value for tag, value in FIXED_INCAR_TAGS},
        "kpoints": {"mode": "A", "length": 40},
        "y_constr_latt": [0, 0, 0, 0, 0, 0, 0],
    }
    path_manifest.write_text(
        json.dumps(dict_manifest, indent=2) + "\n", encoding="utf-8")


def generate_kpar_ncore_dirs(path_root: str = None,
                             lpairs: list[tuple[int, int]] = None,
                             ntasks: int = DEFAULT_NTASKS,
                             srcdir: str = "y_full_relax",
                             outdir: str = "y_kpar_ncore",
                             force: bool = False) -> Path:
    """Generate static KPAR/NCORE benchmark inputs.

    Args:
        path_root (str): Absolute directory containing ``srcdir``.
        lpairs (list[tuple[int, int]]): Pairs to generate. ``None`` uses the
            requested 21-pair default matrix.
        ntasks (int): Fixed MPI rank count for compatibility checks.
        srcdir (str): Relaxed reference directory name.
        outdir (str): Output workflow directory name.
        force (bool): Delete an existing output directory without prompting.

    Returns:
        Path: Generated workflow directory.
    """
    ### check ------------------------------------------------------------------
    if path_root is None:
        fail("path_root is required")
    path_root = Path(path_root)
    if not path_root.is_absolute():
        fail("path_root must be absolute: %s" % path_root)
    path_root = path_root.resolve()
    if shutil.which(FIND_AND_CHANGE) is None:
        fail("%s not on PATH; source the vasp_utils environment first"
             % FIND_AND_CHANGE)

    lpairs = get_default_pairs() if lpairs is None else list(lpairs)
    lpairs = check_pairs(lpairs, ntasks)
    path_source = path_root / srcdir
    path_structure, lsubmit = check_source_inputs(path_source)
    path_out = path_root / outdir
    # existing output: ask before deleting (blank/No/no-tty aborts, --force skips)
    confirm_prepare_outdir(path_out, force=force)
    # to here: every structural prerequisite has been checked without writing

    ### prepare ----------------------------------------------------------------
    path_ydir = path_out / "y_dir"
    path_ydir.mkdir(parents=True)
    print("📁 source       : %s" % path_source)
    print("📄 structure    : %s -> POSCAR" % path_structure.name)
    print("🧮 fixed ranks  : %d" % ntasks)
    print("🔧 test pairs   : %d" % len(lpairs))
    print("📁 output       : %s" % path_out)
    # to here: the empty workflow root is ready

    ### main -------------------------------------------------------------------
    for index, (kpar, ncore) in enumerate(lpairs, start=1):
        name = "kpar_%d_ncore_%d" % (kpar, ncore)
        path_case = path_ydir / name
        print("\n================ ▶️  %02d/%02d %s"
              % (index, len(lpairs), name))
        copy_case_inputs(path_source, path_structure, lsubmit, path_case)
        write_case_control_files(path_case)
        change_case_incar(path_case, kpar, ncore)
        check_case_incar(path_case, kpar, ncore)
        check_case_control_files(path_case)
        print("✅ fixed INCAR/KPOINTS/Y_CONSTR_LATT; KPAR=%d NCORE=%d; "
              "NPAR remains commented" % (kpar, ncore))
        print("📍 %s" % path_case)
    # to here: every benchmark directory contains a validated INCAR

    path_manifest = path_out / "y_kpar_ncore.json"
    write_workflow_manifest(path_manifest, path_source, ntasks, lpairs)
    print("\n================ 📊 Summary")
    print("generated=%d    failed=0    submitted=0" % len(lpairs))
    print("🎉 Done: inputs generated; manifest at %s" % path_manifest)
    return path_out
