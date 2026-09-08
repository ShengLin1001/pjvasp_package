"""Generate VASP KPAR/NCORE (and k-mesh) benchmark inputs from ``y_full_relax``.

The workflow keeps the VASP workload and the 128-rank Slurm launcher fixed while
varying only ``KPAR`` and ``NCORE``.  Every generated INCAR is changed through
``pei_vasp_univ_find_and_change``; ``NPAR`` is explicitly kept commented so it
cannot compete with the NCORE strategy.

A second, optional axis scans the **automatic k-mesh length** (``KPOINTS`` mode
``A``) and, with it, ``ENCUT``: the same cell is then timed at several k-point
densities, which is what answers "how much does the k mesh cost me?" before a
production campaign fixes its KPOINTS. Defaults reproduce the original
single-mesh benchmark (``A 40``, ``ENCUT 250``), so an existing call keeps its
meaning; pass ``lklength`` with the production density (and ``encut`` with the
production cutoff) to time the real workload instead of the cheap default one.
Atom-count scans are simply two runs with different ``srcdir``/``outdir``.

Functions:
    - get_default_pairs: Return the default 21 KPAR/NCORE combinations.
    - check_pairs: Validate benchmark pairs against the fixed MPI rank count.
    - check_klengths: Validate the automatic k-mesh lengths.
    - get_case_name: Build the ``kpar_<K>_ncore_<N>_klen_<L>`` case name.
    - get_kpoints_text: Render the automatic-mesh KPOINTS for one length.
    - write_case_control_files: Write the case KPOINTS and lattice constraint.
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
# Fixed part of the benchmark workload. ENCUT is deliberately NOT here: it is a
# knob (DEFAULT_ENCUT below) because a benchmark meant to predict a production
# campaign has to run at that campaign's cutoff, not at a cheap generic one.
FIXED_INCAR_TAGS = (
    ("nsw", "0"),
    ("ibrion", "-1"),
    ("isif", "2"),
    ("algo", "Normal"),
    ("lwave", "F"),
    ("lcharg", "F"),
    ("lelf", "F"),
)
DEFAULT_ENCUT = 250
DEFAULT_KLENGTH = 40
# None = 沿用源 INCAR 的 EDIFF。它不在 FIXED_INCAR_TAGS 里，是因为「收敛判据松一档能省
# 几步 SCF」本身就是要测的量，而不是 benchmark 的固定负载。
DEFAULT_EDIFF = None
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


def check_klengths(lklength: list[int]) -> list[int]:
    """Validate the automatic k-mesh lengths of the KPOINTS scan.

    Args:
        lklength (list[int]): Requested ``A <length>`` values.

    Returns:
        list[int]: Validated lengths, kept in the caller's order.
    """
    if not lklength:
        fail("at least one automatic k-mesh length is required")
    if len(set(lklength)) != len(lklength):
        fail("duplicate k-mesh lengths are not allowed: %s" % lklength)
    for klength in lklength:
        if int(klength) != klength or klength <= 0:
            fail("k-mesh length must be a positive integer: %s" % klength)
    return [int(klength) for klength in lklength]


def get_case_name(kpar: int, ncore: int, klength: int) -> str:
    """Build the benchmark case directory name.

    The k-mesh length is always part of the name so a finished tree stays
    self-describing: ``mymetal.post.kpar_ncore`` reads the three axes back out
    of it (and still accepts the legacy ``kpar_<K>_ncore_<N>`` name).

    Args:
        kpar (int): K-point parallelization groups.
        ncore (int): Cores working on one orbital.
        klength (int): Automatic k-mesh length.

    Returns:
        str: ``kpar_<K>_ncore_<N>_klen_<L>``.
    """
    return "kpar_%d_ncore_%d_klen_%d" % (kpar, ncore, klength)


def get_kpoints_text(klength: int) -> str:
    """Render the automatic-mesh KPOINTS content for one length.

    Args:
        klength (int): Automatic mesh length ``A <length>``.

    Returns:
        str: Complete KPOINTS file content.
    """
    return "Automatic mesh\n0\nA\n%d\n" % klength


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


def write_case_control_files(path_case: Path,
                             klength: int = DEFAULT_KLENGTH) -> None:
    """Write this case's automatic KPOINTS and the fully constrained lattice file.

    Args:
        path_case (Path): Benchmark case directory.
        klength (int): Automatic k-mesh length for this case.
    """
    # These files define the benchmark workload, so generate them here instead
    # of inheriting values that may vary between y_full_relax sources. KPOINTS
    # is the one that varies across cases when a k-mesh scan is requested.
    (path_case / "KPOINTS").write_text(get_kpoints_text(klength), encoding="utf-8")
    (path_case / "Y_CONSTR_LATT").write_text(
        Y_CONSTR_LATT_TEXT, encoding="utf-8")


def change_case_incar(path_case: Path, kpar: int, ncore: int,
                      encut: int = DEFAULT_ENCUT, ediff: str = DEFAULT_EDIFF) -> None:
    """Apply the fixed benchmark and parallel INCAR settings.

    Args:
        path_case (Path): Benchmark case containing INCAR.
        kpar (int): K-point parallelization groups.
        ncore (int): Cores working on one orbital.
        encut (int): Plane-wave cutoff shared by every case.
        ediff (str): SCF convergence criterion; ``None`` keeps the source value.
    """
    ltag = FIXED_INCAR_TAGS + (("encut", str(encut)),)
    if ediff is not None:
        ltag = ltag + (("ediff", str(ediff)),)
    for tag, value in ltag:
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


def check_case_incar(path_case: Path, kpar: int, ncore: int,
                     encut: int = DEFAULT_ENCUT, ediff: str = DEFAULT_EDIFF) -> None:
    """Verify all fixed and variable INCAR settings.

    Args:
        path_case (Path): Generated benchmark case.
        kpar (int): Expected KPAR.
        ncore (int): Expected NCORE.
        encut (int): Expected ENCUT.
        ediff (str): Expected EDIFF; ``None`` skips the check (value inherited).
    """
    path_incar = path_case / "INCAR"
    ltag = FIXED_INCAR_TAGS + (("encut", str(encut)),)
    if ediff is not None:
        ltag = ltag + (("ediff", str(ediff)),)
    for tag, expected in ltag:
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


def check_case_control_files(path_case: Path,
                             klength: int = DEFAULT_KLENGTH) -> None:
    """Verify the exact KPOINTS and Y_CONSTR_LATT benchmark content.

    Args:
        path_case (Path): Generated benchmark case.
        klength (int): Automatic k-mesh length this case was written with.
    """
    path_kpoints = path_case / "KPOINTS"
    path_constraint = path_case / "Y_CONSTR_LATT"
    if path_kpoints.read_text(encoding="utf-8") != get_kpoints_text(klength):
        fail("unexpected automatic-mesh KPOINTS content: %s" % path_kpoints)
    if path_constraint.read_text(encoding="utf-8") != Y_CONSTR_LATT_TEXT:
        fail("unexpected lattice constraint content: %s" % path_constraint)


def write_workflow_manifest(path_manifest: Path, path_source: Path,
                            ntasks: int,
                            lpairs: list[tuple[int, int]],
                            lklength: list[int] = None,
                            encut: int = DEFAULT_ENCUT,
                            ediff: str = DEFAULT_EDIFF) -> None:
    """Record the complete benchmark contract for later auditing.

    Args:
        path_manifest (Path): JSON output path.
        path_source (Path): Source ``y_full_relax`` directory.
        ntasks (int): Fixed MPI rank count.
        lpairs (list[tuple[int, int]]): Generated KPAR/NCORE pairs.
        lklength (list[int]): Automatic k-mesh lengths of the KPOINTS scan.
        encut (int): Plane-wave cutoff shared by every case.
        ediff (str): SCF convergence criterion, or ``None`` if inherited.
    """
    lklength = [DEFAULT_KLENGTH] if lklength is None else lklength
    dict_manifest = {
        "source_dir": str(path_source),
        "ntasks": ntasks,
        "pairs": [{"kpar": kpar, "ncore": ncore} for kpar, ncore in lpairs],
        "static": True,
        "npar_commented": True,
        "incar": dict({tag.upper(): value for tag, value in FIXED_INCAR_TAGS},
                      ENCUT=str(encut),
                      EDIFF="inherited" if ediff is None else str(ediff)),
        "kpoints": {"mode": "A", "lengths": lklength},
        "y_constr_latt": [0, 0, 0, 0, 0, 0, 0],
    }
    path_manifest.write_text(
        json.dumps(dict_manifest, indent=2) + "\n", encoding="utf-8")


def generate_kpar_ncore_dirs(path_root: str = None,
                             lpairs: list[tuple[int, int]] = None,
                             ntasks: int = DEFAULT_NTASKS,
                             srcdir: str = "y_full_relax",
                             outdir: str = "y_kpar_ncore",
                             lklength: list[int] = None,
                             encut: int = DEFAULT_ENCUT,
                             ediff: str = DEFAULT_EDIFF,
                             force: bool = False) -> Path:
    """Generate static KPAR/NCORE (and optional k-mesh) benchmark inputs.

    The generated tree is the full product ``pairs x k-mesh lengths``: with the
    default single length it is exactly the original KPAR/NCORE matrix, and with
    several lengths every pair is repeated at each k-point density.

    Args:
        path_root (str): Absolute directory containing ``srcdir``.
        lpairs (list[tuple[int, int]]): Pairs to generate. ``None`` uses the
            requested 21-pair default matrix.
        ntasks (int): Fixed MPI rank count for compatibility checks.
        srcdir (str): Relaxed reference directory name.
        outdir (str): Output workflow directory name.
        lklength (list[int]): Automatic k-mesh lengths (``KPOINTS`` mode ``A``)
            to scan. ``None`` keeps the single default length.
        encut (int): Plane-wave cutoff shared by every case; raise it to the
            production value when the benchmark has to predict a real campaign.
        ediff (str): SCF convergence criterion shared by every case; ``None``
            keeps whatever the source INCAR has. Scan it with one run per value
            (different ``outdir``) -- looser EDIFF changes the SCF step count,
            not the workload of one step.
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
    lklength = [DEFAULT_KLENGTH] if lklength is None else list(lklength)
    lklength = check_klengths(lklength)
    path_source = path_root / srcdir
    path_structure, lsubmit = check_source_inputs(path_source)
    path_out = path_root / outdir
    # existing output: ask before deleting (blank/No/no-tty aborts, -force skips)
    confirm_prepare_outdir(path_out, force=force)
    # to here: every structural prerequisite has been checked without writing

    ### prepare ----------------------------------------------------------------
    path_ydir = path_out / "y_dir"
    path_ydir.mkdir(parents=True)
    print("📁 source       : %s" % path_source)
    print("📄 structure    : %s -> POSCAR" % path_structure.name)
    print("🧮 fixed ranks  : %d" % ntasks)
    print("🔧 test pairs   : %d" % len(lpairs))
    print("🔭 k-mesh scan  : A %s" % ", A ".join(str(k) for k in lklength))
    print("⚡ encut        : %d" % encut)
    print("🎯 ediff        : %s" % ("inherited from source INCAR" if ediff is None else ediff))
    print("📁 output       : %d case(s) -> %s" % (len(lpairs) * len(lklength), path_out))
    # to here: the empty workflow root is ready

    ### main -------------------------------------------------------------------
    lcase = [(kpar, ncore, klength)
             for klength in lklength for kpar, ncore in lpairs]
    for index, (kpar, ncore, klength) in enumerate(lcase, start=1):
        name = get_case_name(kpar, ncore, klength)
        path_case = path_ydir / name
        print("\n================ ▶️  %02d/%02d %s"
              % (index, len(lcase), name))
        copy_case_inputs(path_source, path_structure, lsubmit, path_case)
        write_case_control_files(path_case, klength)
        change_case_incar(path_case, kpar, ncore, encut, ediff)
        check_case_incar(path_case, kpar, ncore, encut, ediff)
        check_case_control_files(path_case, klength)
        print("✅ fixed INCAR/KPOINTS/Y_CONSTR_LATT; KPAR=%d NCORE=%d A %d; "
              "NPAR remains commented" % (kpar, ncore, klength))
        print("📍 %s" % path_case)
    # to here: every benchmark directory contains a validated INCAR

    path_manifest = path_out / "y_kpar_ncore.json"
    write_workflow_manifest(path_manifest, path_source, ntasks, lpairs,
                            lklength=lklength, encut=encut, ediff=ediff)
    print("\n================ 📊 Summary")
    print("generated=%d    failed=0    submitted=0" % len(lcase))
    print("🎉 Done: inputs generated; manifest at %s" % path_manifest)
    return path_out
