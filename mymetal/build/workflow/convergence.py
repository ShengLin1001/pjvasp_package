"""Generate VASP parameter-scan inputs from ``y_full_relax``, one axis per subtree.

The contract is the one every ``pei_vasp_run_*`` workflow follows: a case
inherits **everything** from ``y_full_relax`` and changes **only the one
parameter under test**.  Scanning KPAR/NCORE at the production cutoff is
therefore not a flag on this workflow -- it is a ``y_full_relax`` that already
carries the production ENCUT/EDIFF/KPOINTS, plus ``-kpar_ncore``.

Five axes are supported, each landing in its own subtree so a finished
``y_convergence`` is self-describing and the post-processor can dispatch on the
directories that actually exist::

    y_convergence/
        y_convergence_encuts/y_dir/encut_<E>
        y_convergence_kpoints/y_dir/kpoints_<n1>-<n2>-<n3>     explicit mesh
        y_convergence_klength/y_dir/klen_<L>                   automatic mesh, mode A
        y_convergence_ediff/y_dir/ediff_<D>
        y_convergence_kpar_ncore/y_dir/kpar_<K>_ncore_<N>

``kpoints`` and ``klength`` are two spellings of the same physical knob but stay
apart on purpose: the explicit mesh is a 3-vector and the automatic mesh is a
scalar, so they need different x axes when plotted, and keeping them separate
leaves existing ``y_convergence_kpoints`` trees readable by the old parser.

Every case is forced to a static single point (``FIXED_INCAR_TAGS`` plus a fully
constrained ``Y_CONSTR_LATT``); that is the fixed workload a scan is measured
against, not something inherited.  ``NPAR`` is kept commented so it cannot
compete with ``NCORE`` -- VASP derives it as ``ntasks / (KPAR * NCORE)`` anyway.

Functions:
    - get_default_pairs: Return the default KPAR/NCORE matrix.
    - get_axis_dirname: Map an axis name to its ``y_convergence_*`` directory.
    - get_case_name: Build the self-describing case directory name.
    - check_encuts: Validate the ENCUT scan values.
    - check_kpoints: Validate and group the explicit k-mesh triples.
    - check_klengths: Validate the automatic k-mesh lengths.
    - check_ediffs: Validate the EDIFF scan values.
    - check_pairs: Validate KPAR/NCORE pairs against the MPI rank count.
    - get_kpoints_text_auto: Render a mode-``A`` KPOINTS file.
    - get_kpoints_text_grid: Render an explicit Gamma-centred KPOINTS file.
    - build_axis_cases: Turn validated scan values into (case name, changes) pairs.
    - write_case_control_files: Write this case's KPOINTS and lattice constraint.
    - change_case_incar: Apply the fixed workload plus this case's INCAR changes.
    - check_case_incar: Verify every tag the generator claims to have set.
    - write_workflow_manifest: Record the exact scan contract as JSON.
    - generate_convergence_dirs: Build the whole ``y_convergence`` tree.

The command-line front end lives in ``pei_vasp_run_convergence``; this module
stays import-only, exactly like ``mymetal.build.workflow.kpar_ncore``.
"""

import json
import shutil
import subprocess
from pathlib import Path

from mymetal.universal.print.print import fail, warn, confirm_prepare_outdir
# Reused verbatim from the KPAR/NCORE workflow: source validation, the minimal
# restart-free copy, the INCAR reader, the fixed static workload and the
# fully-constrained lattice file are all axis-independent.
from mymetal.build.workflow.kpar_ncore import (
    FIND_AND_CHANGE, FIXED_INCAR_TAGS, Y_CONSTR_LATT_TEXT, DEFAULT_KPAR_NCORE,
    check_source_inputs, copy_case_inputs as copy_base_inputs,
    get_active_incar_values)


DEFAULT_NTASKS = 128
DEFAULT_SRCDIR = "y_full_relax"
DEFAULT_OUTDIR = "y_convergence"
# Per-axis default sweeps. Calling the workflow with no flags at all runs every axis
# below at these values -- it only writes input directories, so a wide default sweep
# costs nothing until the user picks which cases to submit.
#
# kpoints is deliberately absent: an explicit n1 n2 n3 mesh depends on the cell's shape
# and has no sensible universal default, so that axis only ever runs when asked for by
# name. The automatic mesh (klength) covers the same knob cell-independently.
DEFAULT_ENCUTS = [250, 300, 350, 400, 450, 500, 550, 600]
DEFAULT_KLENGTHS = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
DEFAULT_EDIFFS = ["1E-2", "1E-3", "1E-4", "1E-5", "1E-6",
                  "1E-7", "1E-8", "1E-9", "1E-10"]
# Axes run when the workflow is called with no flags at all.
DEFAULT_AXES = ["encuts", "klength", "ediff", "kpar_ncore"]
# One entry per scan axis. `dirname` is what pei_vasp_plot_convergence dispatches
# on, `prefix` builds the case name. Adding an axis means adding a row here plus a
# build_axis_cases branch -- nothing else in this module hard-codes an axis name.
DICT_AXIS = {
    "encuts":     {"dirname": "y_convergence_encuts",     "prefix": "encut"},
    "kpoints":    {"dirname": "y_convergence_kpoints",    "prefix": "kpoints"},
    "klength":    {"dirname": "y_convergence_klength",    "prefix": "klen"},
    "ediff":      {"dirname": "y_convergence_ediff",      "prefix": "ediff"},
    "kpar_ncore": {"dirname": "y_convergence_kpar_ncore", "prefix": "kpar"},
}


def get_default_pairs(ntasks: int = DEFAULT_NTASKS) -> list:
    """Return the default KPAR/NCORE matrix, inherited from the old benchmark.

    ``DEFAULT_KPAR_NCORE`` enumerates, for each of the six plotted KPAR values, the
    NCORE values with ``KPAR * NCORE <= ntasks``. It is reused verbatim rather than
    re-derived so the new workflow times exactly the same grid the old one did.

    Args:
        ntasks (int): MPI rank count; pairs needing more ranks are dropped.

    Returns:
        list[tuple[int, int]]: ``(KPAR, NCORE)`` pairs in plotting order.
    """
    lpair = []
    for kpar, lncore in DEFAULT_KPAR_NCORE.items():
        for ncore in lncore:
            if ntasks % (kpar * ncore) == 0:
                lpair.append((kpar, ncore))
    return lpair


def get_axis_dirname(axis: str = None) -> str:
    """Map an axis name to its ``y_convergence_*`` subtree name.

    Args:
        axis (str): One of the keys of ``DICT_AXIS``.

    Returns:
        str: Subtree directory name.
    """
    if axis not in DICT_AXIS:
        fail("unknown scan axis %r; supported: %s" % (axis, sorted(DICT_AXIS)))
    return DICT_AXIS[axis]["dirname"]


def get_case_name(axis: str = None, value=None) -> str:
    """Build the self-describing case directory name for one scan point.

    Self-describing rather than a bare value (the pre-2026-09 convention) so a
    case directory still says which axis it belongs to once it has been copied,
    quoted in a report, or collected into a cross-tree table.

    Args:
        axis (str): One of the keys of ``DICT_AXIS``.
        value: Scan value -- int for ``encuts``/``klength``, str for ``ediff``,
            3-list for ``kpoints``, ``(kpar, ncore)`` for ``kpar_ncore``.

    Returns:
        str: Case directory name, e.g. ``encut_400`` or ``kpar_4_ncore_4``.
    """
    prefix = DICT_AXIS[get_axis_key(axis)]["prefix"]
    if axis == "kpoints":
        return "%s_%s" % (prefix, "-".join(str(int(n)) for n in value))
    if axis == "kpar_ncore":
        return "%s_%d_ncore_%d" % (prefix, value[0], value[1])
    return "%s_%s" % (prefix, value)


def get_axis_key(axis: str = None) -> str:
    """Validate an axis name and return it unchanged.

    Args:
        axis (str): Candidate axis name.

    Returns:
        str: The same name, once known to be supported.
    """
    if axis not in DICT_AXIS:
        fail("unknown scan axis %r; supported: %s" % (axis, sorted(DICT_AXIS)))
    return axis


def check_encuts(lencut: list = None) -> list:
    """Validate the ENCUT scan values.

    Args:
        lencut (list): Requested plane-wave cutoffs in eV.

    Returns:
        list[int]: Validated cutoffs, in the caller's order.
    """
    if not lencut:
        fail("at least one ENCUT is required")
    if len(set(lencut)) != len(lencut):
        fail("duplicate ENCUT values are not allowed: %s" % lencut)
    for encut in lencut:
        if int(encut) != float(encut) or int(encut) <= 0:
            fail("ENCUT must be a positive integer (eV): %s" % encut)
    return [int(encut) for encut in lencut]


def check_kpoints(lkpoint: list = None) -> list:
    """Validate the explicit k-mesh values and group them into triples.

    Args:
        lkpoint (list): Flat list of mesh divisions, e.g. ``[3, 3, 3, 4, 4, 4]``.

    Returns:
        list[list[int]]: One ``[n1, n2, n3]`` per requested mesh.
    """
    if not lkpoint:
        fail("at least one explicit k mesh is required")
    if len(lkpoint) % 3 != 0:
        fail("explicit k-mesh list length (%d) is not divisible by 3; "
             "give whole n1 n2 n3 triples" % len(lkpoint))
    lgrid = []
    for index in range(0, len(lkpoint), 3):
        grid = []
        for value in lkpoint[index:index + 3]:
            if int(value) != float(value) or int(value) <= 0:
                fail("k-mesh division must be a positive integer: %s" % value)
            grid.append(int(value))
        lgrid.append(grid)
    if len(set(tuple(grid) for grid in lgrid)) != len(lgrid):
        fail("duplicate explicit k meshes are not allowed: %s" % lgrid)
    return lgrid


def check_klengths(lklength: list = None) -> list:
    """Validate the automatic k-mesh lengths (``KPOINTS`` mode ``A``).

    Args:
        lklength (list): Requested ``A <length>`` values.

    Returns:
        list[int]: Validated lengths, in the caller's order.
    """
    if not lklength:
        fail("at least one automatic k-mesh length is required")
    if len(set(lklength)) != len(lklength):
        fail("duplicate k-mesh lengths are not allowed: %s" % lklength)
    for klength in lklength:
        if int(klength) != float(klength) or int(klength) <= 0:
            fail("k-mesh length must be a positive integer: %s" % klength)
    return [int(klength) for klength in lklength]


def check_ediffs(lediff: list = None) -> list:
    """Validate the EDIFF scan values.

    Kept as the caller's own strings so the INCAR and the case name show exactly
    what was asked for (``1E-6``, not ``1e-06``) -- the scan is easier to read
    back that way, and ``pei_vasp_univ_find_and_change`` writes the string
    verbatim.

    Args:
        lediff (list): Requested SCF criteria, e.g. ``["1E-6", "1E-8"]``.

    Returns:
        list[str]: Validated criteria, in the caller's order.
    """
    if not lediff:
        fail("at least one EDIFF is required")
    lvalue = [str(ediff).strip() for ediff in lediff]
    if len(set(lvalue)) != len(lvalue):
        fail("duplicate EDIFF values are not allowed: %s" % lvalue)
    for ediff in lvalue:
        try:
            value = float(ediff)
        except ValueError:
            fail("EDIFF must be a number: %s" % ediff)
        if value <= 0:
            fail("EDIFF must be positive: %s" % ediff)
    return lvalue


def check_pairs(lpair: list = None, ntasks: int = DEFAULT_NTASKS) -> list:
    """Validate KPAR/NCORE pairs against the fixed MPI rank count.

    Args:
        lpair (list[tuple[int, int]]): Requested ``(KPAR, NCORE)`` pairs.
        ntasks (int): MPI rank count used by every case.

    Returns:
        list[tuple[int, int]]: Validated pairs, in the caller's order.
    """
    if ntasks <= 0:
        fail("ntasks must be a positive integer")
    if not lpair:
        fail("at least one KPAR:NCORE pair is required")
    lpair = [(int(kpar), int(ncore)) for kpar, ncore in lpair]
    if len(set(lpair)) != len(lpair):
        fail("duplicate KPAR:NCORE pairs are not allowed: %s" % lpair)

    # KPAR outside the six values the old benchmark plotted is a warning, not an error:
    # KPAR 1 and 2 are perfectly legal VASP and are the only way to reach the very high
    # NCORE end (2:64, 1:128). The plot just gets a curve it has no preset colour for.
    lkpar_plotted = list(DEFAULT_KPAR_NCORE)
    for kpar, ncore in lpair:
        if kpar * ncore > ntasks:
            fail("KPAR=%d, NCORE=%d needs %d ranks but only %d are available; "
                 "KPAR*NCORE must not exceed ntasks"
                 % (kpar, ncore, kpar * ncore, ntasks))
        if ntasks % (kpar * ncore) != 0:
            fail("KPAR=%d, NCORE=%d is incompatible with %d MPI ranks; "
                 "KPAR*NCORE must divide ntasks" % (kpar, ncore, ntasks))
        if kpar not in lkpar_plotted:
            warn("KPAR=%d is outside the six plotted values %s; the case is generated "
                 "but the KPAR/NCORE figure has no preset curve for it"
                 % (kpar, lkpar_plotted))
    return lpair


def get_kpoints_text_auto(klength: int = None) -> str:
    """Render a mode-``A`` (automatic mesh) KPOINTS file.

    Args:
        klength (int): Automatic mesh length.

    Returns:
        str: Complete KPOINTS content.
    """
    return "Automatic mesh\n0\nA\n%d\n" % klength


def get_kpoints_text_grid(lgrid: list = None) -> str:
    """Render an explicit Gamma-centred KPOINTS file.

    Args:
        lgrid (list[int]): ``[n1, n2, n3]`` mesh divisions.

    Returns:
        str: Complete KPOINTS content.
    """
    return ("Regular %d x %d x %d mesh centered at Gamma\n0\nGamma\n"
            "%d %d %d\n0  0  0\n"
            % (lgrid[0], lgrid[1], lgrid[2], lgrid[0], lgrid[1], lgrid[2]))


def build_axis_cases(axis: str = None, lvalue: list = None) -> list:
    """Turn one axis's validated scan values into concrete case descriptions.

    A case description is ``(name, ltag, kpoints_text)``:

    - ``ltag``          INCAR tags this case changes, on top of the fixed workload;
    - ``kpoints_text``  KPOINTS content, or ``None`` to inherit ``y_full_relax``'s.

    Everything not named here is inherited, which is the whole contract.

    Args:
        axis (str): One of the keys of ``DICT_AXIS``.
        lvalue (list): Values already through the matching ``check_*``.

    Returns:
        list[tuple[str, list, str]]: One entry per case.
    """
    axis = get_axis_key(axis)
    lcase = []
    for value in lvalue:
        name = get_case_name(axis, value)
        if axis == "encuts":
            lcase.append((name, [("encut", str(value))], None))
        elif axis == "ediff":
            lcase.append((name, [("ediff", str(value))], None))
        elif axis == "kpoints":
            lcase.append((name, [], get_kpoints_text_grid(value)))
        elif axis == "klength":
            lcase.append((name, [], get_kpoints_text_auto(value)))
        elif axis == "kpar_ncore":
            # NPAR must stay commented while the NCORE strategy is active: the two
            # are the same split written two ways, and an active NPAR would fight it.
            lcase.append((name, [("kpar", str(value[0])), ("ncore", str(value[1])),
                                 ("npar", "comment")], None))
    return lcase


def copy_case_inputs(path_source: Path = None, path_structure: Path = None,
                     lsubmit: list = None, path_case: Path = None) -> None:
    """Copy one restart-free case input set, inheriting the source KPOINTS.

    ``mymetal.build.workflow.kpar_ncore.copy_case_inputs`` brings INCAR, POTCAR,
    the structure and the Slurm scripts; that workflow always generated KPOINTS
    from a flag, so it deliberately did not copy one. Here KPOINTS is inherited
    like everything else and only the two k axes overwrite it afterwards.

    Args:
        path_source (Path): Validated ``y_full_relax`` directory.
        path_structure (Path): CONTCAR or POSCAR used as the static POSCAR.
        lsubmit (list[Path]): Slurm scripts matching ``sub.*``.
        path_case (Path): New case directory.
    """
    copy_base_inputs(path_source, path_structure, lsubmit, path_case)
    path_kpoints = path_source / "KPOINTS"
    if path_kpoints.is_file():
        shutil.copy2(path_kpoints, path_case / "KPOINTS")


def write_case_control_files(path_case: Path = None,
                             kpoints_text: str = None) -> None:
    """Write this case's KPOINTS (when it owns one) and the lattice constraint.

    Args:
        path_case (Path): Case directory.
        kpoints_text (str): KPOINTS content, or ``None`` to keep the inherited file.
    """
    # Y_CONSTR_LATT is part of the fixed static workload, so it is generated for
    # every case rather than inherited -- a y_full_relax may well carry a
    # relaxation constraint that would silently change what is being measured.
    (path_case / "Y_CONSTR_LATT").write_text(Y_CONSTR_LATT_TEXT, encoding="utf-8")
    if kpoints_text is not None:
        (path_case / "KPOINTS").write_text(kpoints_text, encoding="utf-8")


def change_case_incar(path_case: Path = None, ltag: list = None) -> None:
    """Apply the fixed static workload plus this case's own INCAR changes.

    Args:
        path_case (Path): Case directory containing INCAR.
        ltag (list[tuple[str, str]]): ``(tag, value)`` pairs for this case;
            ``("npar", "comment")`` comments the tag out instead of setting it.
    """
    for tag, value in list(FIXED_INCAR_TAGS) + list(ltag):
        subprocess.run(
            [FIND_AND_CHANGE, "-" + tag, value], cwd=path_case, check=True)


def check_case_incar(path_case: Path = None, ltag: list = None) -> None:
    """Verify every INCAR tag this generator claims to have set.

    Only generated tags are checked; inherited ones are the source's business and
    checking them here would turn an intentional y_full_relax choice into an error.

    Args:
        path_case (Path): Generated case directory.
        ltag (list[tuple[str, str]]): This case's own ``(tag, value)`` pairs.
    """
    path_incar = path_case / "INCAR"
    for tag, expected in list(FIXED_INCAR_TAGS) + list(ltag):
        lvalue = get_active_incar_values(path_incar, tag)
        if expected == "comment":
            if lvalue:
                fail("%s still has active %s values: %s"
                     % (path_incar, tag.upper(), lvalue))
        elif lvalue != [expected]:
            fail("%s has active %s values %s, expected [%s]"
                 % (path_incar, tag.upper(), lvalue, expected))


def check_case_control_files(path_case: Path = None,
                             kpoints_text: str = None) -> None:
    """Verify the generated KPOINTS and lattice constraint content.

    Args:
        path_case (Path): Generated case directory.
        kpoints_text (str): KPOINTS content this case was written with, or ``None``.
    """
    path_constraint = path_case / "Y_CONSTR_LATT"
    if path_constraint.read_text(encoding="utf-8") != Y_CONSTR_LATT_TEXT:
        fail("unexpected lattice constraint content: %s" % path_constraint)
    if kpoints_text is not None:
        path_kpoints = path_case / "KPOINTS"
        if path_kpoints.read_text(encoding="utf-8") != kpoints_text:
            fail("unexpected KPOINTS content: %s" % path_kpoints)


def write_workflow_manifest(path_manifest: Path = None, path_source: Path = None,
                            ntasks: int = DEFAULT_NTASKS,
                            dict_scan: dict = None) -> None:
    """Record the exact scan contract for later auditing.

    Args:
        path_manifest (Path): JSON output path.
        path_source (Path): Source ``y_full_relax`` directory.
        ntasks (int): Fixed MPI rank count.
        dict_scan (dict): Axis name -> validated values actually generated.
    """
    dict_manifest = {
        "source_dir": str(path_source),
        "ntasks": ntasks,
        "static": True,
        "inherits_from_source": "everything not listed under 'scan' or 'incar_fixed'",
        "incar_fixed": {tag.upper(): value for tag, value in FIXED_INCAR_TAGS},
        "y_constr_latt": [0, 0, 0, 0, 0, 0, 0],
        "scan": {axis: [list(value) if isinstance(value, (list, tuple)) else value
                        for value in lvalue]
                 for axis, lvalue in dict_scan.items()},
    }
    path_manifest.write_text(
        json.dumps(dict_manifest, indent=2) + "\n", encoding="utf-8")


def generate_convergence_dirs(path_root: str = None,
                              lencut: list = None,
                              lkpoint: list = None,
                              lklength: list = None,
                              lediff: list = None,
                              lpair: list = None,
                              ntasks: int = DEFAULT_NTASKS,
                              srcdir: str = DEFAULT_SRCDIR,
                              outdir: str = DEFAULT_OUTDIR,
                              force: bool = False,
                              if_append: bool = False) -> Path:
    """Generate a static parameter-scan tree from ``y_full_relax``.

    Each requested axis gets its own ``y_convergence_*`` subtree; axes left as
    ``None`` are simply not created, which is what lets the post-processor
    dispatch on the directories that exist.

    Args:
        path_root (str): Absolute directory containing ``srcdir``.
        lencut (list): ENCUT values to scan, or ``None``.
        lkpoint (list): Flat ``n1 n2 n3 ...`` explicit meshes, or ``None``.
        lklength (list): Automatic mesh lengths (mode ``A``), or ``None``.
        lediff (list): EDIFF values to scan, or ``None``.
        lpair (list[tuple[int, int]]): ``(KPAR, NCORE)`` pairs, or ``None``.
        ntasks (int): Fixed MPI rank count, used to reject impossible pairs.
        srcdir (str): Relaxed reference directory name.
        outdir (str): Output workflow directory name.
        force (bool): Delete an existing output directory without prompting.
        if_append (bool): Add cases to an existing tree instead of replacing it.
            Case directories that already exist are left untouched, so adding a
            few points to a finished scan never re-runs the ones already done.

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

    # Validate every requested axis before writing anything, so a typo in the last
    # flag cannot leave half a tree on disk.
    dict_scan = {}
    if lencut is not None:
        dict_scan["encuts"] = check_encuts(lencut)
    if lkpoint is not None:
        dict_scan["kpoints"] = check_kpoints(lkpoint)
    if lklength is not None:
        dict_scan["klength"] = check_klengths(lklength)
    if lediff is not None:
        dict_scan["ediff"] = check_ediffs(lediff)
    if lpair is not None:
        dict_scan["kpar_ncore"] = check_pairs(lpair, ntasks)
    if not dict_scan:
        fail("no scan axis requested; give at least one of "
             "-encuts / -kpoints / -klengths / -ediffs / -kpar_ncore")

    path_source = path_root / srcdir
    path_structure, lsubmit = check_source_inputs(path_source)
    # Only the two k axes write their own KPOINTS; every other axis inherits the
    # source's. Catch a missing source KPOINTS here rather than letting VASP fall
    # over one case at a time -- the old kpar_ncore workflow always generated a
    # mesh from a flag, so a y_full_relax carried over from it may not have one.
    linherit = sorted(axis for axis in dict_scan
                      if axis not in ("kpoints", "klength"))
    if linherit and not (path_source / "KPOINTS").is_file():
        fail("%s has no KPOINTS, but the %s axis/axes inherit it. Put the "
             "production KPOINTS in %s first -- under this workflow the source "
             "IS the reference setting."
             % (path_source, ", ".join(linherit), srcdir))
    path_out = path_root / outdir
    if if_append:
        if not path_out.is_dir():
            fail("-if_append needs an existing %s to add to; there is none at %s"
                 % (outdir, path_out))
    else:
        # existing output: ask before deleting (blank/No/no-tty aborts, -force skips)
        confirm_prepare_outdir(path_out, force=force)
    # to here: every structural prerequisite has been checked without writing

    ### prepare ----------------------------------------------------------------
    path_out.mkdir(parents=True, exist_ok=if_append)
    dict_axis_cases = {axis: build_axis_cases(axis, lvalue)
                       for axis, lvalue in dict_scan.items()}
    ncase = sum(len(lcase) for lcase in dict_axis_cases.values())
    print("📁 source       : %s" % path_source)
    print("📄 structure    : %s -> POSCAR" % path_structure.name)
    print("🧮 fixed ranks  : %d" % ntasks)
    print("🔭 scan axes    : %d (%s)"
          % (len(dict_scan), ", ".join(sorted(dict_scan))))
    for axis, lvalue in sorted(dict_scan.items()):
        print("   ➤ %-11s: %s" % (axis, lvalue))
    print("📁 output       : %d case(s) -> %s" % (ncase, path_out))
    print("♻️  inherited    : everything else comes from %s" % srcdir)
    # to here: the empty workflow root is ready

    ### main -------------------------------------------------------------------
    index = 0
    nskip = 0
    for axis in sorted(dict_axis_cases):
        path_axis = path_out / get_axis_dirname(axis)
        path_ydir = path_axis / "y_dir"
        path_ydir.mkdir(parents=True, exist_ok=if_append)
        print("\n================ 📁 %s (%d case(s))"
              % (path_axis.name, len(dict_axis_cases[axis])))
        for name, ltag, kpoints_text in dict_axis_cases[axis]:
            index += 1
            path_case = path_ydir / name
            print("\n================ ▶️  %02d/%02d %s/%s"
                  % (index, ncase, path_axis.name, name))
            if if_append and path_case.is_dir():
                nskip += 1
                print("⏭️  already present, left untouched")
                print("📍 %s" % path_case)
                continue
            copy_case_inputs(path_source, path_structure, lsubmit, path_case)
            write_case_control_files(path_case, kpoints_text)
            change_case_incar(path_case, ltag)
            check_case_incar(path_case, ltag)
            check_case_control_files(path_case, kpoints_text)
            print("✅ static workload fixed; this case changes: %s"
                  % (", ".join("%s=%s" % (tag.upper(), value)
                               for tag, value in ltag)
                     if ltag else "KPOINTS only"))
            print("📍 %s" % path_case)
    # to here: every case directory contains a validated INCAR
    path_manifest = path_out / "y_convergence.json"
    if if_append and path_manifest.is_file():
        # Appending must not erase what the tree already records, or the manifest
        # would claim the tree holds only the points added last.
        dict_old = json.loads(path_manifest.read_text(encoding="utf-8"))
        for axis, lvalue in dict_old.get("scan", {}).items():
            lkeep = [value for value in lvalue
                     if value not in dict_scan.get(axis, [])
                     and tuple(value) not in
                     [tuple(v) if isinstance(v, (list, tuple)) else v
                      for v in dict_scan.get(axis, [])]]
            dict_scan[axis] = lkeep + list(dict_scan.get(axis, []))
    write_workflow_manifest(path_manifest, path_source, ntasks, dict_scan)
    print("\n================ 📊 Summary")
    print("generated=%d    skipped=%d    failed=0    submitted=0"
          % (ncase - nskip, nskip))
    print("🎉 Done: inputs generated (this workflow never submits); "
          "manifest at %s" % path_manifest)
    return path_out
