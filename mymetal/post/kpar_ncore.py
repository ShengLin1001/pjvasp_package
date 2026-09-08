"""Post-process and plot VASP KPAR/NCORE and k-mesh timing/energy benchmarks.

``pei_vasp_univ_post`` first extracts one ``Elapsedtime`` value per benchmark
directory into ``y_post_time.txt`` and ``energy(sigma->0)`` into
``y_post_data.txt``.  This module reads both quantities, adds the per-case cost
drivers scraped from OUTCAR (``NKPTS`` / ``NBANDS`` / SCF steps / mean ``LOOP:``
time), writes a sortable table, and draws timing plus relative-energy panels.

Benchmark cases are keyed by ``(KPAR, NCORE, klength)``. A tree generated
before the k-mesh axis existed has no ``_klen_<L>`` suffix; those labels are
read as ``LEGACY_KLENGTH`` so old trees still post-process unchanged. With a
single k length the outputs are exactly the original ones; with several, each
length also gets its own KPAR/NCORE figure plus one k-mesh scan figure, which
is what answers "what does the k mesh cost, and is the cheaper mesh converged".

Functions:
    - run_univ_post: Generate the standard ``y_post_*`` files.
    - read_case_klength: Read the automatic k-mesh length from a case KPOINTS.
    - parse_case_label: Split a case directory name into its three axes.
    - read_kpar_ncore_times: Parse benchmark labels and elapsed seconds.
    - read_kpar_ncore_energies: Parse benchmark labels and total energies.
    - read_vasp_natoms: Read atom counts without requiring chemical symbols.
    - get_kpar_ncore_natoms: Read the benchmark structure's atom count.
    - read_outcar_cost: Scrape NKPTS/NBANDS/SCF-step cost drivers from OUTCAR.
    - get_kpar_ncore_costs: Collect those cost drivers for every case.
    - get_kpar_ncore_forces: Collect the final forces of every case.
    - read_outcar_forces: Read the final atomic forces of one case.
    - get_force_stats: Force magnitudes and their deviation from a reference case.
    - read_vasp_case: Read one finished static run as a single benchmark row.
    - collect_vasp_cases: Build one cross-tree table of settings vs E/F/time.
    - get_delta_energies: Convert total energies to relative meV/atom.
    - write_kpar_ncore_times: Write timing and energy data plus the fastest pair.
    - post_kpar_ncore: Run the complete post-processing workflow.

The command-line front end lives in ``pei_vasp_plot_kpar_ncore``; this module stays
import-only, exactly like ``mymetal.post.hoec_energy``.
"""

import re
import shutil
import subprocess
from pathlib import Path

from mymetal.universal.print.print import fail, warn
from mymetal.universal.plot.workflow import my_plot_kpar_ncore, my_plot_kpar_ncore_klen


KPAR_ORDER = [128, 64, 32, 16, 8, 4]
EXPECTED_PAIRS = {
    128: [1],
    64: [1, 2],
    32: [1, 2, 4],
    16: [1, 2, 4, 8],
    8: [1, 2, 4, 8, 16],
    4: [1, 2, 4, 8, 16, 32],
}
UNIV_POST = "pei_vasp_univ_post"
# The k-mesh suffix is optional: trees generated before that axis existed carry
# bare kpar_<K>_ncore_<N> names, and their single mesh was A 40.
CASE_PATTERN = re.compile(r"^kpar_(\d+)_ncore_(\d+)(?:_klen_(\d+))?$")
LEGACY_KLENGTH = 40
FLOAT_PATTERN = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
TIME_PATTERN = re.compile(
    r"^(\S+).*?Elapsedtime\s+(" + FLOAT_PATTERN + r")")
ENERGY_PATTERN = re.compile(
    r"^(\S+)\s+(" + FLOAT_PATTERN + r")(?:\s|$)")
# OUTCAR cost drivers: NKPTS/NBANDS set the size of the diagonalisation, the
# LOOP: lines time one electronic step and LOOP+: one ionic step.
NKPTS_PATTERN = re.compile(r"NKPTS\s*=\s*(\d+).*?NBANDS\s*=\s*(\d+)")
LOOP_PATTERN = re.compile(r"^\s*LOOP:.*real time\s+(" + FLOAT_PATTERN + r")")
LOOP_PLUS_PATTERN = re.compile(r"^\s*LOOP\+:.*real time\s+(" + FLOAT_PATTERN + r")")
ELAPSED_PATTERN = re.compile(
    r"Elapsed time \(sec\):\s*(" + FLOAT_PATTERN + r")")
SIGMA0_PATTERN = re.compile(
    r"energy\(sigma->0\)\s*=\s*(" + FLOAT_PATTERN + r")")
# 全自动 KPOINTS：第 3 行 A/Auto，长度可以跟在同一行，也可以单独占第 4 行
KLENGTH_PATTERN = re.compile(
    r"^\s*[Aa](?:uto\w*)?\s*(" + FLOAT_PATTERN + r")?\s*$")


def run_univ_post(path_workflow: Path) -> None:
    """Run ``pei_vasp_univ_post`` in the benchmark workflow directory.

    Args:
        path_workflow (Path): Directory containing ``y_dir``.
    """
    if shutil.which(UNIV_POST) is None:
        fail("%s not on PATH; source the vasp_utils environment first" % UNIV_POST)
    print("▶️  running %s in %s" % (UNIV_POST, path_workflow))
    subprocess.run([UNIV_POST], cwd=path_workflow, check=True)


def read_case_klength(path_case: Path = None) -> int:
    """Read the automatic k-mesh length out of a case's own KPOINTS.

    Args:
        path_case (Path): Case directory.

    Returns:
        int: Mesh length, or ``None`` if there is no mode-``A`` KPOINTS to read.
    """
    path_kpoints = path_case / "KPOINTS"
    if not path_kpoints.is_file():
        return None
    lline = path_kpoints.read_text(encoding="utf-8").splitlines()
    if len(lline) > 3 and lline[2].strip().upper().startswith("A"):
        try:
            return int(lline[3].split()[0])
        except (ValueError, IndexError):
            return None
    return None


def parse_case_label(label: str, path_ydir: Path = None) -> tuple[int, int, int]:
    """Split a benchmark case directory name into its three axes.

    Two naming generations reach here. The middle one encoded the mesh in the name
    (``..._klen_81``). The current convergence workflow does not: KPOINTS is
    inherited from ``y_full_relax`` like everything else, so the mesh is read from
    the case's own KPOINTS. Only a case that offers neither -- the oldest trees,
    whose single mesh really was A 40 -- falls back to ``LEGACY_KLENGTH``.

    Args:
        label (str): Case name, ``kpar_<K>_ncore_<N>[_klen_<L>]``.
        path_ydir (Path): ``y_dir`` holding the case, used to read KPOINTS when the
            name carries no mesh. ``None`` skips straight to the legacy default.

    Returns:
        tuple[int, int, int]: ``(KPAR, NCORE, klength)``, or ``None`` if the label
        is not a benchmark case at all.
    """
    match_case = CASE_PATTERN.fullmatch(label)
    if match_case is None:
        return None
    klength = match_case.group(3)
    if klength is not None:
        klength = int(klength)
    elif path_ydir is not None:
        klength = read_case_klength(path_ydir / label)
    return (int(match_case.group(1)), int(match_case.group(2)),
            LEGACY_KLENGTH if klength is None else klength)


def read_kpar_ncore_times(path_time: Path) -> dict[tuple[int, int], float]:
    """Read elapsed seconds from a standard ``y_post_time.txt``.

    Args:
        path_time (Path): File generated by ``pei_vasp_univ_post``.

    Returns:
        dict[tuple[int, int, int], float]: ``(KPAR, NCORE, klength) -> elapsed
        seconds``.
    """
    if not path_time.is_file():
        fail("post-time file not found: %s" % path_time)

    dict_seconds = {}
    lunknown = []
    for line in path_time.read_text(encoding="utf-8").splitlines():
        match_time = TIME_PATTERN.search(line)
        if match_time is None:
            continue
        label = match_time.group(1)
        case = parse_case_label(label, path_time.parent / "y_dir")
        if case is None:
            lunknown.append(label)
            continue
        if case in dict_seconds:
            fail("duplicate timing row for KPAR=%d, NCORE=%d, A %d in %s"
                 % (case[0], case[1], case[2], path_time))
        dict_seconds[case] = float(match_time.group(2))

    if lunknown:
        warn("ignored timing rows with unrecognized labels: %s"
             % ", ".join(sorted(set(lunknown))))
    if not dict_seconds:
        fail("no completed kpar_<KPAR>_ncore_<NCORE>[_klen_<L>] timing rows "
             "found in %s" % path_time)
    return dict_seconds


def read_kpar_ncore_energies(path_data: Path) -> dict[tuple[int, int], float]:
    """Read ``energy(sigma->0)`` values from ``y_post_data.txt``.

    Args:
        path_data (Path): File generated by ``pei_vasp_univ_post``.

    Returns:
        dict[tuple[int, int, int], float]: ``(KPAR, NCORE, klength) -> total
        energy`` in eV.
    """
    if not path_data.is_file():
        fail("post-data file not found: %s" % path_data)

    dict_energy = {}
    lunknown = []
    for line in path_data.read_text(encoding="utf-8").splitlines():
        match_energy = ENERGY_PATTERN.search(line)
        if match_energy is None:
            continue
        label = match_energy.group(1)
        case = parse_case_label(label, path_data.parent / "y_dir")
        if case is None:
            lunknown.append(label)
            continue
        if case in dict_energy:
            fail("duplicate energy row for KPAR=%d, NCORE=%d, A %d in %s"
                 % (case[0], case[1], case[2], path_data))
        dict_energy[case] = float(match_energy.group(2))

    if lunknown:
        warn("ignored energy rows with unrecognized labels: %s"
             % ", ".join(sorted(set(lunknown))))
    if not dict_energy:
        fail("no kpar_<KPAR>_ncore_<NCORE>[_klen_<L>] energy rows found in %s"
             % path_data)
    return dict_energy


def read_vasp_natoms(path_structure: Path) -> int:
    """Read the atom count from a VASP 4/5 POSCAR-style file.

    Args:
        path_structure (Path): POSCAR or CONTCAR path.

    Returns:
        int: Total number of atoms.

    Raises:
        ValueError: If the atom-count line cannot be parsed.
    """
    lline = path_structure.read_text(encoding="utf-8").splitlines()
    if len(lline) < 7:
        raise ValueError("fewer than seven lines")

    # VASP 4 places counts on line 6; VASP 5 inserts element names there and
    # moves counts to line 7. Only the count is needed, so no POTCAR lookup is
    # required for old structures whose comment does not contain element names.
    for index in (5, 6):
        ltoken = lline[index].split("!", maxsplit=1)[0].split()
        try:
            lcount = [int(token) for token in ltoken]
        except ValueError:
            continue
        if lcount and sum(lcount) > 0 and all(count >= 0 for count in lcount):
            return sum(lcount)
    raise ValueError("no positive atom-count line at POSCAR line 6 or 7")


def get_kpar_ncore_natoms(path_workflow: Path) -> int:
    """Find one reference structure and return its atom count.

    The relaxed source is preferred because it is the workflow's physical
    reference. A generated case is an equivalent fallback: every case copied
    that same structure to ``POSCAR``.

    Args:
        path_workflow (Path): KPAR/NCORE workflow directory.

    Returns:
        int: Number of atoms in each benchmark calculation.
    """
    path_relax = path_workflow.parent / "y_full_relax"
    lcandidate = [path_relax / "CONTCAR", path_relax / "POSCAR"]
    for path_case in sorted((path_workflow / "y_dir").iterdir()):
        if path_case.is_dir():
            lcandidate.extend([path_case / "POSCAR", path_case / "CONTCAR"])

    linvalid = []
    for path_structure in lcandidate:
        if not path_structure.is_file():
            continue
        try:
            natoms = read_vasp_natoms(path_structure)
        except (OSError, UnicodeError, ValueError) as exc:
            linvalid.append("%s (%s)" % (path_structure, exc))
            continue
        print("🧮 atoms        : %d from %s" % (natoms, path_structure))
        return natoms

    if linvalid:
        warn("unreadable reference structures: %s" % "; ".join(linvalid))
    fail("no usable CONTCAR/POSCAR found in %s or its y_dir cases"
         % path_relax)


def read_outcar_cost(path_outcar: Path) -> dict:
    """Scrape the cost drivers of one finished static run from OUTCAR.

    ``NKPTS`` and ``NBANDS`` are what a cost model like ``t ~ NKPTS*NBANDS^2``
    is built on, and the ``LOOP:`` lines separate "how expensive is one
    electronic step" from "how many steps did SCF need" -- two effects that a
    single elapsed time cannot tell apart.

    Args:
        path_outcar (Path): OUTCAR of a benchmark case.

    Returns:
        dict: ``nkpts``, ``nbands``, ``n_scf``, ``loop_s`` (mean electronic-step
        real time) and ``loop_plus_s`` (ionic-step real time); missing values
        are ``None``.
    """
    dict_cost = {"nkpts": None, "nbands": None, "n_scf": 0,
                 "loop_s": None, "loop_plus_s": None}
    if not path_outcar.is_file():
        return dict_cost

    lloop = []
    for line in path_outcar.read_text(encoding="utf-8", errors="replace").splitlines():
        if dict_cost["nkpts"] is None:
            match_dim = NKPTS_PATTERN.search(line)
            if match_dim is not None:
                dict_cost["nkpts"] = int(match_dim.group(1))
                dict_cost["nbands"] = int(match_dim.group(2))
                continue
        match_loop = LOOP_PATTERN.search(line)
        if match_loop is not None:
            lloop.append(float(match_loop.group(1)))
            continue
        match_plus = LOOP_PLUS_PATTERN.search(line)
        if match_plus is not None:
            dict_cost["loop_plus_s"] = float(match_plus.group(1))
    if lloop:
        dict_cost["n_scf"] = len(lloop)
        dict_cost["loop_s"] = sum(lloop) / len(lloop)
    return dict_cost


def get_kpar_ncore_costs(path_workflow: Path,
                         lcase: list = None) -> dict[tuple[int, int, int], dict]:
    """Collect the OUTCAR cost drivers of every benchmark case.

    Args:
        path_workflow (Path): KPAR/NCORE workflow directory.
        lcase (list): Cases to read; ``None`` reads every case directory found.

    Returns:
        dict[tuple[int, int, int], dict]: ``case -> read_outcar_cost`` mapping.
    """
    dict_cost = {}
    for path_case in sorted((path_workflow / "y_dir").iterdir()):
        if not path_case.is_dir():
            continue
        case = parse_case_label(path_case.name, path_case.parent)
        if case is None or (lcase is not None and case not in lcase):
            continue
        dict_cost[case] = read_outcar_cost(path_case / "OUTCAR")
    return dict_cost


def get_kpar_ncore_forces(path_workflow: Path,
                          lcase: list = None) -> dict[tuple[int, int, int], object]:
    """Read the final forces of every benchmark case.

    Args:
        path_workflow (Path): KPAR/NCORE workflow directory.
        lcase (list): Cases to read; ``None`` reads every case directory found.

    Returns:
        dict[tuple[int, int, int], np.ndarray]: ``case -> (natoms, 3)`` forces.
    """
    dict_forces = {}
    for path_case in sorted((path_workflow / "y_dir").iterdir()):
        if not path_case.is_dir():
            continue
        case = parse_case_label(path_case.name, path_case.parent)
        if case is None or (lcase is not None and case not in lcase):
            continue
        dict_forces[case] = read_outcar_forces(path_case / "OUTCAR")
    return dict_forces


def read_outcar_forces(path_outcar: Path = None):
    """Read the final-step atomic forces of one finished static run.

    Energy alone cannot say whether a cheaper setting is safe: a potential is
    fitted to forces too, and the force error is what a looser k mesh or SCF
    criterion shows first.

    Args:
        path_outcar (Path): OUTCAR of a benchmark case.

    Returns:
        np.ndarray: ``(natoms, 3)`` forces in eV/Angstrom, or ``None`` when the
        run has no readable force block yet.
    """
    if not path_outcar.is_file():
        return None
    try:
        # ase 的 OUTCAR 读取器已经处理了 VASP 各版本的排版差异，别自己数行
        from ase.io import read as ase_read
        atoms = ase_read(str(path_outcar), index=-1, format="vasp-out")
        return atoms.get_forces()
    except Exception as exc:                     # 未完成/被截断的 OUTCAR 很常见
        warn("could not read forces from %s (%s)" % (path_outcar, exc))
        return None


def get_force_stats(dict_forces: dict, case_ref: tuple = None) -> dict:
    """Per-case force magnitude and, against a reference case, its deviation.

    Args:
        dict_forces (dict): ``case -> (natoms, 3)`` force array.
        case_ref (tuple): Case whose forces are treated as the truth (e.g. the
            densest k mesh). ``None`` reports magnitudes only.

    Returns:
        dict: ``case -> {'fmax', 'frms', 'dfmax', 'dfrms'}`` in eV/Angstrom;
        the deviation entries are ``None`` without a reference or when the two
        runs hold different atom counts.
    """
    import numpy as np

    force_ref = None if case_ref is None else dict_forces.get(case_ref)
    dict_stat = {}
    for case, force in dict_forces.items():
        if force is None:
            continue
        stat = {"fmax": float(np.abs(force).max()),
                "frms": float(np.sqrt((force ** 2).sum(axis=1).mean())),
                "dfmax": None, "dfrms": None}
        if force_ref is not None and force_ref.shape == force.shape:
            ldiff = force - force_ref
            # dfmax 取「单个原子的力矢量差的模」的最大值：训练里真正要紧的是
            # 最差的那个原子，不是逐分量平均
            stat["dfmax"] = float(np.linalg.norm(ldiff, axis=1).max())
            stat["dfrms"] = float(np.sqrt((ldiff ** 2).sum(axis=1).mean()))
        dict_stat[case] = stat
    return dict_stat


def read_vasp_case(path_case: Path = None, label: str = None) -> dict:
    """Read one finished static run as a single row of a benchmark table.

    A parameter study spreads its axes over separate workflow trees (one per
    ENCUT, EDIFF, cell size ...), so the comparison needs a reader that takes
    a bare VASP directory and reports both what was asked for (ENCUT, EDIFF,
    k-mesh length, KPAR/NCORE, atom count) and what it cost and produced
    (NKPTS, SCF steps, elapsed time, energy, forces).

    Args:
        path_case (Path): A VASP run directory holding INCAR/KPOINTS/OUTCAR.
        label (str): Row name; ``None`` uses the directory name.

    Returns:
        dict: One row. ``finished`` is False when OUTCAR carries no elapsed
        time yet, in which case the result columns stay ``None``; ``forces``
        holds the ``(natoms, 3)`` array for later deviation statistics.
    """
    path_case = Path(path_case)
    row = {"label": path_case.name if label is None else label,
           "path": str(path_case), "finished": False,
           "encut": None, "ediff": None, "klength": None,
           "kpar": None, "ncore": None, "natoms": None,
           "nkpts": None, "nbands": None, "n_scf": None, "loop_s": None,
           "elapsed_s": None, "energy_eV": None, "energy_per_atom_eV": None,
           "fmax": None, "frms": None, "forces": None}

    path_incar = path_case / "INCAR"
    if path_incar.is_file():
        for line in path_incar.read_text(encoding="utf-8", errors="replace").splitlines():
            line = line.split("#")[0].split("!")[0]
            if "=" not in line:
                continue
            tag, _, value = line.partition("=")
            tag, value = tag.strip().upper(), value.strip().split()[0] if value.strip() else ""
            if tag in ("ENCUT", "EDIFF") and value:
                row[tag.lower()] = float(value)
            elif tag in ("KPAR", "NCORE") and value:
                row[tag.lower()] = int(float(value))

    path_kpoints = path_case / "KPOINTS"
    if path_kpoints.is_file():
        lline = path_kpoints.read_text(encoding="utf-8", errors="replace").splitlines()
        for index, line in enumerate(lline):
            match_k = KLENGTH_PATTERN.match(line)
            if match_k is None:
                continue
            value = match_k.group(1)
            if value is None and index + 1 < len(lline):
                value = lline[index + 1].strip().split()[0] if lline[index + 1].strip() else None
            if value is not None:
                row["klength"] = float(value)
            break

    for name in ("POSCAR", "CONTCAR"):
        if (path_case / name).is_file():
            row["natoms"] = read_vasp_natoms(path_case / name)
            break

    path_outcar = path_case / "OUTCAR"
    if not path_outcar.is_file():
        return row
    row.update(read_outcar_cost(path_outcar))
    text = path_outcar.read_text(encoding="utf-8", errors="replace")
    lelapsed = ELAPSED_PATTERN.findall(text)
    lenergy = SIGMA0_PATTERN.findall(text)
    if not lelapsed:                       # 还在跑或被 kill：只回报设置，不编结果
        return row
    row["finished"] = True
    row["elapsed_s"] = float(lelapsed[-1])
    if lenergy:
        row["energy_eV"] = float(lenergy[-1])
        if row["natoms"]:
            row["energy_per_atom_eV"] = row["energy_eV"] / row["natoms"]
    force = read_outcar_forces(path_outcar)
    if force is not None:
        import numpy as np
        row["forces"] = force
        row["fmax"] = float(np.abs(force).max())
        row["frms"] = float(np.sqrt((force ** 2).sum(axis=1).mean()))
    return row


def collect_vasp_cases(dict_case: dict = None, label_ref: str = None):
    """Collect several VASP runs into one table, optionally referenced to one.

    Args:
        dict_case (dict): ``label -> case directory``. Insertion order is kept.
        label_ref (str): Row treated as the converged answer; the table then
            also carries ``delta_energy_meV_per_atom``, ``dfmax`` and ``dfrms``
            against it. Cases with a different atom count keep ``nan`` there --
            energies and forces of different cells are not comparable one to one.

    Returns:
        pandas.DataFrame: One row per case, ``forces`` dropped, sorted as given.
    """
    import numpy as np
    import pandas as pd

    lrow = [read_vasp_case(path_case, label=label)
            for label, path_case in dict_case.items()]
    dict_force = {row["label"]: row["forces"] for row in lrow}
    dict_stat = get_force_stats({key: value for key, value in dict_force.items()
                                 if value is not None},
                                case_ref=label_ref)
    row_ref = next((row for row in lrow if row["label"] == label_ref), None)
    for row in lrow:
        stat = dict_stat.get(row["label"], {})
        row["dfmax"] = stat.get("dfmax", None)
        row["dfrms"] = stat.get("dfrms", None)
        row["delta_energy_meV_per_atom"] = None
        if (row_ref is not None and row["energy_per_atom_eV"] is not None
                and row_ref["energy_per_atom_eV"] is not None
                and row["natoms"] == row_ref["natoms"]):
            row["delta_energy_meV_per_atom"] = 1000.0 * (
                row["energy_per_atom_eV"] - row_ref["energy_per_atom_eV"])
        row.pop("forces")
        # 缺测值直接写成 nan：留 None 会让整列退化成 object 列
        for key, value in row.items():
            if value is None:
                row[key] = np.nan
    return pd.DataFrame(lrow)


def get_delta_energies(
        dict_energy: dict[tuple[int, int, int], float],
        natoms: int,
        energy0: float = None) -> dict[tuple[int, int, int], float]:
    """Convert total energies to ``(E - E0)`` in meV/atom.

    Args:
        dict_energy (dict[tuple[int, int, int], float]): Total energies in eV.
        natoms (int): Number of atoms in each benchmark case.
        energy0 (float): Reference energy in eV; ``None`` uses the minimum
            measured one. The k-mesh scan passes the densest mesh instead, so
            its panel reads as a convergence curve rather than a spread.

    Returns:
        dict[tuple[int, int, int], float]: Relative energies in meV/atom.
    """
    if not dict_energy:
        fail("at least one completed energy is required")
    if natoms <= 0:
        fail("natoms must be positive")
    energy0 = min(dict_energy.values()) if energy0 is None else energy0
    return {
        case: (energy - energy0) * 1000.0 / natoms
        for case, energy in dict_energy.items()
    }


def get_missing_pairs(dict_seconds: dict[tuple[int, int, int], float],
                      klength: int = None) -> list[tuple[int, int]]:
    """Return default benchmark pairs that have no elapsed time yet.

    Args:
        dict_seconds (dict[tuple[int, int, int], float]): Parsed elapsed seconds.
        klength (int): Only look at this k-mesh length; ``None`` uses the single
            length present (and reports nothing when several are mixed, where
            "the default matrix" is no longer a meaningful expectation).
    """
    lklength = sorted({case[2] for case in dict_seconds})
    if klength is None:
        if len(lklength) != 1:
            return []
        klength = lklength[0]
    lmissing = []
    for kpar in KPAR_ORDER:
        for ncore in EXPECTED_PAIRS[kpar]:
            if (kpar, ncore, klength) not in dict_seconds:
                lmissing.append((kpar, ncore))
    return lmissing


def write_kpar_ncore_times(
        dict_seconds: dict[tuple[int, int, int], float],
        dict_energy: dict[tuple[int, int, int], float],
        dict_delta_energy: dict[tuple[int, int, int], float],
        natoms: int,
        path_save: Path,
        dict_cost: dict[tuple[int, int, int], dict] = None,
        dict_force: dict[tuple[int, int, int], dict] = None) -> None:
    """Write timing/energy data and the measured fastest combination.

    Args:
        dict_seconds (dict[tuple[int, int, int], float]): Parsed elapsed seconds.
        dict_energy (dict[tuple[int, int, int], float]): Completed total energies in eV.
        dict_delta_energy (dict[tuple[int, int, int], float]): Relative energies in
            meV/atom.
        natoms (int): Number of atoms in each benchmark case.
        path_save (Path): Output text path.
        dict_cost (dict[tuple[int, int, int], dict]): OUTCAR cost drivers per
            case; ``None`` writes the cost columns as ``nan``.
        dict_force (dict[tuple[int, int, int], dict]): Force statistics per case
            from :func:`get_force_stats`; ``None`` writes them as ``nan``.
    """
    dict_cost = {} if dict_cost is None else dict_cost
    dict_force = {} if dict_force is None else dict_force
    best_case = min(dict_seconds, key=dict_seconds.get)
    energy0 = min(dict_energy.values())
    lline = [
        "# KPAR/NCORE/k-mesh VASP timing/energy benchmark",
        "# best: KPAR=%d NCORE=%d A %d time_min=%.6f"
        % (best_case[0], best_case[1], best_case[2],
           dict_seconds[best_case] / 60.0),
        "# natoms=%d E0_eV=%.10f" % (natoms, energy0),
        "# klen      automatic k-mesh length (KPOINTS mode A)",
        "# nkpts/nbands/n_scf/loop_s  cost drivers scraped from OUTCAR",
        "# fmax/frms   force magnitude, eV/A; dfmax/dfrms  deviation from the "
        "reference case (densest k mesh), eV/A",
        "# KPAR NCORE klen elapsed_seconds time_min energy_eV "
        "delta_energy_meV_per_atom nkpts nbands n_scf loop_s fmax frms dfmax dfrms",
    ]
    for klength in sorted({case[2] for case in dict_seconds}):
        for kpar in KPAR_ORDER:
            lcase = sorted(
                [case for case in dict_seconds
                 if case[0] == kpar and case[2] == klength],
                key=lambda case: case[1],
            )
            for case in lcase:
                seconds = dict_seconds[case]
                energy = dict_energy.get(case, float("nan"))
                delta_energy = dict_delta_energy.get(case, float("nan"))
                cost = dict_cost.get(case, {})
                force = dict_force.get(case, {})
                lline.append(
                    "%6d %6d %6d %15.6f %12.6f %18.10f %25.10f %8s %8s %6s %12s "
                    "%10s %10s %10s %10s"
                    % (case[0], case[1], case[2], seconds, seconds / 60.0,
                       energy, delta_energy,
                       cost.get("nkpts", "nan"), cost.get("nbands", "nan"),
                       cost.get("n_scf", "nan"),
                       "nan" if cost.get("loop_s") is None
                       else "%.3f" % cost["loop_s"],
                       *["nan" if force.get(key) is None else "%.6f" % force[key]
                         for key in ("fmax", "frms", "dfmax", "dfrms")]))
    path_save.write_text("\n".join(lline) + "\n", encoding="utf-8")


def post_kpar_ncore(path_workflow: str = None, run_post: bool = True,
                     save_fig_path: str = None,
                     save_txt_path: str = None) -> dict[tuple[int, int, int], float]:
    """Run standard post-processing and create timing/energy outputs.

    With one k-mesh length the outputs are the historical pair: one table and
    one KPAR/NCORE figure. With several, each length also gets its own figure
    (``..._klen<L>.pdf``) and a k-mesh scan figure (``..._klen.pdf``) showing
    time and energy convergence against the densest mesh.

    Args:
        path_workflow (str): Directory containing ``y_dir``.
        run_post (bool): Run ``pei_vasp_univ_post`` before reading times.
        save_fig_path (str): Figure path. Defaults to
            ``p_post_kpar_ncore.pdf`` inside the workflow.
        save_txt_path (str): Table path. Defaults to
            ``p_post_kpar_ncore.txt`` inside the workflow.

    Returns:
        dict[tuple[int, int, int], float]: Parsed elapsed seconds.
    """
    ### check ------------------------------------------------------------------
    if path_workflow is None:
        fail("path_workflow is required")
    path_workflow = Path(path_workflow).resolve()
    if not path_workflow.is_dir():
        fail("workflow directory not found: %s" % path_workflow)
    if not (path_workflow / "y_dir").is_dir():
        fail("y_dir not found under %s" % path_workflow)
    # to here: the standard workflow layout is available

    ### prepare ----------------------------------------------------------------
    if run_post:
        run_univ_post(path_workflow)
    path_time = path_workflow / "y_post_time.txt"
    path_data = path_workflow / "y_post_data.txt"
    dict_seconds = read_kpar_ncore_times(path_time)
    dict_energy_all = read_kpar_ncore_energies(path_data)
    lmissing_energy = sorted(set(dict_seconds) - set(dict_energy_all))
    lextra_energy = sorted(set(dict_energy_all) - set(dict_seconds))
    if lmissing_energy:
        warn("completed timing rows without energy: %s"
             % ", ".join("(%d,%d,A%d)" % case for case in lmissing_energy))
    if lextra_energy:
        warn("ignored energies without completed timing rows: %s"
             % ", ".join("(%d,%d,A%d)" % case for case in lextra_energy))
    dict_energy = {
        case: dict_energy_all[case]
        for case in dict_seconds if case in dict_energy_all
    }
    if not dict_energy:
        fail("no completed benchmark has both elapsed time and energy")
    natoms = get_kpar_ncore_natoms(path_workflow)
    dict_cost = get_kpar_ncore_costs(path_workflow, lcase=list(dict_seconds))
    dict_forces = get_kpar_ncore_forces(path_workflow, lcase=list(dict_seconds))
    lklength = sorted({case[2] for case in dict_seconds})
    # A single mesh keeps the historical reference (the lowest measured energy);
    # a scan is a convergence question, so the densest mesh is the reference.
    energy0 = (None if len(lklength) == 1 else
               min([dict_energy[case] for case in dict_energy
                    if case[2] == max(lklength)]))
    dict_delta_energy = get_delta_energies(dict_energy, natoms, energy0=energy0)
    # 力的参照与能量同一个：最密网格那一档里最快的算例（并行设置不改物理）
    lcase_ref = [case for case in dict_seconds if case[2] == max(lklength)]
    case_ref = min(lcase_ref, key=dict_seconds.get) if lcase_ref else None
    dict_force = get_force_stats(dict_forces, case_ref=case_ref)
    path_fig = (Path(save_fig_path).resolve() if save_fig_path
                else path_workflow / "p_post_kpar_ncore.pdf")
    path_txt = (Path(save_txt_path).resolve() if save_txt_path
                else path_workflow / "p_post_kpar_ncore.txt")
    # to here: matched timing/energy data and output paths are ready

    ### main -------------------------------------------------------------------
    write_kpar_ncore_times(
        dict_seconds, dict_energy, dict_delta_energy, natoms, path_txt,
        dict_cost=dict_cost, dict_force=dict_force)
    dict_minutes = {
        case: seconds / 60.0 for case, seconds in dict_seconds.items()
    }
    lxtick = sorted({
        ncore for lncore in EXPECTED_PAIRS.values() for ncore in lncore
    })
    lpath_fig = []
    for klength in lklength:
        # one KPAR/NCORE figure per mesh; the single-mesh case keeps the plain
        # file name so existing trees and pei_vasp_plot_all are unaffected
        path_fig_k = (path_fig if len(lklength) == 1 else
                      path_fig.with_name("%s_klen%d%s"
                                         % (path_fig.stem, klength, path_fig.suffix)))
        my_plot_kpar_ncore(
            dict_time={(case[0], case[1]): value
                       for case, value in dict_minutes.items() if case[2] == klength},
            dict_delta_energy={(case[0], case[1]): value
                               for case, value in dict_delta_energy.items()
                               if case[2] == klength},
            lkpar=KPAR_ORDER,
            lncore=lxtick,
            if_save=True,
            savefile=path_fig_k,
            if_close=True,
        )
        lpath_fig.append(path_fig_k)
    path_fig_scan = None
    if len(lklength) > 1:
        path_fig_scan = path_fig.with_name("%s_klen%s"
                                           % (path_fig.stem, path_fig.suffix))
        my_plot_kpar_ncore_klen(
            dict_time=dict_minutes,
            dict_delta_energy=dict_delta_energy,
            dict_cost=dict_cost,
            dict_force=dict_force,
            lklength=lklength,
            if_save=True,
            savefile=path_fig_scan,
            if_close=True,
        )
    lmissing = get_missing_pairs(dict_seconds)
    best_case = min(dict_seconds, key=dict_seconds.get)
    print("\n================ 📊 Summary")
    print("completed=%d    missing=%d    k-mesh=%s"
          % (len(dict_seconds), len(lmissing),
             ", ".join("A %d" % klength for klength in lklength)))
    if lmissing:
        warn("missing elapsed times: %s"
             % ", ".join("(%d,%d)" % pair for pair in lmissing))
    print("✅ best measured: KPAR=%d NCORE=%d A %d time=%.6f min"
          % (best_case[0], best_case[1], best_case[2],
             dict_seconds[best_case] / 60.0))
    print("✅ energy E0   : %.10f eV; max ΔE=%.6f meV/atom"
          % (min(dict_energy.values()), max(dict_delta_energy.values())))
    ldfmax = [stat["dfmax"] for stat in dict_force.values() if stat["dfmax"] is not None]
    if ldfmax:
        print("✅ force       : max deviation %.6f eV/A from KPAR=%d NCORE=%d A %d"
              % (max(ldfmax), case_ref[0], case_ref[1], case_ref[2]))
    if len(lklength) > 1:
        # the number the campaign actually needs: cost and accuracy of the
        # cheaper meshes relative to the densest one
        for klength in lklength:
            lcase_k = [case for case in dict_seconds if case[2] == klength]
            best_k = min(lcase_k, key=dict_seconds.get)
            print("   A %-4d: best %8.2f min (KPAR=%d NCORE=%d), "
                  "NKPTS=%s, ΔE=%+.3f meV/atom, dF=%s eV/A"
                  % (klength, dict_seconds[best_k] / 60.0, best_k[0], best_k[1],
                     dict_cost.get(best_k, {}).get("nkpts", "?"),
                     dict_delta_energy.get(best_k, float("nan")),
                     "n/a" if dict_force.get(best_k, {}).get("dfmax") is None
                     else "%.4f" % dict_force[best_k]["dfmax"]))
    print("📄 table : %s" % path_txt)
    for path_one in lpath_fig:
        print("📈 figure: %s" % path_one)
    if path_fig_scan is not None:
        print("📈 figure: %s" % path_fig_scan)
    print("🎉 Done: KPAR/NCORE timing/energy post-processing complete")
    return dict_seconds
