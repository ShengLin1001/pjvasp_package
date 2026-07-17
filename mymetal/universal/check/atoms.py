"""
atoms submodule

This submodule provides functions for analyzing atomic structures using the Common Neighbor Analysis (CNA) method. It includes
functions for counting different crystal structures (FCC, HCP, BCC, ICO, OTHER) from OVITO files, and for comparing the CNA
phase populations of two structures to detect a phase transition. These functions are designed to streamline the analysis of
atomic configurations in materials science simulations.

Functions:
    - get_cna_count: Get the count of different crystal structures (FCC, HCP, BCC, ICO, OTHER) from an OVITO file.
    - get_cna_count_dict: Get the CNA counts of an Atoms object as a phase-keyed dictionary.
    - get_cna_count_vasp: Get the CNA counts of a VASP structure file as a phase-keyed dictionary.
    - check_phase_transition: Compare the CNA phase populations of an initial and a final VASP structure.
"""

from ovito.modifiers import CommonNeighborAnalysisModifier
from pathlib import Path
from ovito.io.ase import ase_to_ovito
from ase import Atoms
from ase.io.vasp import read_vasp
from ovito.pipeline import StaticSource, Pipeline

# CNA phases in the order OVITO reports them; used as the dict key order everywhere below,
# so callers can zip/print them without re-deriving the ordering.
LPHASE = ("FCC", "HCP", "BCC", "ICO", "OTHER")

# check_phase_transition() status values. Kept as module constants rather than bare strings
# so downstream CLIs map them to exit codes by name instead of re-spelling the literals.
STATUS_SAME = "same"
STATUS_CHANGED = "changed"
STATUS_SKIPPED = "skipped"
STATUS_FAILED = "failed"


def get_cna_count(atoms: Atoms = None ) -> tuple[int, int, int, int, int]:
        """Get the count of different crystal structures (FCC, HCP, BCC, ICO, OTHER) from an OVITO file.

        Args:
            atoms (Atoms): The ASE Atoms object containing atomic data.

        Returns:
            tuple[int, int, int, int, int]: A tuple containing the counts of FCC, HCP, BCC, ICO, and OTHER structures, respectively.
        """

        if atoms is None:
            raise ValueError("The 'atoms' argument cannot be None.")

        # Convert ASE Atoms to OVITO format
        atoms_ovito = ase_to_ovito(atoms)

        # Read in OVITO for CNA analysis
        pipeline = Pipeline(source=StaticSource(data=atoms_ovito))
        cna_modifier = CommonNeighborAnalysisModifier(mode=CommonNeighborAnalysisModifier.Mode.AdaptiveCutoff)
        pipeline.modifiers.append(cna_modifier)
        results = pipeline.compute()

        fcc_count = results.attributes['CommonNeighborAnalysis.counts.FCC']
        hcp_count = results.attributes['CommonNeighborAnalysis.counts.HCP']
        bcc_count = results.attributes['CommonNeighborAnalysis.counts.BCC']
        ico_count = results.attributes['CommonNeighborAnalysis.counts.ICO']
        other_count = results.attributes['CommonNeighborAnalysis.counts.OTHER']

        return fcc_count, hcp_count, bcc_count, ico_count, other_count


def get_cna_count_dict(atoms: Atoms = None) -> dict[str, int]:
    """Get the CNA counts of an Atoms object as a phase-keyed dictionary.

    Args:
        atoms (Atoms): The ASE Atoms object containing atomic data.

    Returns:
        dict[str, int]: Counts keyed by phase name, plus ``atoms_num`` for the total.
    """
    lcount = get_cna_count(atoms)
    dict_count = {str_phase: int(int_count) for str_phase, int_count in zip(LPHASE, lcount)}
    dict_count["atoms_num"] = len(atoms)
    return dict_count


def get_cna_count_vasp(path_structure: Path = None) -> dict[str, int]:
    """Get the CNA counts of a VASP structure file as a phase-keyed dictionary.

    Reading through ASE rather than handing the file to OVITO directly keeps the
    element names that a CONTCAR carries on its species line, which the adaptive
    cutoff needs to classify neighbours sensibly.

    Args:
        path_structure (Path): POSCAR/CONTCAR-like file path.

    Returns:
        dict[str, int]: Counts keyed by phase name, plus ``atoms_num`` for the total.

    Raises:
        FileNotFoundError: When the structure is missing or empty.
    """
    path_structure = Path(path_structure)
    if not path_structure.is_file() or path_structure.stat().st_size == 0:
        raise FileNotFoundError("structure missing or empty: " + str(path_structure))

    atoms = read_vasp(str(path_structure))
    return get_cna_count_dict(atoms)


def check_phase_transition(path_dir: Path = None,
                           str_initial: str = "POSCAR",
                           str_final: str = "CONTCAR") -> dict:
    """Compare the CNA phase populations of an initial and a final VASP structure.

    A relaxation that keeps every phase population identical did not change phase;
    any per-phase count that moved means atoms were reclassified, i.e. a transition.

    The check is skipped whenever either structure is missing or empty. That is the
    intended behaviour, not an error: it is what makes the check apply only to real
    calculation directories (``parallel`` mode jobs, or the per-subdirectory jobs of
    ``each-subdir``) and silently pass over an orchestrating parent job, whose
    directory holds no structures at all.

    Args:
        path_dir (Path): Directory holding the two structures.
        str_initial (str): Reference structure file name. Default: POSCAR.
        str_final (str): Current structure file name. Default: CONTCAR.

    Returns:
        dict: Result with keys ``status`` (one of STATUS_SAME / STATUS_CHANGED /
            STATUS_SKIPPED / STATUS_FAILED), ``reason``, ``dict_initial``,
            ``dict_final``, ``ldiff`` (phase keys that differ), and the two paths
            as ``path_initial`` / ``path_final``.
    """
    if path_dir is None:
        raise ValueError("The 'path_dir' argument cannot be None.")

    path_dir = Path(path_dir).resolve()
    path_initial = path_dir / str_initial
    path_final = path_dir / str_final

    dict_result = {
        "status": STATUS_SKIPPED,
        "reason": "",
        "dict_initial": None,
        "dict_final": None,
        "ldiff": [],
        "path_initial": path_initial,
        "path_final": path_final,
    }

    # Both must be present and non-empty. VASP truncates CONTCAR to 0 bytes while it
    # rewrites it after an ionic step, so an empty file is a normal transient race,
    # not a broken run -- skip this round and let the next one catch it.
    lmissing = [
        path_item.name
        for path_item in (path_initial, path_final)
        if not path_item.is_file() or path_item.stat().st_size == 0
    ]
    if lmissing:
        dict_result["reason"] = "missing or empty: " + ", ".join(lmissing)
        return dict_result

    try:
        dict_result["dict_initial"] = get_cna_count_vasp(path_initial)
        dict_result["dict_final"] = get_cna_count_vasp(path_final)
    except Exception as exc:
        dict_result["status"] = STATUS_FAILED
        dict_result["reason"] = "CNA failed: " + repr(exc)
        return dict_result

    ldiff = [
        str_key
        for str_key in ("atoms_num",) + LPHASE
        if dict_result["dict_initial"][str_key] != dict_result["dict_final"][str_key]
    ]
    dict_result["ldiff"] = ldiff
    dict_result["status"] = STATUS_CHANGED if ldiff else STATUS_SAME
    dict_result["reason"] = ("counts differ: " + ", ".join(ldiff)) if ldiff else "all phase counts identical"
    return dict_result
