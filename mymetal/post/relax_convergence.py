"""
relax_convergence post-processing submodule

Reads the per-ionic-step convergence data files produced by the bash helper
``pei_vasp_univ_extract_convergence`` and renders one plot per y_dir
sub-directory showing the ionic-relaxation trajectory:

    * energy(sigma->0) relative to the last frame, in meV/atom
    * max force norm, in eV/Ang, together with the |EDIFFG| force criterion

Data files live in ``y_post_convergence/y_post_convergence_<name>.txt`` and the
plots are written next to them as ``y_post_convergence_<name>.pdf``. This
mirrors the layout of the other ``yin_vasp_plot_*`` post-processing scripts
(a ``y_post_<x>/`` folder holding ``y_post_<x>_<jobn>`` outputs).

Note: this is distinct from :mod:`mymetal.post.convergence`, which handles
ENCUT / KPOINTS convergence *testing*; here we track the convergence of a
single ionic relaxation.

Functions:
    - my_univ_post_convergence: driver, plot every data file in the folder.
    - read_univ_post_convergence: parse one bash-generated data file.

Change log:
    - Written by J. P. on 2026.07.05.
"""

from pathlib import Path

from mymetal.universal.plot.workflow import my_plot_relax_convergence


def my_univ_post_convergence(workflow: str = 'y_post_convergence',
                             yscale: str = 'log') -> None:
    """Plot every convergence data file found in ``workflow``.

    Run from the directory that contains ``y_dir`` (and thus ``workflow``).
    For each ``<workflow>/<workflow>_<name>.txt`` produced by the bash helper
    ``pei_vasp_univ_extract_convergence`` a ``<...>.pdf`` is written alongside.

    Both panels plot absolute values so a large early force/energy does not
    hide the near-converged detail on the log axis.

    Args:
        workflow (str): output folder holding the data files (default
            ``'y_post_convergence'``).
        yscale (str): y-axis scale for both panels, ``'log'`` (default) or
            ``'linear'``.
    """
    outdir = Path(workflow)
    if not outdir.is_dir():
        print(f"==> no {outdir}/ folder, run pei_vasp_univ_extract_convergence first.")
        return

    datafiles = sorted(outdir.glob(f"{workflow}_*.txt"))
    if not datafiles:
        print(f"==> no {workflow}_*.txt data files in {outdir}/.")
        return

    for datafile in datafiles:
        lframe, lenergy, lforce, natoms, ediffg, isif = read_univ_post_convergence(datafile)
        if len(lframe) == 0:
            print(f"    ⏭️  {datafile.name}: no ionic steps, skip.")
            continue

        output_file = datafile.with_suffix('.pdf')
        my_plot_relax_convergence(
            lframe=lframe,
            lenergy=lenergy,
            lforce=lforce,
            natoms=natoms,
            ediffg=ediffg,
            isif=isif,
            yscale=yscale,
            if_save=True,
            savefile=output_file,
            if_close=True,
        )
        print(f"    ✅ {output_file}")


def read_univ_post_convergence(path) -> tuple:
    """Parse a convergence data file written by ``pei_vasp_univ_extract_convergence``.

    Args:
        path: path to a ``y_post_convergence_<name>.txt`` file.

    Returns:
        tuple: (lframe, lenergy, lforce, natoms, ediffg, isif)
            - lframe  (list[int]):   ionic-step indices, starting at 1.
            - lenergy (list[float]): energy(sigma->0) per step, eV.
            - lforce  (list[float]): max force norm per step, eV/Ang.
            - natoms  (int | None):  number of atoms.
            - ediffg  (float | None): EDIFFG (eV/Ang when < 0).
            - isif    (int | None):  ISIF tag.
    """
    natoms, ediffg, isif = None, None, None
    lframe, lenergy, lforce = [], [], []

    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            # header metadata: "# natoms 3", "# ediffg -.1E-02", "# isif 4"
            if parts[0] == '#':
                if len(parts) >= 3 and parts[1] == 'natoms':
                    natoms = _to_int(parts[2])
                elif len(parts) >= 3 and parts[1] == 'ediffg':
                    ediffg = _to_float(parts[2])
                elif len(parts) >= 3 and parts[1] == 'isif':
                    isif = _to_int(parts[2])
                continue
            # data rows: frame energy force (column-header line is skipped by the try)
            try:
                frame = int(parts[0])
                energy = float(parts[1])
                force = float(parts[2])
            except (ValueError, IndexError):
                continue
            lframe.append(frame)
            lenergy.append(energy)
            lforce.append(force)

    return lframe, lenergy, lforce, natoms, ediffg, isif


def _to_float(s):
    try:
        return float(s)
    except (ValueError, TypeError):
        return None


def _to_int(s):
    try:
        return int(float(s))
    except (ValueError, TypeError):
        return None
