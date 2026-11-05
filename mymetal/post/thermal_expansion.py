"""
thermal_expansion post module

Averages of a finished ``y_thermal_expansion`` tree (built by ``pei_vasp_run_thermal_expansion``,
cases rendered from ``lmp_utils/template/npt_case_template.in``) into one table, one row per
temperature plus the 0 K row of the cell relaxation every case starts with.

Per case it reads:
    POSCAR               comment line ``T:<T> size:<nx>x<ny>x<nz> seed:<seed>``, atom count
    y_lmp_relax0.txt     0 K box and energy
    y_lmp_md_prod.txt    fix ave/time block rows of the production run
    log.lammps           ``### NNP EW SUMMARY ###`` lines of the production run (hdnnp only)

Lattice lengths are those of the unit cell: the averaged box vectors divided by the supercell
multipliers, |a1| = lx/nx, |a2| = |(xy, ly)|/ny, |a3| = |(xz, yz, lz)|/nz. Error bars are
standard errors of the mean over :data:`NBLOCK` equal consecutive blocks of the production
rows, so correlations shorter than a block do not shrink them.

Functions:
    - read_case_meta: Temperature, supercell size and atom count of a case.
    - read_lmp_ave_time: Rows of a fix ave/time file as a DataFrame.
    - get_lattice_lengths: Unit-cell |a1| |a2| |a3| from box columns.
    - get_block_sem: Standard error of the mean from consecutive blocks.
    - count_lmp_ew_prod: Extrapolation warnings of the last run in a log.lammps.
    - post_thermal_expansion_case: One row of the table.
    - get_relax0_row: The 0 K row.
    - post_thermal_expansion: Whole tree -> p_post_thermal_expansion.txt (+ figure).
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd

from mymetal.universal.plot.workflow import my_plot_thermal_expansion
from mymetal.universal.print.print import fail, warn

DONE_MARKER = 'y_lmp_done.txt'
RELAX0_FILE = 'y_lmp_relax0.txt'
PROD_FILE = 'y_lmp_md_prod.txt'
LOG_FILE = 'log.lammps'
LCOL_MD = ['step', 'T', 'P', 'lx', 'ly', 'lz', 'xy', 'xz', 'yz', 'vol', 'pe', 'msd',
           'f_fcc', 'f_hcp']
LCOL_RELAX0 = ['lx', 'ly', 'lz', 'xy', 'xz', 'yz', 'vol', 'pe']
NBLOCK = 10
# a solid's MSD from its 0 K sites saturates at a few 0.1 A^2 even near melting; a liquid
# diffuses tens of A^2 in a 100 ps production run
MSD_MELT = 1.0
TXT_FILE = 'p_post_thermal_expansion.txt'
FIG_FILE = 'p_post_thermal_expansion.pdf'
RE_EW = re.compile(r'### NNP EW SUMMARY ### TS:\s*(\d+)\s+EW\s+(\d+)')
RE_RUN = re.compile(r'^run\s', re.M)


def read_case_meta(path_case: Path = None) -> dict:
    """Temperature, supercell size and atom count from the case POSCAR.

    Args:
        path_case: ``y_dir/T-<T>``.

    Returns:
        dict: ``{'T': int, 'size': [nx, ny, nz], 'natoms': int}``.
    """
    llines = (path_case / 'POSCAR').read_text().splitlines()
    dict_tag = dict(token.split(':', 1) for token in llines[0].split() if ':' in token)
    if 'T' not in dict_tag or 'size' not in dict_tag:
        fail('POSCAR comment of %s lacks T:/size: (%r)' % (path_case, llines[0]))
    natoms = sum(int(n) for n in llines[6].split())
    return {'T': int(dict_tag['T']), 'size': [int(n) for n in dict_tag['size'].split('x')],
            'natoms': natoms}


def read_lmp_ave_time(path_file: Path = None, lcol: list = None) -> pd.DataFrame:
    """Rows of a ``fix ave/time`` file (``#`` header lines skipped).

    Args:
        path_file: The file.
        lcol: Column names, step first.

    Returns:
        pd.DataFrame: One row per output step.
    """
    data = np.loadtxt(path_file, comments='#', ndmin=2)
    if data.shape[1] != len(lcol):
        fail('%s has %d columns, expected %d (%s)'
             % (path_file, data.shape[1], len(lcol), ' '.join(lcol)))
    return pd.DataFrame(data, columns=lcol)


def get_lattice_lengths(df: pd.DataFrame = None, lsize: list = None) -> np.ndarray:
    """Unit-cell lattice lengths per row.

    Args:
        df: Columns lx ly lz xy xz yz (one row per sample or block).
        lsize: Supercell multipliers ``[nx, ny, nz]``.

    Returns:
        np.ndarray: ``(nrow, 3)`` of |a1| |a2| |a3|, Angstrom.
    """
    la1 = df['lx'].to_numpy()
    la2 = np.hypot(df['xy'].to_numpy(), df['ly'].to_numpy())
    la3 = np.sqrt(df['xz'].to_numpy() ** 2 + df['yz'].to_numpy() ** 2 + df['lz'].to_numpy() ** 2)
    return np.column_stack([la1, la2, la3]) / np.array(lsize, dtype=float)


def get_block_sem(values: np.ndarray = None, nblock: int = NBLOCK) -> float:
    """Standard error of the mean over ``nblock`` equal consecutive blocks.

    Args:
        values: Time series (rows beyond a whole number of blocks are dropped).
        nblock: Number of blocks.

    Returns:
        float: SEM, nan with fewer rows than blocks.
    """
    values = np.asarray(values, dtype=float)
    nper = len(values) // nblock
    if nper == 0:
        return float('nan')
    lmean = values[:nper * nblock].reshape(nblock, nper).mean(axis=1)
    return float(lmean.std(ddof=1) / np.sqrt(nblock))


def count_lmp_ew_prod(path_log: Path = None) -> float:
    """Sum of the hdnnp extrapolation warnings logged after the last ``run`` command.

    Args:
        path_log: ``log.lammps`` (LAMMPS echoes every input line into it).

    Returns:
        float: Warning count of the production run; nan if the log is missing or the pair
        style printed no EW summary (non-hdnnp, or ``showewsum 0``).
    """
    if not path_log.is_file():
        return float('nan')
    text = path_log.read_text(errors='replace')
    lmatch = list(RE_RUN.finditer(text))
    if not lmatch:
        return float('nan')
    lew = [int(ew) for _, ew in RE_EW.findall(text[lmatch[-1].start():])]
    return float(sum(lew)) if lew else float('nan')


def post_thermal_expansion_case(path_case: Path = None) -> dict:
    """Production averages of one finished case.

    Args:
        path_case: ``y_dir/T-<T>``.

    Returns:
        dict: One row of the table (see :func:`post_thermal_expansion`).
    """
    dict_meta = read_case_meta(path_case)
    df = read_lmp_ave_time(path_case / PROD_FILE, LCOL_MD)
    natoms = dict_meta['natoms']
    lattice = get_lattice_lengths(df, dict_meta['size'])
    vol_atom = df['vol'].to_numpy() / natoms
    pe_atom = df['pe'].to_numpy() / natoms
    dict_row = {'T_set': dict_meta['T'], 'T': df['T'].mean(), 'P_bar': df['P'].mean(),
                'nrow': len(df)}
    for i, name in enumerate(['a1', 'a2', 'a3']):
        dict_row[name] = lattice[:, i].mean()
        dict_row[name + '_sem'] = get_block_sem(lattice[:, i])
    ratio = lattice[:, 2] / lattice[:, 0]
    dict_row.update({'a3_a1': ratio.mean(), 'a3_a1_sem': get_block_sem(ratio),
                     'V_atom': vol_atom.mean(), 'V_atom_sem': get_block_sem(vol_atom),
                     'E_atom': pe_atom.mean(), 'E_atom_sem': get_block_sem(pe_atom),
                     'msd_end': df['msd'].iloc[-1], 'f_fcc': df['f_fcc'].mean(),
                     'f_hcp': df['f_hcp'].mean(),
                     'ew_prod': count_lmp_ew_prod(path_case / LOG_FILE)})
    return dict_row


def get_relax0_row(path_case: Path = None) -> dict:
    """The 0 K row from a case's cell relaxation (identical in every case of a tree).

    Args:
        path_case: Any finished ``y_dir/T-<T>``.

    Returns:
        dict: Row with T_set = T = 0 and zero error bars.
    """
    dict_meta = read_case_meta(path_case)
    data = np.loadtxt(path_case / RELAX0_FILE, ndmin=2)
    df = pd.DataFrame(data, columns=LCOL_RELAX0)
    lattice = get_lattice_lengths(df, dict_meta['size'])[0]
    natoms = dict_meta['natoms']
    return {'T_set': 0, 'T': 0.0, 'P_bar': float('nan'), 'nrow': 0,
            'a1': lattice[0], 'a1_sem': 0.0, 'a2': lattice[1], 'a2_sem': 0.0,
            'a3': lattice[2], 'a3_sem': 0.0, 'a3_a1': lattice[2] / lattice[0],
            'a3_a1_sem': 0.0, 'V_atom': df['vol'].iloc[0] / natoms, 'V_atom_sem': 0.0,
            'E_atom': df['pe'].iloc[0] / natoms, 'E_atom_sem': 0.0, 'msd_end': 0.0,
            'f_fcc': float('nan'), 'f_hcp': float('nan'), 'ew_prod': float('nan')}


def post_thermal_expansion(path_work: Path = None, if_plot: bool = True) -> tuple:
    """Table (and figure) of a ``y_thermal_expansion`` tree.

    Args:
        path_work: The tree (contains ``y_dir``).
        if_plot: Also write :data:`FIG_FILE` via ``my_plot_thermal_expansion``.

    Returns:
        tuple: ``(df, lunfinished)``; ``df`` has one row per finished case plus the 0 K row,
        sorted by T_set; it is also written to :data:`TXT_FILE` in ``path_work``.
    """
    path_ydir = Path(path_work) / 'y_dir'
    if not path_ydir.is_dir():
        fail('no y_dir in %s' % path_work)
    lpath_case = sorted(p for p in path_ydir.iterdir() if p.is_dir())
    ldone = [p for p in lpath_case if (p / DONE_MARKER).is_file()]
    lunfinished = [p.name for p in lpath_case if p not in ldone]
    if not ldone:
        fail('no finished case (%s) in %s' % (DONE_MARKER, path_ydir))
    lrow = [get_relax0_row(ldone[0])]
    for path_case in ldone:
        print('📁 %s' % path_case.name)
        lrow.append(post_thermal_expansion_case(path_case))
    df = pd.DataFrame(lrow).sort_values('T_set').reset_index(drop=True)
    path_txt = Path(path_work) / TXT_FILE
    with open(path_txt, 'w') as file:
        file.write('# thermal expansion, production averages of every finished case; T_set = 0 is '
                   'the 0 K cell relaxation\n'
                   '# a1 a2 a3: unit-cell lattice lengths (A); V_atom (A^3); E_atom: potential '
                   'energy (eV/atom); P_bar (bar)\n'
                   '# *_sem: standard error over %d production blocks; msd_end: MSD from the 0 K '
                   'sites at the end (A^2, > %g suggests melting)\n'
                   '# f_fcc f_hcp: CNA fractions; ew_prod: hdnnp extrapolation warnings in the '
                   'production run (nan: not reported)\n' % (NBLOCK, MSD_MELT))
        file.write(df.to_string(index=False, float_format=lambda x: '%.8g' % x) + '\n')
    print('📄 %s' % path_txt)
    if if_plot:
        my_plot_thermal_expansion(df, if_save=True, savefile=str(Path(path_work) / FIG_FILE))
        print('📄 %s' % (Path(path_work) / FIG_FILE))
    lmelt = df.loc[df['msd_end'] > MSD_MELT, 'T_set'].tolist()
    if lmelt:
        warn('MSD > %g A^2 (melted or transformed?) at T = %s K' % (MSD_MELT, lmelt))
    return df, lunfinished
