"""
lammps_input workflow module

Render the per-case LAMMPS inputs that stand in for INCAR / KPOINTS / POTCAR when a
``y_full_relax``-style workflow runs with ``-calc lammps``. The rendered ``lmp.in`` is a calc
file like INCAR: put it in ``y_full_relax/`` (or ``y_full_relax_neb/ini/`` as ``lmp.in.neb``)
and the workflow copies it into every case.

Templates live in ``lmp_utils/template/``. ``general_init.mod`` and ``general_potential.mod``
are inlined verbatim, so their settings stay single-sourced and no ``*.mod`` file has to be
copied into tens of thousands of case directories.

Functions:
    - get_lmp_template_path: Locate a file in ``lmp_utils/template``.
    - get_lmp_potential_block: general_potential.mod with pair and mass lines filled in.
    - render_lmp_template: Fill block (whole-line) and inline placeholders of a template.
    - generate_lmp_relax_input: Write a single-structure relaxation ``lmp.in``.
    - generate_lmp_neb_input: Write a climbing-image NEB ``lmp.in``.
"""

from pathlib import Path

from mymetal.universal.print.print import fail

PATH_LMP_TEMPLATE = Path(__file__).resolve().parents[3] / 'lmp_utils' / 'template'

# cell constraint of a relaxation case; mirrors the VASP ISIF / Y_CONSTR_LATT choices in use
DICT_BOX_RELAX = {
    'none': '',                                  # fixed cell, ions relaxed (ISIF=2)
    # sigma_zz -> 0 only (ISIF=3 + modeid 3); fixedpoint at the origin, as in
    # stretch_constrained_relax.mod, keeps zlo = 0 and the fractional coordinates intact
    'z':    'fix br all box/relax z 0.0 fixedpoint 0 0 0',
}

TEMPLATE_TOKEN = '_template'


def get_lmp_template_path(name: str = None) -> Path:
    """Return ``lmp_utils/template/<name>``, failing when it does not exist.

    Args:
        name: Template file name, e.g. ``'relax_case_template.in'``.

    Returns:
        Path: Absolute template path.
    """
    path_template = PATH_LMP_TEMPLATE / name
    if not path_template.is_file():
        fail('LAMMPS template not found: %s' % path_template)
    return path_template


def get_lmp_potential_block(pair_style: str = None, pair_coeff: str = None,
                            lmass: list = None) -> str:
    """general_potential.mod with the pair lines and the masses filled in.

    ``include general_mass.mod`` is replaced by the mass lines themselves, which is what
    ``pei_lmp_run_properties`` writes into general_mass.mod at run time.

    Args:
        pair_style: Everything after ``pair_style``, e.g. ``'hdnnp 6.01 dir /abs/potential ...'``.
        pair_coeff: Everything after ``pair_coeff * *``, e.g. ``'Au'``.
        lmass: Full mass commands, e.g. ``['mass 1 196.9666']``.

    Returns:
        str: The rendered block.
    """
    if not pair_style or not pair_coeff or not lmass:
        fail('pair_style, pair_coeff and lmass are all required')
    text = get_lmp_template_path('general_potential.mod').read_text()
    llines = []
    for line in text.splitlines():
        if line.strip() == 'include general_mass.mod':
            llines.extend(lmass)
            continue
        llines.append(line.replace('pair_style_template', pair_style)
                          .replace('pair_coeff_template', pair_coeff))
    return '\n'.join(llines)


def render_lmp_template(text: str = None, dict_block: dict = None,
                        dict_inline: dict = None) -> str:
    """Fill a template: block tokens replace a whole line, inline tokens a word.

    Comment lines are left alone so the documentation of a token in a template header is
    never expanded into live commands.

    Args:
        text: Template text.
        dict_block: ``{token: block}`` for lines consisting of the token alone.
        dict_inline: ``{token: value}`` replaced inside command lines.

    Returns:
        str: Rendered text; fails if a ``*_template`` token survives in a command line.
    """
    dict_block = dict_block or {}
    dict_inline = dict_inline or {}
    llines = []
    for line in text.splitlines():
        if line.lstrip().startswith('#'):
            llines.append(line)
            continue
        if line.strip() in dict_block:
            llines.append(dict_block[line.strip()])
            continue
        for token, value in dict_inline.items():
            line = line.replace(token, str(value))
        if TEMPLATE_TOKEN in line:
            fail('unfilled placeholder in rendered LAMMPS input: %r' % line)
        llines.append(line)
    return '\n'.join(llines) + '\n'


def generate_lmp_relax_input(path_out: Path = None, pair_style: str = None,
                             pair_coeff: str = None, lmass: list = None,
                             box_relax: str = 'none') -> Path:
    """Write a single-structure relaxation input (one minimize, VASP-format outputs).

    Args:
        path_out: Output file, usually ``y_full_relax/lmp.in``.
        pair_style: See :func:`get_lmp_potential_block`.
        pair_coeff: See :func:`get_lmp_potential_block`.
        lmass: See :func:`get_lmp_potential_block`.
        box_relax: Key of :data:`DICT_BOX_RELAX` (``'none'`` or ``'z'``).

    Returns:
        Path: ``path_out``.
    """
    if box_relax not in DICT_BOX_RELAX:
        fail('box_relax must be one of %s, got %r' % (sorted(DICT_BOX_RELAX), box_relax))
    text = get_lmp_template_path('relax_case_template.in').read_text()
    dict_block = {
        'general_init_template': get_lmp_template_path('general_init.mod').read_text().rstrip(),
        'general_potential_template': get_lmp_potential_block(pair_style, pair_coeff, lmass),
        'box_relax_template': DICT_BOX_RELAX[box_relax],
    }
    path_out = Path(path_out)
    path_out.write_text(render_lmp_template(text, dict_block=dict_block))
    return path_out


def generate_lmp_neb_input(path_out: Path = None, pair_style: str = None,
                           pair_coeff: str = None, lmass: list = None,
                           spring: float = 5.0, ftol: float = 0.01, n1: int = 20000,
                           n2: int = 20000, nevery: int = 100,
                           timestep: float = 0.002) -> Path:
    """Write a climbing-image NEB input (fixed cell; replicas from ``data.lmp`` + ``coords.<i>``).

    Args:
        path_out: Output file, usually ``y_full_relax_neb/ini/lmp.in.neb``.
        pair_style: See :func:`get_lmp_potential_block`.
        pair_coeff: See :func:`get_lmp_potential_block`.
        lmass: See :func:`get_lmp_potential_block`.
        spring: Parallel spring constant, eV/A^2 (``fix neb ... parallel neigh``).
        ftol: Force tolerance, eV/A (two-norm of the 3N force per replica).
        n1: Max iterations of the plain NEB stage.
        n2: Max iterations of the climbing-image stage.
        nevery: Output interval of the universe log.
        timestep: FIRE time step, ps.

    Returns:
        Path: ``path_out``.
    """
    text = get_lmp_template_path('neb_case_template.in').read_text()
    dict_block = {
        'general_init_template': get_lmp_template_path('general_init.mod').read_text().rstrip(),
        'general_potential_template': get_lmp_potential_block(pair_style, pair_coeff, lmass),
    }
    dict_inline = {
        'neb_spring_template': '%g' % spring,
        'neb_ftol_template': '%g' % ftol,
        'neb_n1_template': '%d' % n1,
        'neb_n2_template': '%d' % n2,
        'neb_nevery_template': '%d' % nevery,
        'neb_timestep_template': '%g' % timestep,
    }
    path_out = Path(path_out)
    path_out.write_text(render_lmp_template(text, dict_block=dict_block, dict_inline=dict_inline))
    return path_out
