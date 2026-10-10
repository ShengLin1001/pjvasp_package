"""Property structures exported from lmp_utils, without running DFT.

Functions:
    - read_property_dump: Read metal-unit dumps with exact triclinic origins.
    - get_property_frame_rules: Describe attribute-specific frame selection.
    - export_property_structures: Export POSCARs, provenance and coverage states.
"""

from pathlib import Path
import hashlib
import json

import numpy as np
from ase import Atoms
from ase.data import atomic_numbers
from ase.geometry import get_distances
from ase.io import read, write

from mymetal.build.bulk.planar import LIST_SEQSF_FAULT, LIST_SSF_FAULT, LIST_SURF_TYPE
from mymetal.build.bulk.stacking import DICT_STACKING_OFFSET, get_stacking_sequence
from mymetal.post.Cij_energy import read_cij_energy
from mymetal.post.planar import read_ssf, read_surface
from mymetal.post.stretch import my_read_stretch

LSTRAIN_MODE = ['c11', 'c12', 'c13', 'c33', 'c44']
# 数值健康护栏：体相最近邻 2.885 Å，物理构型不会低于 ~2.3 Å；层错能是 meV/atom 量级，
# 每原子 -0.2 eV 只可能是势外推的伪深谷（run 003 prism1w_SF2 实测 -1.6 eV/atom）。
DMIN_OVERLAP_A = 1.2
DE_PER_ATOM_SPURIOUS_EV = -0.2
DICT_FIGURE_SSF = {
    'HCP_basal_I2': 'Basal I2', 'HCP_prism1w_SF2': 'Prism-SF2',
    'HCP_pyr1n_SF1': 'Pyr. I-SF1', 'HCP_pyr1w_SF2': 'Pyr. I',
    'HCP_pyr1w_SF3': 'Pyr. I-SF3', 'HCP_pyr2_SF1': 'Pyr. II',
}
DICT_MODE_PROPERTY = {
    'c11': ['C11', 'C12', 'C13'], 'c12': ['C12'],
    'c13': ['C13'], 'c33': ['C33', 'C13'], 'c44': ['C44'],
}


def fail(msg):
    print('❌ ERROR: ' + str(msg))
    raise SystemExit(1)


def get_file_hash(path_file):
    """Return a SHA256 without loading large trajectories into memory."""
    digest = hashlib.sha256()
    with Path(path_file).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read_property_dump(path_dump: Path = None, lele: list = None) -> list:
    """Read id/type/x/y/z dumps and normalize the true box origin to zero.

    Parse the documented LAMMPS bounds directly, retaining atom ID order
    and the original origin. POSCAR cannot store a displaced box origin,
    so Cartesian positions are translated by that origin, not deformed.

    Args:
        path_dump: Custom text dump in LAMMPS metal units (Å).
        lele: Explicit type-to-element map, e.g. ['Au'].

    Returns:
        list: ASE Atoms with timestep, box_origin_A and atom_ids metadata.
    """
    if not lele:
        raise ValueError('Explicit element mapping is required.')
    lframe = []
    with Path(path_dump).open(encoding='utf-8') as stream:
        while True:
            line = stream.readline()
            if not line:
                break
            if line.strip() != 'ITEM: TIMESTEP':
                raise ValueError('Unexpected dump header: ' + line.strip())
            timestep = int(stream.readline())
            if stream.readline().strip() != 'ITEM: NUMBER OF ATOMS':
                raise ValueError('Missing atom count header.')
            natoms = int(stream.readline())
            lbox = stream.readline().split()[3:]
            if 'abc' in lbox:
                raise ValueError('General triclinic dumps are not supported here.')
            if_tri = lbox[:3] == ['xy', 'xz', 'yz']
            lbound = [[float(value) for value in stream.readline().split()]
                      for _ in range(3)]
            xy, xz, yz = [row[2] for row in lbound] if if_tri else [0., 0., 0.]
            xlo = lbound[0][0] - min(0., xy, xz, xy + xz)
            xhi = lbound[0][1] - max(0., xy, xz, xy + xz)
            ylo = lbound[1][0] - min(0., yz)
            yhi = lbound[1][1] - max(0., yz)
            zlo, zhi = lbound[2][:2]
            cell = np.array([[xhi - xlo, 0., 0.], [xy, yhi - ylo, 0.],
                             [xz, yz, zhi - zlo]])
            origin = np.array([xlo, ylo, zlo])
            lcol = stream.readline().split()[2:]
            dict_col = {name: index for index, name in enumerate(lcol)}
            if not set(['id', 'type', 'x', 'y', 'z']).issubset(dict_col):
                raise ValueError('Expected id/type/x/y/z columns in ' + str(path_dump))
            arr = np.array([[float(value) for value in stream.readline().split()]
                            for _ in range(natoms)])
            arr = arr[np.argsort(arr[:, dict_col['id']])]
            lid = arr[:, dict_col['id']].astype(int)
            ltype = arr[:, dict_col['type']].astype(int)
            if len(set(lid)) != natoms or min(ltype) < 1 or max(ltype) > len(lele):
                raise ValueError('Invalid atom IDs/types in ' + str(path_dump))
            positions = arr[:, [dict_col[key] for key in ['x', 'y', 'z']]] - origin
            atoms = Atoms([lele[value - 1] for value in ltype], positions=positions,
                          cell=cell, pbc=[value == 'pp' for value in lbox[-3:]])
            atoms.info.update(timestep=timestep, box_origin_A=origin.tolist(),
                              atom_ids=lid.tolist())
            check_property_atoms(atoms)
            lframe.append(atoms)
    if not lframe:
        raise ValueError('Empty trajectory: ' + str(path_dump))
    return lframe


def check_property_atoms(atoms):
    """Reject nonperiodic, empty, nonfinite or inverted export geometries."""
    if not len(atoms) or not np.all(atoms.pbc):
        raise ValueError('Expected a nonempty, three-periodic cell (including vacuum).')
    if not np.all(np.isfinite(atoms.positions)) or not np.all(np.isfinite(atoms.cell)):
        raise ValueError('Nonfinite geometry.')
    if np.linalg.det(atoms.cell) <= 0:
        raise ValueError('Nonpositive cell volume.')


def get_property_frame_rules(strain_points: int = 5) -> list:
    """Map each attribute family to its actual output, not its last frame.

    Args:
        strain_points: Odd number of uniformly selected strain points; None
            exports all points. Raw scans always remain in their source directory.

    Returns:
        list: Records with phase/group/model/path, states and normalization.
    """
    if strain_points is not None and (strain_points < 3 or strain_points % 2 != 1):
        raise ValueError('strain_points must be odd and >= 3, or None.')
    lrule = []
    for phase in ['fcc', 'hcp']:
        lrule.append(dict(phase=phase, group='stretch', model='bulk',
                          source='y_stretch/' + phase + '/dump/data.full_relax',
                          states=['equilibrium'], format='lammps-data'))
        for mode in LSTRAIN_MODE:
            lrule.append(dict(phase=phase, group='cij', model=mode,
                              source='y_Cij_energy/' + phase + '/dump/y_cij_energy_' + mode + '/movie.lammpstrj',
                              states=['strain'], strain_points=strain_points,
                              format='lammps-dump-text'))
    for phase, model in LIST_SURF_TYPE:
        lrule.append(dict(phase=phase, group='surface', model=model,
                          source='y_surface/' + phase + '/' + model + '/dump/movie.lammpstrj',
                          states=['reference', 'unrelaxed', 'relaxed'], nsurface=2,
                          format='lammps-dump-text'))
    for phase, model, fault, t1, t2 in LIST_SSF_FAULT:
        name = model + '_' + fault
        lrule.append(dict(phase=phase, group='ssf', model=name,
                          source='y_ssf/' + phase + '/' + name + '/dump/movie.lammpstrj',
                          ideal='y_ssf/' + phase + '/' + name + '/dump/ideal.lammpstrj',
                          states=['reference', 'relaxed'], nfault=1,
                          initial_t_cartesian_fraction=[t1, t2], format='lammps-dump-text'))
    for phase, name, n_repeat in LIST_SEQSF_FAULT:
        _, seq_ref, seq_fault = get_stacking_sequence(name, n_repeat)
        lrule.append(dict(phase=phase, group='seqsf', model=name,
                          source='y_seqsf/' + phase + '/' + name + '/dump/movie.lammpstrj',
                          ideal='y_seqsf/' + phase + '/' + name + '/dump/ideal.lammpstrj',
                          states=['reference', 'relaxed'], nfault=2,
                          sequence_reference=seq_ref, sequence_fault=seq_fault,
                          format='lammps-dump-text'))
    return lrule


def read_property_data(path_data, lele):
    """Read atomic-style data, preserving IDs and subtracting the true origin."""
    atoms = read(str(path_data), format='lammps-data', style='atomic', sort_by_id=True,
                 Z_of_type={index + 1: atomic_numbers[element]
                            for index, element in enumerate(lele)})
    dict_origin = {}
    for line in Path(path_data).read_text(encoding='utf-8').splitlines():
        lfield = line.split('#', 1)[0].split()
        if len(lfield) == 4 and lfield[-2:] in [['xlo', 'xhi'], ['ylo', 'yhi'], ['zlo', 'zhi']]:
            dict_origin[lfield[-2][0]] = float(lfield[0])
    origin = np.array([dict_origin[key] for key in ['x', 'y', 'z']])
    atoms.positions -= origin
    atoms.info.update(box_origin_A=origin.tolist(), atom_ids=atoms.arrays['id'].tolist())
    check_property_atoms(atoms)
    return atoms


def read_rule_frames(path_root, dict_rule, lele):
    """Read each family with strict frame count and state-index validation."""
    path_source = path_root / dict_rule['source']
    if dict_rule['format'] == 'lammps-data':
        atoms = read_property_data(path_source, lele)
        return [('equilibrium', 0, atoms, path_source, None)]
    lframe = read_property_dump(path_source, lele)
    path_data = path_source.parent / 'y_post_data.txt'
    arr = np.loadtxt(path_data, skiprows=1, usecols=(0, 1))
    if arr.shape != (len(lframe), 2) or not np.all(np.isfinite(arr)):
        raise ValueError('Trajectory/energy table mismatch: ' + str(path_source))
    lselected = []
    if dict_rule['group'] == 'cij':
        strain_points = dict_rule['strain_points']
        if len(lframe) < 3 or len(lframe) % 2 != 1:
            raise ValueError('Expected odd, symmetric elastic scan.')
        lindex = (list(range(len(lframe))) if strain_points is None else
                  sorted(set(np.rint(np.linspace(0, len(lframe) - 1,
                                                 min(strain_points, len(lframe)))).astype(int))))
        if not np.allclose(arr[:, 0] - 1., -(arr[::-1, 0] - 1.), atol=1e-10):
            raise ValueError('Asymmetric strain scan.')
        for index in lindex:
            lselected.append(('strain', int(index), lframe[index], path_source,
                              float(arr[index, 0] - 1.)))
    else:
        if len(lframe) != len(dict_rule['states']):
            raise ValueError('Wrong state count: ' + str(path_source))
        if not np.allclose(arr[:, 0], np.arange(len(lframe)), atol=1e-10):
            raise ValueError('State IDs must be 0/1[/2].')
        if not all(np.allclose(atoms.cell[:2], lframe[0].cell[:2], atol=1e-8, rtol=0.)
                   for atoms in lframe):
            raise ValueError('In-plane cell changed.')
        if dict_rule['group'] == 'surface':
            if not all(len(atoms) == len(lframe[0]) for atoms in lframe):
                raise ValueError('Surface atom counts changed.')
            # run 0 may wrap atoms by an in-plane lattice vector; that is not relaxation.
            delta = (lframe[1].positions - lframe[0].positions) @ np.linalg.inv(lframe[0].cell)
            delta[:, :2] -= np.rint(delta[:, :2])
            if not np.allclose(delta @ lframe[0].cell[:], 0., atol=1e-8, rtol=0.):
                raise ValueError('Unrelaxed surface atoms moved during cleavage.')
            if not np.allclose(lframe[1].cell, lframe[2].cell, atol=1e-8, rtol=0.):
                raise ValueError('Surface cell relaxed.')
        lselected = [(state, index, lframe[index], path_source, None)
                     for index, state in enumerate(dict_rule['states'])]
    if dict_rule['group'] == 'cij':
        path_equil = path_root / 'y_stretch' / dict_rule['phase'] / 'dump/data.full_relax'
        atoms_equil = read(str(path_equil), format='lammps-data', style='atomic')
        if not np.allclose(lframe[len(lframe) // 2].cell, atoms_equil.cell, atol=1e-5, rtol=0.):
            raise ValueError('Elastic zero strain does not use this potential\'s equilibrium cell.')
    if 'ideal' in dict_rule:
        path_ideal = path_root / dict_rule['ideal']
        lideal = read_property_dump(path_ideal, lele)
        if len(lideal) != 1 or len(lideal[0]) != len(lframe[-1]):
            raise ValueError('Ideal/relaxed fault atom count mismatch.')
        lselected.append(('ideal', 0, lideal[0], path_ideal, None))
    return lselected


def get_stacking_letters(atoms):
    """Identify one-atom-per-layer ABC sequences, modulo rigid in-plane drift."""
    lindex = np.argsort(atoms.positions[:, 2])
    frac = atoms.get_scaled_positions(wrap=False)[lindex, :2]
    loffset = np.array(list(DICT_STACKING_OFFSET.values()))
    lletter = list(DICT_STACKING_OFFSET)
    delta = frac - frac[0]
    lsequence = []
    for position in delta:
        diff = position - loffset
        diff -= np.rint(diff)
        ldistance = np.linalg.norm(diff, axis=1)
        if min(ldistance) > 0.05:
            raise ValueError('Layer registry no longer matches an ABC sublattice.')
        lsequence.append(lletter[int(np.argmin(ldistance))])
    return ''.join(lsequence)


def get_periodic_fault_count(sequence, phase):
    """Count local fcc layers in HCP, or local hcp layers in FCC."""
    local = [sequence[(index - 1) % len(sequence)] == sequence[(index + 1) % len(sequence)]
             for index in range(len(sequence))]
    return sum(local) if phase == 'fcc' else len(sequence) - sum(local)


def get_min_distance(atoms) -> float:
    """Smallest PBC interatomic distance; catches collapsed/overlapped exports."""
    if len(atoms) < 2:
        return float('inf')
    dist = get_distances(atoms.positions, cell=atoms.cell, pbc=True)[1]
    np.fill_diagonal(dist, np.inf)
    return float(dist.min())


def check_poscar_roundtrip(atoms, path_poscar):
    """Check cell, species, PBC and all fractional coordinates after conversion."""
    atoms_back = read(str(path_poscar), format='vasp')
    if atoms.get_chemical_symbols() != atoms_back.get_chemical_symbols():
        raise ValueError('POSCAR species/atom order changed.')
    if not np.array_equal(atoms.pbc, atoms_back.pbc):
        raise ValueError('POSCAR PBC changed.')
    dcell = float(np.max(np.abs(atoms.cell[:] - atoms_back.cell[:])))
    dfrac = atoms.get_scaled_positions(wrap=False) - atoms_back.get_scaled_positions(wrap=False)
    dfrac -= np.rint(dfrac)
    dpos = float(np.max(np.abs(dfrac @ atoms.cell[:])))
    if dcell > 1e-9 or dpos > 1e-9:
        raise ValueError('POSCAR geometry round-trip failed.')
    return dict(cell_error_A=dcell, position_error_A=dpos)


def generate_candidate(atoms, path_out, candidate_id, dict_meta, path_vasp_static=None):
    """Write only new files; existing candidates must have identical provenance."""
    path_candidate = path_out / 'y_dir' / candidate_id
    path_candidate.mkdir(parents=True, exist_ok=True)
    path_poscar = path_candidate / 'POSCAR'
    if path_poscar.exists():
        dict_old = json.loads((path_candidate / 'source.json').read_text(encoding='utf-8'))
        if dict_old['source_sha256'] != dict_meta['source_sha256']:
            raise ValueError('Candidate source changed; use a new output: ' + candidate_id)
    else:
        write(str(path_poscar), atoms, format='vasp', direct=True, sort=False, vasp5=True)
    dict_check = check_poscar_roundtrip(atoms, path_poscar)
    dict_geometry = dict(symbols=atoms.get_chemical_symbols(),
                          cell_A=np.round(atoms.cell[:], 10).tolist(),
                          fractional=np.round(atoms.get_scaled_positions(), 10).tolist())
    fingerprint = hashlib.sha256(json.dumps(dict_geometry, sort_keys=True).encode('utf-8')).hexdigest()
    dict_meta.update(geometry_sha256=fingerprint, poscar_sha256=get_file_hash(path_poscar),
                     min_distance_A=get_min_distance(atoms),
                     candidate_id=candidate_id, poscar=str(path_poscar.resolve()),
                     natoms=len(atoms), symbols=atoms.get_chemical_symbols(),
                     cell_A=atoms.cell[:].tolist(), volume_A3=float(atoms.get_volume()),
                     area_A2=float(np.linalg.norm(np.cross(atoms.cell[0], atoms.cell[1]))),
                     pbc=atoms.pbc.tolist(), box_origin_A=atoms.info.get('box_origin_A', [0., 0., 0.]),
                     atom_ids=atoms.info.get('atom_ids', list(range(1, len(atoms) + 1))),
                     timestep=atoms.info.get('timestep'), units='angstrom', **dict_check)
    if path_vasp_static is not None:
        for name in ['INCAR', 'KPOINTS', 'POTCAR']:
            path_source = path_vasp_static / name
            path_dest = path_candidate / name
            if path_dest.exists():
                if get_file_hash(path_dest) != get_file_hash(path_source):
                    raise ValueError('Existing VASP input differs: ' + str(path_dest))
            else:
                path_dest.write_bytes(path_source.read_bytes())
        dict_meta['dft_status'] = 'input_only_not_submitted'
    (path_candidate / 'source.json').write_text(json.dumps(dict_meta, ensure_ascii=False, indent=2) + '\n',
                                               encoding='utf-8')
    return dict_meta


def get_rule_property_keys(dict_rule, state):
    """Use explicit phase/property keys so the 18 figure attributes are traceable."""
    phase, group, model = [dict_rule[key] for key in ['phase', 'group', 'model']]
    if group == 'stretch':
        return [phase + '_' + key for key in ['a', 'c', 'c_over_a', 'E0']]
    if group == 'cij':
        return [phase + '_' + key for key in DICT_MODE_PROPERTY[model]]
    if group == 'surface':
        lkey = ['surf_unrel_' + model, 'surf_rel_' + model]
        return lkey if state == 'reference' else [lkey[0 if state == 'unrelaxed' else 1]]
    return ['ssf_' + model]


def get_rule_properties(path_root, dict_rule):
    """Read property scalars and retain unstable states explicitly."""
    phase, group, model = [dict_rule[key] for key in ['phase', 'group', 'model']]
    if group == 'stretch':
        result = my_read_stretch(path_root / 'y_stretch' / phase / 'p_post_stretch.txt')
        return dict(a_fit_A=float(result[-1][-1][0]), c_fit_A=float(result[-1][-1][2]),
                    energy_fit_eV_atom=float(result[-1][2]))
    if group == 'cij':
        return read_cij_energy(path_root / 'y_Cij_energy' / phase / 'y_post_cij_energy.txt')
    path_case = path_root / ('y_' + group) / phase / model
    if group == 'surface':
        return read_surface(path_case / 'y_post_surface.txt')
    return read_ssf(path_case / 'y_post_ssf.txt')


def check_static_template(path_vasp_static):
    """Reject relaxation templates: labeling must not replace the candidate geometry."""
    for name in ['INCAR', 'KPOINTS', 'POTCAR']:
        if not (path_vasp_static / name).is_file():
            raise ValueError('Missing static template: ' + name)
    dict_tag = {}
    for line in (path_vasp_static / 'INCAR').read_text(encoding='utf-8').splitlines():
        line = line.split('#', 1)[0].split('!', 1)[0]
        for field in line.split(';'):
            if '=' in field:
                key, value = field.split('=', 1)
                dict_tag[key.strip().upper()] = value.strip()
    if int(dict_tag.get('NSW', '-99')) != 0 or int(dict_tag.get('IBRION', '-99')) != -1:
        raise ValueError('VASP labeling requires NSW=0 and IBRION=-1.')


def export_property_structures(path_root: Path = None, path_out: Path = None,
                               lele: list = None, run: str = None, epoch: int = None,
                               round_id: str = None, strain_points: int = 5,
                               path_vasp_static: Path = None) -> dict:
    """Export the full FCC/HCP property set and report every failed family.

    Args:
        path_root: One potential's completed lmp_utils calculation directory.
        path_out: Independent export directory; source outputs are never modified.
        lele: Explicit LAMMPS element order.
        run: Committee member identity.
        epoch: Common weight snapshot number, not an implicit final epoch.
        round_id: Source training round identity.
        strain_points: Number of selected points per mode (raw scans retained).
        path_vasp_static: Optional INCAR/KPOINTS/POTCAR, input smoke test only.

    Returns:
        dict: Manifest records and per-family successful/unstable/failed states.
    """
    ### check ###
    for path_value in [path_root, path_out]:
        if path_value is None or not Path(path_value).is_absolute():
            raise ValueError('Property source/output paths must be absolute.')
    path_root, path_out = Path(path_root).resolve(), Path(path_out).resolve()
    if path_vasp_static is not None:
        path_vasp_static = Path(path_vasp_static).resolve()
        check_static_template(path_vasp_static)
    lrule = get_property_frame_rules(strain_points)
    dict_potential_hash = {path_file.name: get_file_hash(path_file)
                           for path_file in sorted((path_root / 'potential').iterdir())
                           if path_file.is_file()}
    lseed = [line.split('#', 1)[0].split()[1] for line in
             (path_root / 'potential/input.nn').read_text(encoding='utf-8').splitlines()
             if line.split('#', 1)[0].split()[:1] == ['random_seed']]
    seed = int(lseed[0]) if lseed else None
    # to here

    ### prepare ###
    path_out.mkdir(parents=True, exist_ok=True)
    dict_identity = dict(source=str(path_root), potential_sha256=dict_potential_hash,
                         run=run, epoch=epoch, round=round_id, strain_points=strain_points,
                         elements=lele, vasp_template=str(path_vasp_static))
    path_identity = path_out / 'identity.json'
    if path_identity.exists() and json.loads(path_identity.read_text(encoding='utf-8')) != dict_identity:
        raise ValueError('Export identity changed; use an independent output directory.')
    path_identity.write_text(json.dumps(dict_identity, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    lrecord, lstatus, lfailed = [], [], []
    # to here

    ### main ###
    for dict_rule in lrule:
        phase, group, model = [dict_rule[key] for key in ['phase', 'group', 'model']]
        label = phase + '/' + group + '/' + model
        try:
            dict_values = get_rule_properties(path_root, dict_rule)
            if not all(np.isfinite(value) for value in dict_values.values()):
                raise ValueError('Nonfinite/missing property result.')
            lselected = read_rule_frames(path_root, dict_rule, lele)
            lmodel_record = []
            for state, index, atoms, path_source, strain in lselected:
                candidate_id = label + '/' + state + ('_%03d' % index if state == 'strain' else '')
                dict_meta = dict(round=round_id, run=run, epoch=epoch, seed=seed, phase=phase,
                                 property_group=group, model=model, state=state, frame=index,
                                 source_path=str(path_source), source_sha256=get_file_hash(path_source),
                                 property_keys=get_rule_property_keys(dict_rule, state),
                                 potential_sha256=dict_potential_hash, strain=strain,
                                 stable=dict_values.get('stable'),
                                 selection_reason='attribute_specific_state_not_last_frame',
                                 normalization='per_atom_bulk_reference' if group in ['ssf', 'seqsf'] else None)
                for key in ['nfault', 'nsurface', 'sequence_reference', 'sequence_fault',
                            'initial_t_cartesian_fraction']:
                    if key in dict_rule:
                        dict_meta[key] = dict_rule[key]
                if group == 'cij':
                    dict_meta['constrains'] = DICT_MODE_PROPERTY[model]
                if group == 'seqsf':
                    sequence = get_stacking_letters(atoms)
                    expected = dict_rule['sequence_reference' if state == 'reference' else 'sequence_fault']
                    if state in ['reference', 'ideal'] and sequence != expected:
                        raise ValueError('Stacking sequence does not match constructed model.')
                    if state == 'relaxed' and sequence not in (expected * 2):
                        dict_meta['geometry_status'] = 'stacking_changed'
                    else:
                        dict_meta['geometry_status'] = 'stacking_preserved'
                    dict_meta['sequence_actual'] = sequence
                    dict_meta['local_nonbulk_layers'] = get_periodic_fault_count(sequence, phase)
                lmodel_record.append(generate_candidate(atoms, path_out, candidate_id, dict_meta, path_vasp_static))
            lrecord.extend(lmodel_record)
            # 状态优先级：伪深谷 > 原子重叠 > 层序改变 > 无稳定层错 > 正常。
            # 塌缩/外推构型是独立诊断项，不能因为 stable=yes 就当正常层错能读。
            if group in ('ssf', 'seqsf'):
                natoms_fault = len(next(atoms for state, _, atoms, _, _ in lselected if state == 'relaxed'))
                dE_per_atom = float(dict_values.get('dE', float('nan'))) / natoms_fault
                dict_values['dE_per_atom'] = dE_per_atom
            else:
                dE_per_atom = None
            if dE_per_atom is not None and dE_per_atom < DE_PER_ATOM_SPURIOUS_EV:
                status = 'spurious_basin'
            elif any(record['min_distance_A'] < DMIN_OVERLAP_A for record in lmodel_record):
                status = 'overlap'
            elif any(record.get('geometry_status') == 'stacking_changed' for record in lmodel_record):
                status = 'stacking_changed'
            else:
                status = 'unstable' if dict_values.get('stable') is False else 'ok'
            lstatus.append(dict(phase=phase, group=group, model=model, status=status,
                                values=dict_values, candidates=len(lmodel_record)))
            print(('✅ ' if status == 'ok' else '⚠️ ') + label + ': ' + str(len(lmodel_record)) + ' candidates; ' + status)
        except (Exception, SystemExit) as error:
            lfailed.append(label)
            lstatus.append(dict(phase=phase, group=group, model=model, status='failed', error=str(error)))
            print('❌ ' + label + ': ' + str(error))
    (path_out / 'frame_rules.json').write_text(json.dumps(lrule, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    (path_out / 'manifest.jsonl').write_text(''.join(json.dumps(record, ensure_ascii=False) + '\n' for record in lrecord), encoding='utf-8')
    dict_summary = dict(identity=dict_identity, candidates=len(lrecord), families=lstatus,
                        failed=lfailed, unstable=[item['model'] for item in lstatus if item['status'] == 'unstable'],
                        stacking_changed=[item['model'] for item in lstatus if item['status'] == 'stacking_changed'],
                        spurious_basin=[item['model'] for item in lstatus if item['status'] == 'spurious_basin'],
                        overlap=[item['model'] for item in lstatus if item['status'] == 'overlap'])
    (path_out / 'summary.json').write_text(json.dumps(dict_summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print('================ 📊 Summary: %d families, %d candidates, %d failures, %d unstable, '
          '%d spurious_basin, %d overlap' %
          (len(lrule), len(lrecord), len(lfailed), len(dict_summary['unstable']),
           len(dict_summary['spurious_basin']), len(dict_summary['overlap'])))
    if lfailed:
        fail('Failed families: ' + ', '.join(lfailed))
    return dict_summary
    # to here
