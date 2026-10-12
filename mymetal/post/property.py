"""Property structures exported from lmp_utils, without running DFT.

Functions:
    - read_property_dump: Read metal-unit dumps with exact triclinic origins.
    - get_property_frame_rules: Describe attribute-specific frame selection.
    - get_frame_health: Flag spurious-basin / overlapping frames.
    - export_property_structures: Export POSCARs, provenance and coverage states.
    - check_same_geometry: Same cell and atoms up to permutation/rigid shift.
    - check_same_crystal: Same infinite crystal even in a differently tilted cell.
    - get_translation_count: Pure lattice translations mapping a frame onto itself.
    - get_frame_action: Decide whether one exported frame is sent to DFT.
    - select_property_candidates: Deduplicate exports into one DFT to-do list.
    - find_fault_lost_selected: Re-check a frozen to-do list for vanished faults.
    - get_static_outcar_status: Classify one static VASP OUTCAR (finished, EDIFF, single step).
    - read_static_dft_frame: POSCAR geometry + OUTCAR energy/forces as one labelled frame.
"""

from pathlib import Path
import hashlib
import json

import numpy as np
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.data import atomic_numbers
from ase.geometry import get_distances
from ase.io import read, write

from mymetal.build.bulk.planar import LIST_SEQSF_FAULT, LIST_SSF_FAULT, LIST_SURF_TYPE
from mymetal.build.bulk.stacking import DICT_STACKING_OFFSET, get_stacking_sequence
from mymetal.post.Cij_energy import read_cij_energy
from mymetal.post.planar import read_ssf, read_surface
from mymetal.post.stretch import my_read_stretch

LSTRAIN_MODE = ['c11', 'c12', 'c13', 'c33', 'c44']
# 数值健康护栏，分界取自 2026-10-10 六势 epoch-350 实测：
# - 弛豫帧（stretch 平衡、表面/层错 relaxed）：正常 d_min >= 2.63 Å、E/atom >= 本势基态；
#   伪深谷 d_min 1.33–1.88 Å、比本势基态低 1.6–2.9 eV/atom。故弛豫帧 d_min < 2.3 Å 或
#   低于基态 0.2 eV/atom 即判伪深谷。按每帧总能量判，不按层错 dE 摊到全胞（会被稀释）。
# - 构造帧（reference/unrelaxed/strain/ideal）是刚性搭出来的，ideal 最短可到 1.95 Å
#   （pyr1n -0.11），不能套弛豫帧阈值；只设硬下限：Au PAW RCORE 1.32 Å，再近即核重叠。
DMIN_RELAXED_A = 2.3
DMIN_HARD_A = 1.5
DE_BELOW_GROUND_SPURIOUS_EV = 0.2
LSTATE_RELAXED = ['equilibrium', 'relaxed']
# 选帧去重：同一条势的零应变帧/共用参考帧实测偏差 <= 3e-10 Å；不同势的胞最近只差
# 1.4e-5 Å（epoch-350 六势，FCC 平衡胞）。取 1e-6 Å：只并真重复，不会把不同势、
# 不同胞应变或不同层错位置误并。
DTOL_DEDUP_A = 1e-6
# 层错是否还在：弛豫帧比同族未弛豫帧多出纯平移 = 层错面消失、整胞变成另一种晶体。
# 层错/表面的原子错排是 Å 量级，1e-3 Å 只认 LAMMPS 弛豫留下的数值噪声。2026-10-11 round 0
# 实测：prism1w_SF2 relaxed（5/6 势）与 pyr2_SF1 relaxed（2/6 势）从 1 个平移变成 20 个
# （1 原子原胞 Bravais 晶格，CN=10），其余层错族 relaxed 与 ideal 平移数相同。
DTOL_TRANSLATION_A = 1e-3
# 2026-10-10 PJ：Pyr. I-SF1 的 t2 初猜存疑（-0.11 下 4/6 回到完整晶格、2/6 伪深谷），
# 先用 DFT 定位层错位置再重跑，本族的**弛豫帧**暂缓送算（round 0 用的是 ['HCP_pyr1n_SF1']）。
# 2026-10-11 DFT 定位在 t2≈-0.17，planar.py 初猜已改为 -0.17，PJ 同意从 round 1 起取消暂缓。
# 暂缓机制保留：以后有族要暂缓时在这里加，或调用处显式传 ldefer_model。
LDEFER_MODEL = []
# DFT 单点回读：OUTCAR 的离子坐标只印 5 位小数，1e-4 Å 足以区分「同一几何的舍入」与
# 「原子顺序 / 胞对不上」（后者是 Å 量级）。round 0 实测 POSCAR 与 OUTCAR 坐标最大差 4.7e-6 Å。
DTOL_OUTCAR_POS_A = 1e-4
# get_static_outcar_status 的返回值；只有 'ok' 可以进训练集。
DICT_OUTCAR_STATUS_NOTE = {
    'ok': '正常结束、电子步收敛、单个离子步',
    'missing': '没有 OUTCAR 或为空',
    'unfinished': '没有 General timing（中途被杀或还在跑）',
    'ediff_not_reached': '正常结束但没有 EDIFF is reached（NELM 用尽，能量/力不可信）',
    'not_static': 'TOTAL-FORCE 块不是恰好一个（不是单点，末步几何已不是 POSCAR）',
}
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
    # 每帧 NNP 总能量随帧走，供伪深谷「相对本势基态」判据用；ideal 帧没有能量记录。
    for atoms, energy in zip(lframe, arr[:, 1]):
        atoms.info['energy_nnp_eV'] = float(energy)
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


def get_frame_health(state, energy_atom, min_distance, energy_ground) -> list:
    """Flag one exported frame against the spurious-basin and overlap guards.

    Args:
        state: Exported state name (``relaxed``/``ideal``/``strain``...).
        energy_atom: NNP energy per atom (eV), or None when not recorded (ideal).
        min_distance: Smallest PBC interatomic distance (Å).
        energy_ground: This potential's lowest bulk E0 (eV/atom), or None.

    Returns:
        list: Empty when healthy; otherwise ``below_ground`` / ``dmin_relaxed`` /
        ``dmin_hard`` flags.
    """
    lflag = []
    if (energy_atom is not None and energy_ground is not None
            and energy_atom < energy_ground - DE_BELOW_GROUND_SPURIOUS_EV):
        lflag.append('below_ground')
    if state in LSTATE_RELAXED and min_distance < DMIN_RELAXED_A:
        lflag.append('dmin_relaxed')
    elif min_distance < DMIN_HARD_A:
        lflag.append('dmin_hard')
    return lflag


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
    # 伪深谷按「本势自己的基态」判，须在逐族循环前拿到两相 E0（stretch 族排在 cij 之间）。
    # 读不到时只关掉能量判据，stretch 族本身会在主循环里记 failed，整批最终 fail。
    lenergy_e0 = []
    for dict_rule in lrule:
        if dict_rule['group'] == 'stretch':
            try:
                lenergy_e0.append(get_rule_properties(path_root, dict_rule)['energy_fit_eV_atom'])
            except (Exception, SystemExit) as error:
                print('⚠️  ground-state energy unavailable (' + dict_rule['phase'] + '): ' + str(error))
    energy_ground = min(lenergy_e0) if len(lenergy_e0) == 2 else None
    print('📍 ground-state energy: ' + ('%.8f eV/atom' % energy_ground if energy_ground is not None
                                        else 'unavailable, energy guard disabled'))
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
                # stretch 平衡帧来自 data 文件、无逐帧能量，取同一次弛豫的拟合 E0。
                energy = (dict_values['energy_fit_eV_atom'] * len(atoms) if state == 'equilibrium'
                          else atoms.info.get('energy_nnp_eV'))
                dict_meta['energy_nnp_eV'] = energy
                dict_meta['energy_nnp_eV_atom'] = None if energy is None else energy / len(atoms)
                dict_meta['energy_ground_eV_atom'] = energy_ground
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
                if group == 'ssf' and state == 'ideal':
                    atoms_ref = next(item[2] for item in lselected if item[0] == 'reference')
                    dict_meta['geometry_status'] = ('perfect_lattice' if check_same_crystal(atoms_ref, atoms)
                                                    else 'fault')
                if group in ('ssf', 'seqsf') and state == 'relaxed' \
                        and dict_meta.get('geometry_status') != 'stacking_changed':
                    atoms_ideal = next((item[2] for item in lselected if item[0] == 'ideal'), None)
                    if atoms_ideal is not None \
                            and get_translation_count(atoms) > get_translation_count(atoms_ideal):
                        dict_meta['geometry_status'] = 'fault_lost'
                dict_meta['health'] = get_frame_health(state, dict_meta['energy_nnp_eV_atom'],
                                                       get_min_distance(atoms), energy_ground)
                lmodel_record.append(generate_candidate(atoms, path_out, candidate_id, dict_meta, path_vasp_static))
            lrecord.extend(lmodel_record)
            # 状态优先级：伪深谷 > 原子重叠 > 层序改变 > 无稳定层错 > 层错面消失 > 正常。
            # 塌缩/外推构型是独立诊断项，不能因为 stable=yes 就当正常层错能读。
            if group in ('ssf', 'seqsf'):
                natoms_fault = len(next(atoms for state, _, atoms, _, _ in lselected if state == 'relaxed'))
                dict_values['dE_per_atom'] = float(dict_values.get('dE', float('nan'))) / natoms_fault
            lflag = [flag for record in lmodel_record for flag in record['health']]
            if 'below_ground' in lflag or 'dmin_relaxed' in lflag:
                status = 'spurious_basin'
            elif 'dmin_hard' in lflag:
                status = 'overlap'
            elif any(record.get('geometry_status') == 'stacking_changed' for record in lmodel_record):
                status = 'stacking_changed'
            elif dict_values.get('stable') is False:
                status = 'unstable'
            elif any(record.get('geometry_status') == 'fault_lost' for record in lmodel_record):
                status = 'fault_lost'
            else:
                status = 'ok'
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
                        energy_ground_eV_atom=energy_ground,
                        failed=lfailed, unstable=[item['model'] for item in lstatus if item['status'] == 'unstable'],
                        stacking_changed=[item['model'] for item in lstatus if item['status'] == 'stacking_changed'],
                        fault_lost=[item['model'] for item in lstatus if item['status'] == 'fault_lost'],
                        spurious_basin=[item['model'] for item in lstatus if item['status'] == 'spurious_basin'],
                        overlap=[item['model'] for item in lstatus if item['status'] == 'overlap'])
    (path_out / 'summary.json').write_text(json.dumps(dict_summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print('================ 📊 Summary: %d families, %d candidates, %d failures, %d unstable, '
          '%d spurious_basin, %d overlap, %d fault_lost' %
          (len(lrule), len(lrecord), len(lfailed), len(dict_summary['unstable']),
           len(dict_summary['spurious_basin']), len(dict_summary['overlap']),
           len(dict_summary['fault_lost'])))
    if lfailed:
        fail('Failed families: ' + ', '.join(lfailed))
    return dict_summary
    # to here


DICT_ACTION_NOTE = {
    'selected': '送 DFT 单点',
    'deferred': '弛豫结果暂缓：层错位置存疑，先用 DFT 定位再重跑',
    'spurious_basin': '弛豫帧落入伪深谷，同族未弛豫帧（ideal/unrelaxed）代替送算',
    'overlap': '构造帧原子核重叠，不能送 DFT',
    'returned_to_lattice': '无稳定层错：弛豫帧回到完整晶格，与 reference 等价',
    'perfect_lattice': '未弛豫帧的错排矢量是晶格平移，即完整晶格换了胞，与 reference 等价',
    'stacking_changed': '弛豫后层序改变，同族未弛豫帧代替送算',
    'fault_lost': '弛豫后层错面消失（比同族未弛豫帧多出纯平移，整胞变成另一种晶体），'
                  '不是层错构型，同族未弛豫帧代替送算',
}
LSOURCE_KEY = ['round', 'run', 'epoch', 'seed', 'phase', 'property_group', 'model', 'state',
               'frame', 'strain', 'property_keys', 'stable', 'candidate_id', 'natoms',
               'min_distance_A', 'energy_nnp_eV_atom', 'energy_ground_eV_atom', 'health',
               'geometry_sha256', 'poscar', 'poscar_sha256']


def check_same_geometry(atoms_a, atoms_b, tol_A: float = DTOL_DEDUP_A) -> bool:
    """Whether two frames are one DFT calculation: same cell, same atoms.

    Atoms may be permuted and rigidly translated (energy and forces are
    invariant), but the cell matrix must match: different strains or fault
    positions are never merged.

    Args:
        atoms_a: First frame.
        atoms_b: Second frame.
        tol_A: Cartesian tolerance for cell vectors and atom positions (Å).

    Returns:
        bool: True when every atom of ``atoms_a`` maps one-to-one onto ``atoms_b``.
    """
    if len(atoms_a) != len(atoms_b):
        return False
    lsymbol_a = np.array(atoms_a.get_chemical_symbols())
    lsymbol_b = np.array(atoms_b.get_chemical_symbols())
    if sorted(lsymbol_a) != sorted(lsymbol_b):
        return False
    cell = atoms_a.cell[:]
    if not np.allclose(cell, atoms_b.cell[:], atol=tol_A, rtol=0.):
        return False
    frac_a = atoms_a.get_scaled_positions(wrap=False)
    frac_b = atoms_b.get_scaled_positions(wrap=False)
    mask_species = lsymbol_a[:, None] != lsymbol_b[None, :]
    for index_b in np.where(lsymbol_b == lsymbol_a[0])[0]:
        diff = frac_a[:, None, :] + (frac_b[index_b] - frac_a[0]) - frac_b[None, :, :]
        diff -= np.rint(diff)
        dist = np.linalg.norm(diff @ cell, axis=2)
        dist[mask_species] = np.inf
        lmatch = dist.argmin(axis=1)
        if np.all(dist[np.arange(len(lmatch)), lmatch] < tol_A) and len(set(lmatch)) == len(lmatch):
            return True
    return False


def check_same_crystal(atoms_a, atoms_b, tol_A: float = DTOL_DEDUP_A) -> bool:
    """Whether two cells describe one infinite crystal, without any rigid shift.

    Catches tilted-cell "faults" whose tilt is a lattice translation, e.g.
    FCC (100) at t = (0.5, 0.5) shifts by a/2<011>: the perfect crystal in
    another cell. check_same_geometry cannot see this (the cells differ).

    Args:
        atoms_a: Reference frame.
        atoms_b: Frame to test (same atom count and volume).
        tol_A: Cartesian tolerance (Å).

    Returns:
        bool: True when every atom of ``atoms_b`` is a site of ``atoms_a`` and
        every cell vector of ``atoms_b`` is a translation symmetry of ``atoms_a``.
    """
    if len(atoms_a) != len(atoms_b) or not np.isclose(atoms_a.get_volume(), atoms_b.get_volume(),
                                                       rtol=1e-8, atol=0.):
        return False
    cell = atoms_a.cell[:]
    frac_a = atoms_a.get_scaled_positions(wrap=False)
    lsymbol_a = np.array(atoms_a.get_chemical_symbols())

    def check_on_sites(frac_point, lsymbol_point):
        diff = frac_point[:, None, :] - frac_a[None, :, :]
        diff -= np.rint(diff)
        dist = np.linalg.norm(diff @ cell, axis=2)
        dist[lsymbol_point[:, None] != lsymbol_a[None, :]] = np.inf
        return bool(np.all(dist.min(axis=1) < tol_A))

    frac_b = np.linalg.solve(cell.T, atoms_b.positions.T).T
    if not check_on_sites(frac_b, np.array(atoms_b.get_chemical_symbols())):
        return False
    return all(check_on_sites(frac_a + np.linalg.solve(cell.T, vector), lsymbol_a)
               for vector in atoms_b.cell[:])


def get_translation_count(atoms, tol_A: float = DTOL_TRANSLATION_A) -> int:
    """Count pure lattice translations (identity included) that map a frame onto itself.

    Equals natoms / natoms_primitive. A planar defect pins the cell's period
    normal to the plane, so a relaxed fault frame that gains translations over
    its unrelaxed frame no longer contains the fault.

    Args:
        atoms: Frame to test.
        tol_A: Cartesian tolerance (Å).

    Returns:
        int: Number of vectors ``r_j - r_0`` under which every atom lands on a like atom.
    """
    cell = atoms.cell[:]
    frac = atoms.get_scaled_positions(wrap=False)
    lsymbol = np.array(atoms.get_chemical_symbols())
    mask_species = lsymbol[:, None] != lsymbol[None, :]
    count = 0
    for index in np.where(lsymbol == lsymbol[0])[0]:
        diff = frac[:, None, :] + (frac[index] - frac[0]) - frac[None, :, :]
        diff -= np.rint(diff)
        dist = np.linalg.norm(diff @ cell, axis=2)
        dist[mask_species] = np.inf
        if np.all(dist.min(axis=1) < tol_A):
            count += 1
    return count


def get_frame_action(dict_record, ldefer_model: list = None) -> str:
    """Decide whether one exported frame is the property configuration to label.

    Order matters: a deferred family loses its relaxed frame; a frame that
    the health guards flag never goes to DFT; a relaxed frame that slid back
    to the perfect lattice and an unrelaxed frame that *is* the perfect
    lattice duplicate the reference; a relaxed fault that turned into another
    crystal (fault_lost) or changed its stacking is not the property
    configuration. Genuine unrelaxed frames are always sent.

    Args:
        dict_record: One manifest record written by export_property_structures.
        ldefer_model: Families whose relaxed result is postponed (default LDEFER_MODEL).

    Returns:
        str: ``selected`` or one exclusion key of DICT_ACTION_NOTE.
    """
    ldefer_model = LDEFER_MODEL if ldefer_model is None else ldefer_model
    lflag = dict_record['health']
    state = dict_record['state']
    if dict_record['model'] in ldefer_model and state == 'relaxed':
        return 'deferred'
    if 'below_ground' in lflag or 'dmin_relaxed' in lflag:
        return 'spurious_basin'
    if 'dmin_hard' in lflag:
        return 'overlap'
    if dict_record.get('stable') is False and state == 'relaxed':
        return 'returned_to_lattice'
    if dict_record.get('geometry_status') == 'perfect_lattice':
        return 'perfect_lattice'
    if dict_record.get('geometry_status') == 'fault_lost':
        return 'fault_lost'
    if dict_record.get('geometry_status') == 'stacking_changed':
        return 'stacking_changed'
    return 'selected'


def read_export_records(path_export):
    """Read one export's manifest; refuse failed or pre-health-guard exports."""
    dict_summary = json.loads((path_export / 'summary.json').read_text(encoding='utf-8'))
    if dict_summary['failed']:
        raise ValueError('Export has failed families: ' + str(path_export))
    lrecord = [json.loads(line) for line in
               (path_export / 'manifest.jsonl').read_text(encoding='utf-8').splitlines() if line.strip()]
    if not lrecord or any('health' not in record or 'energy_nnp_eV_atom' not in record
                          for record in lrecord):
        raise ValueError('Export lacks energy/health fields; re-export with current mymetal: '
                         + str(path_export))
    return lrecord


def get_ground_offset_meV(dict_source) -> float:
    """E_nnp/atom above this potential's ground E0 (meV/atom); nan if unknown."""
    energy, energy_ground = dict_source['energy_nnp_eV_atom'], dict_source['energy_ground_eV_atom']
    if energy is None or energy_ground is None:
        return float('nan')
    return 1000. * (energy - energy_ground)


def get_select_table(lunique, lexcluded, dict_identity):
    """Render the human-readable to-do list (pandas-readable, '#' comments)."""
    lline = ['# Property candidates for DFT single-point labelling (deduplicated)',
             '# exports       ' + ' '.join(Path(item['path']).name for item in dict_identity['exports']),
             '# deferred      ' + (' '.join(dict_identity['ldefer_model']) or '-'),
             '# tol_dedup_A   %g' % dict_identity['tol_dedup_A'],
             '# frames %d selected -> %d unique; %d excluded' %
             (sum(len(item['sources']) for item in lunique), len(lunique), len(lexcluded)),
             '# dEg = E_nnp/atom - this potential\'s ground E0 (meV/atom); nan = ideal (no NNP energy)',
             '%-6s %6s %4s %7s %9s  %-4s %-8s %-16s %-12s %s' %
             ('id', 'natoms', 'nsrc', 'dmin_A', 'dEg_meV', 'ph', 'group', 'model', 'state', 'runs')]
    for item in lunique:
        source = item['sources'][0]
        deg = get_ground_offset_meV(source)
        lrun = sorted(set(src['run'] for src in item['sources']))
        lline.append('%-6s %6d %4d %7.3f %9.3f  %-4s %-8s %-16s %-12s %s' %
                     (item['dft_id'], item['natoms'], len(item['sources']), item['min_distance_A'],
                      deg, source['phase'], source['property_group'], source['model'],
                      source['state'] + ('_%03d' % source['frame'] if source['state'] == 'strain' else ''),
                      ','.join(lrun)))
    lline += ['', '# Excluded frames (not sent to DFT)',
              '%-4s %-4s %-8s %-16s %-12s %7s %9s  %s' %
              ('run', 'ph', 'group', 'model', 'state', 'dmin_A', 'dEg_meV', 'action')]
    for source in lexcluded:
        deg = get_ground_offset_meV(source)
        lline.append('%-4s %-4s %-8s %-16s %-12s %7.3f %9.1f  %s' %
                     (source['run'], source['phase'], source['property_group'], source['model'],
                      source['state'], source['min_distance_A'], deg, source['action']))
    return '\n'.join(lline) + '\n'


def select_property_candidates(lpath_export: list = None, path_out: Path = None,
                               ldefer_model: list = None, tag_prefix: str = 'P0',
                               tol_dedup_A: float = DTOL_DEDUP_A) -> dict:
    """Merge per-potential exports into one deduplicated DFT to-do list.

    Applies get_frame_action to every exported frame, then merges frames
    that are the same calculation (check_same_geometry) while keeping all
    sources. Families whose relaxed frame is a spurious basin keep their
    unrelaxed frame as the stand-in. Nothing here writes VASP inputs.

    Args:
        lpath_export: Absolute export directories (one per potential/epoch).
        path_out: Absolute output directory for the to-do list.
        ldefer_model: Families whose relaxed result is postponed (default LDEFER_MODEL).
        tag_prefix: n2p2 tag prefix, e.g. ``P0`` for property round 0.
        tol_dedup_A: Deduplication tolerance (Å).

    Returns:
        dict: Summary with unique/selected/excluded counts per action and run.
    """
    ### check ###
    if not lpath_export:
        raise ValueError('lpath_export must list at least one export directory.')
    for path_value in list(lpath_export) + [path_out]:
        if path_value is None or not Path(path_value).is_absolute():
            raise ValueError('Export/output paths must be absolute: ' + str(path_value))
    lpath_export = [Path(path_value).resolve() for path_value in lpath_export]
    path_out = Path(path_out).resolve()
    ldefer_model = list(LDEFER_MODEL if ldefer_model is None else ldefer_model)
    dict_lrecord = {path_export: read_export_records(path_export) for path_export in lpath_export}
    # to here

    ### prepare ###
    path_out.mkdir(parents=True, exist_ok=True)
    dict_identity = dict(exports=[dict(path=str(path_export),
                                       manifest_sha256=get_file_hash(path_export / 'manifest.jsonl'))
                                  for path_export in lpath_export],
                         ldefer_model=ldefer_model, tag_prefix=tag_prefix, tol_dedup_A=tol_dedup_A)
    path_identity = path_out / 'identity.json'
    if path_identity.exists() and json.loads(path_identity.read_text(encoding='utf-8')) != dict_identity:
        raise ValueError('Selection identity changed; use an independent output directory.')
    path_identity.write_text(json.dumps(dict_identity, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    lunique, lexcluded = [], []
    dict_bucket = {}
    # to here

    ### main ###
    for path_export, lrecord in dict_lrecord.items():
        for record in lrecord:
            dict_source = {key: record.get(key) for key in LSOURCE_KEY}
            dict_source['export'] = str(path_export)
            action = get_frame_action(record, ldefer_model)
            if action != 'selected':
                dict_source['action'] = action
                lexcluded.append(dict_source)
                continue
            path_poscar = Path(record['poscar'])
            if get_file_hash(path_poscar) != record['poscar_sha256']:
                raise ValueError('POSCAR changed after export: ' + str(path_poscar))
            atoms = read(str(path_poscar), format='vasp')
            lbucket = dict_bucket.setdefault(len(atoms), [])
            match = next((item for item in lbucket
                          if check_same_geometry(item['atoms'], atoms, tol_dedup_A)), None)
            if match is not None:
                match['sources'].append(dict_source)
                continue
            item = dict(dft_id='%05d' % len(lunique), atoms=atoms, natoms=len(atoms),
                        tag='-'.join([tag_prefix, record['phase'], record['property_group'], record['model']]),
                        min_distance_A=record['min_distance_A'], poscar=str(path_poscar),
                        poscar_sha256=record['poscar_sha256'], cell_A=atoms.cell[:].tolist(),
                        sources=[dict_source])
            lunique.append(item)
            lbucket.append(item)
    # 伪深谷/层序改变/层错面消失的族靠同族未弛豫帧代替；代替帧缺席说明该族在本轮完全没有标签。
    lfamily_selected = set((src['export'], src['phase'], src['property_group'], src['model'])
                           for item in lunique for src in item['sources']
                           if src['state'] in ('ideal', 'unrelaxed'))
    lmissing = sorted(set((src['run'], src['model']) for src in lexcluded
                          if src['action'] in ('spurious_basin', 'stacking_changed', 'fault_lost')
                          and (src['export'], src['phase'], src['property_group'], src['model'])
                          not in lfamily_selected))
    for run, model in lmissing:
        print('⚠️  ' + run + '/' + model + ': relaxed frame excluded and no unrelaxed stand-in selected')
    (path_out / 'select.jsonl').write_text(
        ''.join(json.dumps({key: value for key, value in item.items() if key != 'atoms'},
                           ensure_ascii=False) + '\n' for item in lunique), encoding='utf-8')
    (path_out / 'excluded.jsonl').write_text(
        ''.join(json.dumps(src, ensure_ascii=False) + '\n' for src in lexcluded), encoding='utf-8')
    (path_out / 'p_post_select.txt').write_text(get_select_table(lunique, lexcluded, dict_identity),
                                                encoding='utf-8')
    dict_action = {}
    for src in lexcluded:
        dict_action[src['action']] = dict_action.get(src['action'], 0) + 1
    dict_run = {}
    for item in lunique:
        for src in item['sources']:
            dict_run[src['run']] = dict_run.get(src['run'], 0) + 1
    dict_natoms = {}
    for item in lunique:
        dict_natoms[item['natoms']] = dict_natoms.get(item['natoms'], 0) + 1
    dict_summary = dict(identity=dict_identity, unique=len(lunique),
                        selected_frames=sum(len(item['sources']) for item in lunique),
                        excluded_frames=len(lexcluded), excluded_by_action=dict_action,
                        selected_frames_by_run=dict_run,
                        unique_by_natoms={str(key): dict_natoms[key] for key in sorted(dict_natoms)},
                        missing_stand_in=['/'.join(pair) for pair in lmissing])
    (path_out / 'summary.json').write_text(json.dumps(dict_summary, ensure_ascii=False, indent=2) + '\n',
                                           encoding='utf-8')
    for action, count in sorted(dict_action.items()):
        print('⚠️  excluded %-20s %4d  (%s)' % (action, count, DICT_ACTION_NOTE[action]))
    print('📊 unique by natoms: ' + ', '.join('%s:%d' % (key, value)
                                              for key, value in dict_summary['unique_by_natoms'].items()))
    print('================ 📊 select_property_candidates: %d exports, %d frames selected -> %d unique '
          'DFT jobs, %d excluded -> %s' % (len(lpath_export), dict_summary['selected_frames'],
                                          len(lunique), len(lexcluded), path_out))
    return dict_summary
    # to here


def find_fault_lost_selected(path_select: Path = None, tol_A: float = DTOL_TRANSLATION_A) -> list:
    """Re-check a frozen to-do list for relaxed faults that turned into another crystal.

    A list selected before the fault_lost guard existed (round 0) cannot be
    re-selected: dft_id would be renumbered under the already prepared VASP
    dirs. This returns the items to drop instead; the unrelaxed stand-in of
    the same family is already in the list (unrelaxed frames are always sent).
    An item is dropped only when every one of its sources is such a frame.

    Args:
        path_select: Absolute select.jsonl (the ``data_dft/<k>/`` copy).
        tol_A: Translation tolerance (Å).

    Returns:
        list: One dict per dropped item with dft_id, action, translation counts
        (relaxed vs each source's unrelaxed frame) and source identities.
    """
    ### check ###
    if path_select is None or not Path(path_select).is_absolute():
        raise ValueError('path_select must be absolute: ' + str(path_select))
    lunique = [json.loads(line) for line in
               Path(path_select).read_text(encoding='utf-8').splitlines() if line.strip()]
    # to here

    ### main ###
    dict_lrecord, ldrop = {}, []
    for item in lunique:
        lsource = item['sources']
        if not all(src['property_group'] in ('ssf', 'seqsf') and src['state'] == 'relaxed'
                   for src in lsource):
            continue
        n_relaxed = get_translation_count(read(item['poscar'], format='vasp'), tol_A)
        lcount = []
        for src in lsource:
            path_export = Path(src['export'])
            if path_export not in dict_lrecord:
                dict_lrecord[path_export] = read_export_records(path_export)
            record_ideal = next((record for record in dict_lrecord[path_export]
                                 if record['model'] == src['model'] and record['state'] == 'ideal'), None)
            lcount.append(None if record_ideal is None else
                          get_translation_count(read(record_ideal['poscar'], format='vasp'), tol_A))
        if all(count is not None and n_relaxed > count for count in lcount):
            ldrop.append(dict(dft_id=item['dft_id'], action='fault_lost', translations_relaxed=n_relaxed,
                              translations_ideal=lcount,
                              sources=[{key: src.get(key) for key in ('run', 'phase', 'property_group',
                                                                      'model', 'state', 'candidate_id')}
                                       for src in lsource]))
            print('⚠️  %s fault_lost: %s relaxed has %d translations vs ideal %s'
                  % (item['dft_id'], '/'.join(src['run'] + ':' + src['model'] for src in lsource),
                     n_relaxed, lcount))
    print('================ 📊 find_fault_lost_selected: %d of %d item(s) to drop (%s)'
          % (len(ldrop), len(lunique), DICT_ACTION_NOTE['fault_lost']))
    return ldrop
    # to here


def get_static_outcar_status(path_outcar: Path = None) -> str:
    """Classify one static (NSW=0) VASP OUTCAR before its labels are used.

    ``pei_vasp_univ_sbatch`` only checks that a static run wrote its energy and
    timing, so a run that exhausted NELM without reaching EDIFF also looks
    finished; here that case is rejected. A static run must also carry exactly
    one force block, otherwise its last geometry is not the POSCAR one.

    Args:
        path_outcar: OUTCAR path.

    Returns:
        str: One key of ``DICT_OUTCAR_STATUS_NOTE`` ('ok' when usable).
    """
    path_outcar = Path(path_outcar)
    if not path_outcar.is_file() or path_outcar.stat().st_size == 0:
        return 'missing'
    text = path_outcar.read_text(encoding='utf-8', errors='ignore')
    if 'General timing and accounting' not in text:
        return 'unfinished'
    if 'aborting loop because EDIFF is reached' not in text:
        return 'ediff_not_reached'
    if text.count('TOTAL-FORCE (eV/Angst)') != 1:
        return 'not_static'
    return 'ok'


def read_static_dft_frame(path_poscar: Path = None, path_outcar: Path = None,
                          tol_pos_A: float = DTOL_OUTCAR_POS_A) -> Atoms:
    """Return the POSCAR geometry labelled with the OUTCAR energy and forces.

    The geometry is taken from POSCAR (full precision, SHA-checked by the
    caller) rather than parsed from OUTCAR: VASP prints lattice vectors as
    fixed-width ``%13.9f``, so a vector >= 100 Å fuses with its neighbour
    (``0.000000000103.734655055``) and ASE's OUTCAR reader fails on such long
    seqSF cells. Energy is the last ``energy(sigma->0)`` (the value ASE reads
    as ``energy``, same as the MD rounds); forces are the last TOTAL-FORCE block,
    whose positions must match POSCAR (minimum image) within ``tol_pos_A``.

    Args:
        path_poscar: POSCAR the static run started from.
        path_outcar: OUTCAR of that run.
        tol_pos_A: Max POSCAR/OUTCAR position mismatch (Å).

    Returns:
        Atoms: POSCAR atoms with a SinglePointCalculator (energy, free_energy, forces).

    Raises:
        ValueError: Missing energy/force block, wrong atom count, or a position mismatch.
    """
    atoms = read(path_poscar, format='vasp')
    natoms = len(atoms)
    lline = Path(path_outcar).read_text(encoding='utf-8', errors='ignore').splitlines()
    lenergy = [float(line.split('energy(sigma->0) =')[1]) for line in lline
               if 'energy(sigma->0) =' in line]
    lfree = [float(line.split()[4]) for line in lline if 'free  energy   TOTEN' in line]
    lforce_head = [index for index, line in enumerate(lline) if 'TOTAL-FORCE (eV/Angst)' in line]
    if not lenergy or not lfree or not lforce_head:
        raise ValueError('No energy/force block in ' + str(path_outcar))
    # 表头下一行是虚线，再往下 natoms 行为 x y z fx fy fz，紧接着又是虚线。
    lrow = lline[lforce_head[-1] + 2:lforce_head[-1] + 2 + natoms + 1]
    if len(lrow) != natoms + 1 or not lrow[-1].strip().startswith('---'):
        raise ValueError('Force block of %s does not hold %d atoms' % (path_outcar, natoms))
    arr = np.array([[float(value) for value in row.split()[:6]] for row in lrow[:-1]])
    dfrac = np.linalg.solve(atoms.cell[:].T, (arr[:, :3] - atoms.get_positions()).T).T
    dpos = np.linalg.norm((dfrac - np.round(dfrac)) @ atoms.cell[:], axis=1).max()
    if dpos > tol_pos_A:
        raise ValueError('OUTCAR positions differ from %s by %.2e A (> %.0e)'
                         % (path_poscar, dpos, tol_pos_A))
    atoms.calc = SinglePointCalculator(atoms, energy=lenergy[-1], free_energy=lfree[-1],
                                       forces=arr[:, 3:])
    return atoms
