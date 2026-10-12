"""Property candidate selection/DFT-prep tests: no real solvers or scheduler calls."""

import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np
from ase import Atoms
from ase.build import bulk
from ase.io import write

from mymetal.ml.n2p2.workflow_property import PeiN2p2Property
from mymetal.post.property import (check_same_crystal, check_same_geometry, find_fault_lost_selected,
                                   get_file_hash, get_frame_action, get_frame_health,
                                   get_translation_count, select_property_candidates)

CELL_TILT = np.array([[2.9, 0., 0.], [-1.45, 2.51, 0.], [0.3, 0.2, 4.7]])
POSITIONS = np.array([[0., 0., 0.], [1.45, .84, 2.35], [.2, 1.6, 1.1]])


def get_atoms(cell=CELL_TILT, positions=POSITIONS):
    return Atoms('Au%d' % len(positions), positions=positions, cell=cell, pbc=True)


def get_fcc_pair(shift=(0.7, 0.3, 0.)):
    """8-atom FCC (1x1x2 cubic): perfect crystal vs. top half rigidly slipped (a fault)."""
    atoms_perfect = bulk('Au', 'fcc', a=4.08, cubic=True).repeat((1, 1, 2))
    atoms_fault = atoms_perfect.copy()
    mask_top = atoms_fault.positions[:, 2] > atoms_fault.cell[2, 2] / 2 - 0.1
    atoms_fault.positions[mask_top] += shift
    return atoms_perfect, atoms_fault


def write_export(path_export, lframe):
    """Fake export: one POSCAR + manifest record per (atoms, record overrides)."""
    lrecord = []
    for index, (atoms, dict_override) in enumerate(lframe):
        path_poscar = path_export / 'y_dir' / ('%03d' % index) / 'POSCAR'
        path_poscar.parent.mkdir(parents=True)
        write(str(path_poscar), atoms, format='vasp', direct=True, vasp5=True)
        record = dict(round='r', run=path_export.name[:3], epoch=350, seed=1, phase='hcp',
                      property_group='ssf', model='HCP_basal_I2', state='reference', frame=0,
                      strain=None, property_keys=['ssf_HCP_basal_I2'], stable=True,
                      candidate_id='c%03d' % index, natoms=len(atoms), min_distance_A=2.9,
                      energy_nnp_eV_atom=-3.9, energy_ground_eV_atom=-3.92, health=[],
                      geometry_sha256='g', poscar=str(path_poscar),
                      poscar_sha256=get_file_hash(path_poscar))
        record.update(dict_override)
        lrecord.append(record)
    (path_export / 'manifest.jsonl').write_text(''.join(json.dumps(record) + '\n' for record in lrecord),
                                                encoding='utf-8')
    (path_export / 'summary.json').write_text(json.dumps(dict(failed=[])), encoding='utf-8')


class TestPropertySelect(unittest.TestCase):

    def test_frame_health_relaxed_vs_constructed(self):
        # 伪深谷：弛豫帧近距离或低于基态；构造帧 1.95 Å 合法，核重叠才拒绝。
        self.assertEqual(get_frame_health('relaxed', -3.9, 1.88, -3.92), ['dmin_relaxed'])
        self.assertEqual(get_frame_health('relaxed', -5.5, 2.9, -3.92), ['below_ground'])
        self.assertEqual(get_frame_health('relaxed', -3.91, 2.63, -3.92), [])
        self.assertEqual(get_frame_health('ideal', None, 1.95, -3.92), [])
        self.assertEqual(get_frame_health('ideal', None, 1.3, -3.92), ['dmin_hard'])
        self.assertEqual(get_frame_health('strain', -3.9, 2.88, None), [])

    def test_same_geometry_permutation_shift_but_not_strain(self):
        atoms = get_atoms()
        atoms_moved = get_atoms(positions=(POSITIONS + [.7, -.3, 1.9])[[2, 0, 1]])
        self.assertTrue(check_same_geometry(atoms, atoms_moved))
        # 单个原子平移一个格矢 = 同一周期结构。
        atoms_wrapped = get_atoms(positions=POSITIONS + np.array([[0., 0., 0.], CELL_TILT[2], [0., 0., 0.]]))
        self.assertTrue(check_same_geometry(atoms, atoms_wrapped))
        atoms_strain = get_atoms(cell=CELL_TILT * [[1.003], [1.], [1.]])
        atoms_strain.set_scaled_positions(atoms.get_scaled_positions())
        self.assertFalse(check_same_geometry(atoms, atoms_strain))
        positions_fault = POSITIONS.copy()
        positions_fault[1, 0] += .1
        self.assertFalse(check_same_geometry(atoms, get_atoms(positions=positions_fault)))

    def test_same_crystal_lattice_translation_tilt(self):
        # FCC(100) 正交胞：c 倾斜 a/2[110] 是晶格平移 = 完整晶格；倾斜 a/4[110] 才是真错排。
        a = 4.08
        atoms_ref = bulk('Au', 'fcc', a=a, cubic=True).repeat((1, 1, 2))
        for shift, expected in [(.5, True), (.25, False)]:
            atoms_tilt = atoms_ref.copy()
            cell = atoms_ref.cell[:].copy()
            cell[2, :2] += shift * a
            atoms_tilt.set_cell(cell, scale_atoms=False)
            atoms_tilt.wrap()
            self.assertEqual(check_same_crystal(atoms_ref, atoms_tilt), expected)
        # HCP 一个原子挪到另一套子格位置不算同一晶体。
        atoms_hcp = bulk('Au', 'hcp', a=2.88, c=4.8)
        atoms_bad = atoms_hcp.copy()
        atoms_bad.positions[1] += [0., 0., .3]
        self.assertTrue(check_same_crystal(atoms_hcp, atoms_hcp.copy()))
        self.assertFalse(check_same_crystal(atoms_hcp, atoms_bad))

    def test_frame_action_order(self):
        record = dict(model='HCP_basal_I2', state='relaxed', stable=True, health=[])
        self.assertEqual(get_frame_action(record), 'selected')
        # 暂缓机制：默认表 2026-10-11 起为空，这里显式传 round 0 的暂缓族来测机制本身
        ldefer = ['HCP_pyr1n_SF1']
        self.assertEqual(get_frame_action(dict(record, model='HCP_pyr1n_SF1', health=['dmin_relaxed']),
                                          ldefer), 'deferred')
        self.assertEqual(get_frame_action(dict(record, model='HCP_pyr1n_SF1')), 'selected')
        # 暂缓只扣弛豫帧：reference 与未弛豫帧照送（PJ 2026-10-10：未弛豫帧要送）。
        self.assertEqual(get_frame_action(dict(record, model='HCP_pyr1n_SF1', state='ideal',
                                               stable=False, geometry_status='fault')), 'selected')
        self.assertEqual(get_frame_action(dict(record, model='HCP_pyr1n_SF1', state='reference')),
                         'selected')
        self.assertEqual(get_frame_action(dict(record, health=['dmin_relaxed'], stable=False)),
                         'spurious_basin')
        self.assertEqual(get_frame_action(dict(record, stable=False)), 'returned_to_lattice')
        self.assertEqual(get_frame_action(dict(record, state='ideal', stable=False,
                                               geometry_status='perfect_lattice')), 'perfect_lattice')
        self.assertEqual(get_frame_action(dict(record, state='reference', stable=False)), 'selected')
        self.assertEqual(get_frame_action(dict(record, geometry_status='stacking_changed')),
                         'stacking_changed')
        self.assertEqual(get_frame_action(dict(record, geometry_status='fault_lost')), 'fault_lost')
        # 无稳定层错（回到完整晶格）优先于层错面消失：两者都多出平移，前者语义更具体。
        self.assertEqual(get_frame_action(dict(record, stable=False, geometry_status='fault_lost')),
                         'returned_to_lattice')
        # 整胞相变：post_ssf 给 transformed=yes 且 stable=no，归入层错面消失而非回到完整晶格；
        # 未弛豫帧照送。
        self.assertEqual(get_frame_action(dict(record, stable=False, transformed=True)), 'fault_lost')
        self.assertEqual(get_frame_action(dict(record, state='ideal', stable=False, transformed=True,
                                               geometry_status='fault')), 'selected')

    def test_translation_count_fault_lost(self):
        atoms_perfect, atoms_fault = get_fcc_pair()
        # 完美 FCC 是 1 原子原胞：8 原子胞有 8 个纯平移；滑移后只剩 identity 与面内 (a/2, a/2, 0)。
        self.assertEqual(get_translation_count(atoms_perfect), 8)
        self.assertEqual(get_translation_count(atoms_fault), 2)
        # 弛豫噪声在 1e-3 Å 以内仍算同一平移；整体平移不改变计数。
        atoms_noisy = atoms_perfect.copy()
        atoms_noisy.positions += np.random.default_rng(0).uniform(-1e-4, 1e-4, atoms_noisy.positions.shape)
        atoms_noisy.positions += [1e-5, 1e-5, 1e-5]
        self.assertEqual(get_translation_count(atoms_noisy), 8)

    def test_find_fault_lost_selected_and_prepare_skip(self):
        with tempfile.TemporaryDirectory() as tmp:
            path_root = Path(tmp)
            (path_root / 'file').mkdir()
            (path_root / 'file' / 'input.nn.nosf').write_text('number_of_elements 1\nelements Au\n',
                                                              encoding='utf-8')
            path_template = path_root / 'file' / 'vasp_static'
            path_template.mkdir()
            (path_template / 'INCAR').write_text('NSW = 0\nIBRION = -1\n', encoding='utf-8')
            (path_template / 'KPOINTS').write_text('k\n0\nA\n81\n', encoding='utf-8')
            (path_template / 'POTCAR').write_text('potcar\n', encoding='utf-8')
            path_export = path_root / 'data_prop' / '0' / 'candidates_001_000350'
            path_export.mkdir(parents=True)
            atoms_perfect, atoms_fault = get_fcc_pair()
            _, atoms_fault_b = get_fcc_pair(shift=(0.4, 0.9, 0.))
            atoms_fault_b_relaxed = atoms_fault_b.copy()
            atoms_fault_b_relaxed.positions[0] += [0., 0., 0.05]
            # 模拟 round 0：导出时还没有 fault_lost 标记，弛豫后变回完整晶体的帧照常进了清单。
            write_export(path_export, [(atoms_fault, dict(state='ideal')),
                                       (atoms_perfect, dict(state='relaxed')),
                                       (atoms_fault_b, dict(model='HCP_pyr2_SF1', state='ideal')),
                                       (atoms_fault_b_relaxed, dict(model='HCP_pyr2_SF1', state='relaxed'))])
            self.addCleanup(os.chdir, os.getcwd())
            os.chdir(path_root)
            pei = PeiN2p2Property(dir_root=path_root, dir_lmp_utils=path_root,
                                  path_python=Path('/usr/bin/python3'))
            with mock.patch('mymetal.ml.n2p2.workflow_property.os.system', return_value=0) as mock_system:
                pei.select_prop_candidates(iround=0)
                dir_ydir = pei.prepare_prop_dft(iround=0, dir_vasp_template=path_template)
                self.assertEqual(len(list(dir_ydir.iterdir())), 4)
                ldrop = find_fault_lost_selected((path_root / 'data_dft/0/select.jsonl').resolve())
                self.assertEqual([(item['dft_id'], item['translations_relaxed'], item['translations_ideal'])
                                  for item in ldrop], [('00001', 8, [2])])
                path_excluded = pei.exclude_prop_after_select(iround=0)
                self.assertEqual([json.loads(line)['dft_id'] for line in path_excluded.read_text().splitlines()],
                                 ['00001'])
                # 已铺目录保留（不删），但重跑 prepare 不再投它。
                pei.prepare_prop_dft(iround=0, dir_vasp_template=path_template)
                cli_args = mock_system.call_args[0][0]
                self.assertIn('-lsubdir 00000 00002 00003', cli_args)
                self.assertTrue((dir_ydir / '00001' / 'POSCAR').is_file())

    def test_select_dedup_and_stand_in(self):
        with tempfile.TemporaryDirectory() as tmp:
            path_a, path_b = Path(tmp) / '001_000350', Path(tmp) / '002_000350'
            path_a.mkdir()
            path_b.mkdir()
            atoms = get_atoms()
            atoms_dup = get_atoms(positions=(POSITIONS + [.5, .5, .5])[[1, 2, 0]])
            write_export(path_a, [(atoms, {}), (atoms_dup, dict(property_group='cij', model='c11',
                                                                 state='strain', frame=2)),
                                  (atoms, dict(state='relaxed', health=['dmin_relaxed'])),
                                  (get_atoms(positions=POSITIONS * 1.01), dict(state='ideal')),
                                  (atoms, dict(model='HCP_pyr1n_SF1', state='relaxed'))])
            atoms_b = get_atoms(cell=CELL_TILT * 1.001)
            atoms_b.set_scaled_positions(atoms.get_scaled_positions())
            write_export(path_b, [(atoms_b, {})])
            path_out = Path(tmp) / 'y_select'
            dict_summary = select_property_candidates([path_a, path_b], path_out,
                                                      ldefer_model=['HCP_pyr1n_SF1'])
            self.assertEqual(dict_summary['unique'], 3)
            self.assertEqual(dict_summary['selected_frames'], 4)
            self.assertEqual(dict_summary['excluded_by_action'], {'spurious_basin': 1, 'deferred': 1})
            self.assertEqual(dict_summary['missing_stand_in'], [])
            lunique = [json.loads(line) for line in (path_out / 'select.jsonl').read_text().splitlines()]
            self.assertEqual([item['dft_id'] for item in lunique], ['00000', '00001', '00002'])
            self.assertEqual([src['candidate_id'] for src in lunique[0]['sources']], ['c000', 'c001'])
            self.assertEqual(lunique[0]['tag'], 'P0-hcp-ssf-HCP_basal_I2')
            self.assertTrue((path_out / 'p_post_select.txt').is_file())
            # 同一输出目录换了规则 -> 拒绝，避免清单与已铺 DFT 目录错位。
            with self.assertRaises(ValueError):
                select_property_candidates([path_a, path_b], path_out, ldefer_model=[])

    def test_prepare_prop_dft_dirs_and_guards(self):
        with tempfile.TemporaryDirectory() as tmp:
            path_root = Path(tmp)
            (path_root / 'file').mkdir()
            (path_root / 'file' / 'input.nn.nosf').write_text('number_of_elements 1\nelements Au\n',
                                                              encoding='utf-8')
            path_template = path_root / 'file' / 'vasp_static'
            path_template.mkdir()
            (path_template / 'INCAR').write_text('NSW = 0\nIBRION = -1\n', encoding='utf-8')
            (path_template / 'KPOINTS').write_text('k\n0\nA\n81\n', encoding='utf-8')
            (path_template / 'POTCAR').write_text('potcar\n', encoding='utf-8')
            path_export = path_root / 'data_prop' / '0' / 'candidates_001_000350'
            path_export.mkdir(parents=True)
            write_export(path_export, [(get_atoms(), {}),
                                       (get_atoms(positions=POSITIONS * 1.01), dict(state='ideal'))])
            # PeiN2p2 的目录约定相对工作根（驱动脚本先 cd 到根），测试同样先切过去再还原。
            self.addCleanup(os.chdir, os.getcwd())
            os.chdir(path_root)
            pei = PeiN2p2Property(dir_root=path_root, dir_lmp_utils=path_root,
                                  path_python=Path('/usr/bin/python3'))
            with mock.patch('mymetal.ml.n2p2.workflow_property.os.system', return_value=0) as mock_system:
                pei.select_prop_candidates(iround=0)
                dir_ydir = pei.prepare_prop_dft(iround=0, dir_vasp_template=path_template)
                self.assertIn('-if_sbatch False', mock_system.call_args[0][0])
                self.assertEqual(sorted(d.name for d in dir_ydir.iterdir()), ['00000', '00001'])
                lunique = [json.loads(line) for line in
                           (path_root / 'data_dft/0/select.jsonl').read_text().splitlines()]
                self.assertEqual(get_file_hash(dir_ydir / '00001' / 'POSCAR'), lunique[1]['poscar_sha256'])
                # 重跑保留已铺目录；POSCAR 被改则拒绝。
                pei.prepare_prop_dft(iround=0, dir_vasp_template=path_template)
                (dir_ydir / '00000' / 'POSCAR').write_text('changed\n', encoding='utf-8')
                with self.assertRaises(ValueError):
                    pei.prepare_prop_dft(iround=0, dir_vasp_template=path_template)
            (path_template / 'INCAR').write_text('NSW = 50\nIBRION = 2\n', encoding='utf-8')
            with self.assertRaises(ValueError):
                pei.prepare_prop_dft(iround=0, dir_vasp_template=path_template)


if __name__ == '__main__':
    unittest.main()
