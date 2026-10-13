"""Planar-defect post tests: whole-cell transformation criteria and y_post_ssf.txt read-back."""

from pathlib import Path
import tempfile
import unittest

import numpy as np
from ase import Atoms
from ase.build import bulk

from mymetal.post.planar import get_bulk_deviation, get_layer_gap_change, read_ssf, write_ssf_output

# tilted a3（带面内分量），与 SSF 的 tilted cell 同形；4 层波纹：窄 0.8 / 宽 1.6 Å 交替
CELL = np.array([[2.9, 0., 0.], [0., 2.5, 0.], [0.4, 0.3, 4.8]])
LZ_IDEAL = [0.0, 0.8, 2.4, 3.2]


def get_layers(lz=LZ_IDEAL, cell=CELL):
    """每层一个原子，z 给定，面内放在 a3 方向上（分数坐标 (0, 0, z/c_z)）。"""
    lscaled = [[0.1 * i, 0.2 * i, z / cell[2, 2]] for i, z in enumerate(lz)]
    atoms = Atoms('Au%d' % len(lz), cell=cell, pbc=True)
    atoms.set_scaled_positions(lscaled)
    return atoms


class TestLayerGapChange(unittest.TestCase):
    def test_local_relaxation_is_small(self):
        # 只动一层（层错面附近弛豫）：4 个层间距里 2 个变，中位数 = 0.05
        atoms_relaxed = get_layers([0.0, 0.8, 2.5, 3.2])
        self.assertAlmostEqual(get_layer_gap_change(get_layers(), atoms_relaxed), 0.05, places=6)

    def test_whole_cell_flattened(self):
        # 波纹层全部压平成等间距 2.4 Å 的平层（prism1w SF2 的整胞相变形态）：每个间距都变 0.8
        cell = CELL.copy()
        cell[2, 2] = 9.6
        atoms_relaxed = get_layers([0.0, 2.4, 4.8, 7.2], cell=cell)
        self.assertGreater(get_layer_gap_change(get_layers(), atoms_relaxed), 0.2)

    def test_wrap_across_boundary_is_not_a_change(self):
        # 末态有原子折回到 a3 另一侧（分数坐标差 1）：不算层间距变化
        atoms_relaxed = get_layers()
        lscaled = atoms_relaxed.get_scaled_positions(wrap=False)
        lscaled[3, 2] -= 1.0
        atoms_relaxed.set_scaled_positions(lscaled)
        self.assertAlmostEqual(get_layer_gap_change(get_layers(), atoms_relaxed), 0.0, places=6)


class TestBulkDeviation(unittest.TestCase):
    def get_hcp(self):
        return bulk('Au', 'hcp', a=2.86, c=4.80).repeat((3, 3, 4))

    def test_identical_is_zero(self):
        atoms_ref = self.get_hcp()
        self.assertAlmostEqual(get_bulk_deviation(atoms_ref, atoms_ref.copy()), 0.0, places=8)

    def test_local_distortion_keeps_bulk_atoms(self):
        # 只推动一个原子（局域层错的弛豫场）：其余原子仍是体相环境
        atoms_ref = self.get_hcp()
        atoms_relaxed = atoms_ref.copy()
        atoms_relaxed.positions[0] += [0.2, 0.0, 0.0]
        self.assertLess(get_bulk_deviation(atoms_ref, atoms_relaxed), 1e-8)

    def test_whole_cell_change_is_flagged(self):
        # 整胞沿 c 压 10%（层间距处处改变，类比整胞相变 / 薄胞）：没有一个原子像体相
        atoms_ref = self.get_hcp()
        atoms_relaxed = atoms_ref.copy()
        cell = atoms_relaxed.cell[:].copy()
        cell[2] *= 0.9
        atoms_relaxed.set_cell(cell, scale_atoms=True)
        self.assertGreater(get_bulk_deviation(atoms_ref, atoms_relaxed), 0.05)


class TestSsfCell(unittest.TestCase):
    def test_only_stable_keeps_gamma(self):
        from mymetal.ml.n2p2.workflow import PeiN2p2
        self.assertEqual(PeiN2p2._get_ssf_cell({'gamma': 300.0, 'stable': True, 'transformed': False}),
                         (300.0, 1))
        value, state = PeiN2p2._get_ssf_cell({'gamma': 0.01, 'stable': False, 'transformed': False})
        self.assertTrue(np.isnan(value))
        self.assertEqual(state, 0)
        # 整胞相变也是 stable=False，但状态码要标成 transformed
        value, state = PeiN2p2._get_ssf_cell({'gamma': 330.0, 'stable': False, 'transformed': True})
        self.assertTrue(np.isnan(value))
        self.assertEqual(state, -1)


class TestReadSsf(unittest.TestCase):
    def get_dict_out(self, stable, transformed):
        return dict(gamma=323.75, dE=1.2, stable=stable, transformed=transformed, gap_change=0.81, bulk_dev=0.5,
                    Asf=8.0, a11=2.9, a22=2.5, E0bulk=-3.9, nfault=1,
                    t1=0.5, t2=0.5, dt1=0.0, dt2=0.0, da33=0.0, umax=0.1)

    def test_round_trip_transformed(self):
        with tempfile.TemporaryDirectory() as tmp:
            path_txt = Path(tmp) / 'y_post_ssf.txt'
            write_ssf_output(self.get_dict_out(stable=False, transformed=True), save_txt_path=path_txt)
            dict_read = read_ssf(path_txt)
        self.assertFalse(dict_read['stable'])
        self.assertTrue(dict_read['transformed'])
        self.assertAlmostEqual(dict_read['gamma'], 323.75, places=6)
        # property.py 要求返回值全部有限，否则整族判失败
        self.assertTrue(all(np.isfinite(value) for value in dict_read.values()))

    def test_old_txt_without_transformed_line(self):
        with tempfile.TemporaryDirectory() as tmp:
            path_txt = Path(tmp) / 'y_post_ssf.txt'
            write_ssf_output(self.get_dict_out(stable=True, transformed=False), save_txt_path=path_txt)
            lline = [line for line in path_txt.read_text(encoding='utf-8').splitlines()
                     if not line.strip().startswith(('transformed', 'gap_change', 'bulk_dev'))]
            path_txt.write_text('\n'.join(lline) + '\n', encoding='utf-8')
            dict_read = read_ssf(path_txt)
        self.assertTrue(dict_read['stable'])
        self.assertFalse(dict_read['transformed'])
        self.assertTrue(all(np.isfinite(value) for value in dict_read.values()))


class TestReadDftReferencePlanar(unittest.TestCase):
    def test_planar_groups_keep_only_values(self):
        import json
        from mymetal.ml.n2p2.dataset import read_dft_reference
        dict_json = {'meta': {'status': 'final', 'open_items': []},
                     'surface': {'surf_rel_FCC_111': {'value': 973.7}},
                     'ssf': {'ssf_FCC_111_ISF': {'value': 35.6},
                             'ssf_HCP_prism1w_SF2': {'value': None, 'not_converged_reason': 'no_stable_run'}},
                     'seqsf': {'ssf_FCC_111_TB': {'value': 14.88}}}
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / 'A5x_reference.json').write_text(json.dumps(dict_json), encoding='utf-8')
            dict_ref = read_dft_reference(tmp, verbose=False)
            dict_skip = read_dft_reference(tmp, planar_ref=None, verbose=False)
        self.assertEqual(dict_ref['surface'], {'surf_rel_FCC_111': 973.7})
        # 无稳定层错（value=None）不给参考线
        self.assertEqual(dict_ref['ssf'], {'ssf_FCC_111_ISF': 35.6})
        self.assertEqual(dict_ref['seqsf'], {'ssf_FCC_111_TB': 14.88})
        self.assertNotIn('ssf', dict_skip)


if __name__ == '__main__':
    unittest.main()
