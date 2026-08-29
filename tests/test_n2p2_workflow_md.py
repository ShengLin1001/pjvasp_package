"""PeiN2p2MD 里不依赖 Slurm / LAMMPS / VASP 的纯函数单测。

覆盖的都是「错了不会报错、只会安静给出错结果」的部件：三斜盒子换算、dump 解析、
温度区间归属、md_params 渲染、chunks 夹取、静态 INCAR 判据。
"""

import tempfile
import textwrap
import unittest
from pathlib import Path

import numpy as np

from mymetal.ml.n2p2.workflow_md import (PeiN2p2MD, read_lammps_dump, read_temp_ladder,
                                         DICT_MD_PARAMS_DEFAULT)


DUMP_TRICLINIC = textwrap.dedent("""\
    ITEM: TIMESTEP
    100
    ITEM: NUMBER OF ATOMS
    2
    ITEM: BOX BOUNDS xy xz yz pp pp pp
    0.0 2.0 -1.0
    0.0 1.7320508 0.0
    0.0 4.0 0.0
    ITEM: ATOMS id type x y z fx fy fz c_peatom1
    1 1 0.1 0.1 0.1 0.5 -0.5 0.25 -3.5
    2 1 1.0 0.6 2.0 -0.5 0.5 -0.25 -3.25
    ITEM: TIMESTEP
    200
    ITEM: NUMBER OF ATOMS
    2
    ITEM: BOX BOUNDS xy xz yz pp pp pp
    0.0 2.0 -1.0
    0.0 1.7320508 0.0
    0.0 4.0 0.0
    ITEM: ATOMS id type x y z fx fy fz c_peatom1
    1 1 0.2 0.1 0.1 0.1 -0.1 0.0 -3.4
    2 1 1.1 0.6 2.0 -0.1 0.1 0.0 -3.4
    """)


class TestBoundsToCell(unittest.TestCase):
    """LAMMPS dump 的 BOX BOUNDS 是含倾斜量的外接边界，直接相减会把三斜胞算大。"""

    def test_triclinic_bounds_are_unfolded(self):
        from mymetal.ml.n2p2.workflow_md import _bounds_to_cell
        cell = _bounds_to_cell([(0.0, 2.0), (0.0, 1.7320508), (0.0, 4.0)], [-1.0, 0.0, 0.0])
        # xlo = 0 - min(0, -1, 0, -1) = 1 ; xhi = 2 - max(0, -1, 0, -1) = 2 -> lx = 1
        np.testing.assert_allclose(cell[0], [1.0, 0.0, 0.0])
        np.testing.assert_allclose(cell[1], [-1.0, 1.7320508, 0.0])
        np.testing.assert_allclose(cell[2], [0.0, 0.0, 4.0])

    def test_orthogonal_bounds_pass_through(self):
        from mymetal.ml.n2p2.workflow_md import _bounds_to_cell
        cell = _bounds_to_cell([(0.0, 3.0), (0.0, 4.0), (0.0, 5.0)], [0.0, 0.0, 0.0])
        np.testing.assert_allclose(np.diag(cell), [3.0, 4.0, 5.0])


class TestReadLammpsDump(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / 'y_md_hcp.lammpstrj'
        self.path.write_text(DUMP_TRICLINIC, encoding='utf-8')

    def tearDown(self):
        self.tmp.cleanup()

    def test_frames_carry_energy_forces_and_elements(self):
        lframe = read_lammps_dump(self.path, lele=['Au'])
        self.assertEqual([step for step, _ in lframe], [100, 200])
        _, atoms = lframe[0]
        self.assertEqual(atoms.get_chemical_symbols(), ['Au', 'Au'])
        # 总能 = c_peatom1 求和，这是 committee 对比与 input.data 落盘唯一的能量来源
        self.assertAlmostEqual(atoms.get_potential_energy(), -6.75)
        np.testing.assert_allclose(atoms.get_forces()[0], [0.5, -0.5, 0.25])

    def test_missing_energy_column_is_rejected(self):
        # 少了 c_peatom1 时必须报错而不是静默给出 0 能量
        text = DUMP_TRICLINIC.replace(' c_peatom1', '')
        path = Path(self.tmp.name) / 'no_pe.lammpstrj'
        path.write_text(text, encoding='utf-8')
        with self.assertRaises(ValueError):
            read_lammps_dump(path, lele=['Au'])

    def test_unknown_atom_type_is_rejected(self):
        text = DUMP_TRICLINIC.replace('1 1 0.1', '1 3 0.1')
        path = Path(self.tmp.name) / 'bad_type.lammpstrj'
        path.write_text(text, encoding='utf-8')
        with self.assertRaises(ValueError):
            read_lammps_dump(path, lele=['Au'])


class TestTemperatureLadder(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / 'y_md_fcc_temp.txt'
        self.path.write_text('k T_set step_begin step_end\n'
                             '1 25 100 300\n2 50 400 600\n', encoding='utf-8')

    def tearDown(self):
        self.tmp.cleanup()

    def test_ladder_is_parsed_without_header(self):
        self.assertEqual(read_temp_ladder(self.path), [(25.0, 100, 300), (50.0, 400, 600)])

    def test_step_is_assigned_by_interval(self):
        lladder = read_temp_ladder(self.path)
        self.assertEqual(PeiN2p2MD._temperature_of_step(300, lladder), 25.0)
        self.assertEqual(PeiN2p2MD._temperature_of_step(400, lladder), 50.0)
        # 平衡段的 step 落在两档之间，标 nan 而不是硬塞给相邻档
        self.assertTrue(np.isnan(PeiN2p2MD._temperature_of_step(350, lladder)))


class TestRenderMdParams(unittest.TestCase):

    def test_index_and_scalar_variables_are_emitted(self):
        text = PeiN2p2MD._render_md_params(None)
        self.assertIn('variable phase index fcc hcp', text)
        self.assertIn('variable latid index 2 1', text)
        self.assertIn('variable n_temp equal 18', text)

    def test_overrides_win(self):
        text = PeiN2p2MD._render_md_params({'n_temp': 3, 'lphase': ['fcc'], 'llatid': [2],
                                            'laa0': [4.08], 'lnx': [2], 'lny': [2], 'lnz': [2]})
        self.assertIn('variable n_temp equal 3', text)
        self.assertIn('variable phase index fcc', text)

    def test_unequal_per_phase_lists_are_rejected(self):
        # 长度不齐时 md_heating_template.in 的一个 next 推不动全部变量，相循环会错位
        with self.assertRaises(ValueError):
            PeiN2p2MD._render_md_params({'lphase': ['fcc', 'hcp'], 'llatid': [2]})

    def test_default_dict_is_not_mutated(self):
        before = dict(DICT_MD_PARAMS_DEFAULT)
        PeiN2p2MD._render_md_params({'n_temp': 3})
        self.assertEqual(DICT_MD_PARAMS_DEFAULT, before)


class TestClampChunks(unittest.TestCase):
    """提交引擎拒绝 chunks > 作业目录数；补投时只剩几个目录，最需要它跑起来。"""

    def test_chunks_are_clamped_to_job_count(self):
        self.assertEqual(PeiN2p2MD._clamp_chunks({'chunks': 10}, 4)['chunks'], 4)

    def test_chunks_below_job_count_are_untouched(self):
        self.assertEqual(PeiN2p2MD._clamp_chunks({'chunks': 10}, 120)['chunks'], 10)

    def test_original_dict_is_not_mutated(self):
        dict_args = {'chunks': 10}
        PeiN2p2MD._clamp_chunks(dict_args, 2)
        self.assertEqual(dict_args['chunks'], 10)


class TestIncarStaticGuard(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def _write(self, text):
        path = Path(self.tmp.name) / 'INCAR'
        path.write_text(text, encoding='utf-8')
        return path

    def test_static_incar_passes(self):
        PeiN2p2MD._check_incar_static(self._write(' NSW= 0\n IBRION= -1\n ENCUT= 550\n'))

    def test_relaxation_incar_is_rejected(self):
        # MD 帧必须原样标记：带弛豫的模板会让标签对应到另一个构型上
        with self.assertRaises(ValueError):
            PeiN2p2MD._check_incar_static(self._write(' NSW= 150\n IBRION= 2\n'))


class TestInputNnEnergyOnly(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def test_use_short_forces_is_commented_out(self):
        src = Path(self.tmp.name) / 'input.nn'
        dst = Path(self.tmp.name) / 'out.nn'
        src.write_text('epochs 1000\nuse_short_forces   # Use forces.\nrandom_seed 1\n',
                       encoding='utf-8')
        PeiN2p2MD._write_input_nn_energy_only(src, dst)
        ltext = dst.read_text(encoding='utf-8').splitlines()
        self.assertTrue(any(l.startswith('#use_short_forces') for l in ltext))
        self.assertIn('epochs 1000', ltext)
        # 只注释目标关键词，其它行原样保留
        self.assertIn('random_seed 1', ltext)


class TestCountStructures(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def test_begin_lines_are_counted(self):
        path = Path(self.tmp.name) / 'input.data'
        path.write_text('begin set=train\nenergy 1.0\nend\nbegin\nenergy 2.0\nend\n',
                        encoding='utf-8')
        self.assertEqual(PeiN2p2MD._count_structures(path), 2)

    def test_missing_file_counts_zero(self):
        self.assertEqual(PeiN2p2MD._count_structures(Path(self.tmp.name) / 'nope.data'), 0)


class TestPotentialId(unittest.TestCase):

    def test_id_joins_run_and_epoch(self):
        self.assertEqual(PeiN2p2MD._potential_id('001', '001000'), '001_001000')


if __name__ == '__main__':
    unittest.main()
