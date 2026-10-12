"""Property DFT collect/merge tests: fake OUTCARs, no VASP or scheduler calls."""

import contextlib
import io
import json
import os
from pathlib import Path
import tempfile
import unittest

import numpy as np
from ase import Atoms
from ase.build import bulk
from ase.io import write

from mymetal.ml.n2p2.workflow_property import PeiN2p2Property
from mymetal.post.property import get_file_hash, get_static_outcar_status, read_static_dft_frame

LINE_DASH = ' ' + '-' * 83


def write_outcar(path_outcar, atoms, energy, forces, if_finished=True, if_ediff=True,
                 n_force_block=1, shift_A=0.):
    """Minimal static OUTCAR: lattice print (fixed %13.9f like VASP), force block(s), energies."""
    lline = ['      direct lattice vectors                 reciprocal lattice vectors']
    lline += [''.join('%13.9f' % value for value in row) + '     0.000000000  0.000000000  0.000000000'
              for row in atoms.cell[:]]
    if if_ediff:
        lline.append('------------------------ aborting loop because EDIFF is reached ' + '-' * 40)
    for _ in range(n_force_block):
        lline += [' POSITION                                       TOTAL-FORCE (eV/Angst)', LINE_DASH]
        lline += ['%13.5f%13.5f%13.5f%15.6f%14.6f%14.6f' % (*position, *force)
                  for position, force in zip(atoms.get_positions() + shift_A, forces)]
        lline.append(LINE_DASH)
    lline += ['  free  energy   TOTEN  =     %16.8f eV' % (energy + 0.003), '',
              '  energy  without entropy=     %16.8f  energy(sigma->0) =     %16.8f' % (energy + 0.006, energy)]
    if if_finished:
        lline.append(' General timing and accounting informations for this job:')
    Path(path_outcar).write_text('\n'.join(lline) + '\n', encoding='utf-8')


def get_forces(natoms, seed=0):
    return np.round(np.random.default_rng(seed).normal(0., 0.05, (natoms, 3)), 6)


def get_root(tmp):
    """Work root with the files PeiN2p2Property needs, cwd switched to it."""
    path_root = Path(tmp)
    (path_root / 'file').mkdir()
    (path_root / 'file' / 'input.nn.nosf').write_text('number_of_elements 1\nelements Au\n', encoding='utf-8')
    os.chdir(path_root)
    return path_root, PeiN2p2Property(dir_root=path_root, dir_lmp_utils=path_root,
                                      path_python=Path('/usr/bin/python3'))


def write_round(path_root, dict_case, lexcluded=()):
    """data_dft/0 with select.jsonl, excluded list and y_dir/<id>/{POSCAR,OUTCAR}."""
    dir_round = path_root / 'data_dft' / '0'
    lunique = []
    for dft_id, (atoms, tag, dict_outcar) in dict_case.items():
        dir_job = dir_round / 'y_dir' / dft_id
        dir_job.mkdir(parents=True)
        write(str(dir_job / 'POSCAR'), atoms, format='vasp', direct=True, vasp5=True)
        if dict_outcar is not None:
            write_outcar(dir_job / 'OUTCAR', atoms, **dict_outcar)
        lunique.append(dict(dft_id=dft_id, tag=tag, poscar=str(dir_job / 'POSCAR'),
                            poscar_sha256=get_file_hash(dir_job / 'POSCAR'),
                            sources=[dict(run='001', model=tag.split('-')[-1], state='strain')]))
    (dir_round / 'select.jsonl').write_text(''.join(json.dumps(item) + '\n' for item in lunique),
                                            encoding='utf-8')
    (dir_round / 'excluded_after_select.jsonl').write_text(
        ''.join(json.dumps(dict(dft_id=dft_id, action='fault_lost')) + '\n' for dft_id in lexcluded),
        encoding='utf-8')
    return dir_round / 'y_dir'


def get_long_cell():
    """2-atom seqSF-like column with c > 100 Å: VASP fuses the lattice print, ASE cannot parse it."""
    return Atoms('Au2', positions=[[0., 0., 0.6], [1.4427, 0.83296, 2.96361]],
                 cell=[[2.885362799, 0., 0.], [-1.4426814, 2.498797483, 0.], [0., 0., 103.734655055]],
                 pbc=True)


class TestPropertyCollect(unittest.TestCase):

    def setUp(self):
        self.addCleanup(os.chdir, os.getcwd())

    def test_static_outcar_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            atoms, path_outcar = bulk('Au', 'fcc', a=4.08, cubic=True), Path(tmp) / 'OUTCAR'
            self.assertEqual(get_static_outcar_status(path_outcar), 'missing')
            dict_case = {'unfinished': dict(if_finished=False), 'ediff_not_reached': dict(if_ediff=False),
                         'not_static': dict(n_force_block=2), 'ok': {}}
            for status, dict_kw in dict_case.items():
                write_outcar(path_outcar, atoms, -15.7, get_forces(4), **dict_kw)
                self.assertEqual(get_static_outcar_status(path_outcar), status)

    def test_read_static_dft_frame_long_cell_and_mismatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            atoms, forces = get_long_cell(), get_forces(2)
            path_poscar, path_outcar = Path(tmp) / 'POSCAR', Path(tmp) / 'OUTCAR'
            write(str(path_poscar), atoms, format='vasp', direct=True, vasp5=True)
            write_outcar(path_outcar, atoms, -7.8, forces)
            self.assertIn('0.000000000103.734655055', path_outcar.read_text())
            atoms_dft = read_static_dft_frame(path_poscar, path_outcar)
            self.assertAlmostEqual(atoms_dft.get_potential_energy(), -7.8, places=8)
            self.assertAlmostEqual(atoms_dft.calc.get_property('free_energy', atoms_dft), -7.797, places=8)
            np.testing.assert_allclose(atoms_dft.get_forces(), forces, atol=1e-12)
            np.testing.assert_allclose(atoms_dft.cell[:], atoms.cell[:], atol=1e-9)
            write_outcar(path_outcar, atoms, -7.8, forces, shift_A=0.01)
            with self.assertRaises(ValueError):
                read_static_dft_frame(path_poscar, path_outcar)

    def test_collect_resize_tag_exclude_and_valve(self):
        with tempfile.TemporaryDirectory() as tmp:
            path_root, pei = get_root(tmp)
            atoms_fcc = bulk('Au', 'fcc', a=4.08, cubic=True)
            atoms_hcp = bulk('Au', 'hcp', a=2.88, c=4.75)
            forces_fcc = get_forces(4, seed=1)
            dict_case = {'00000': (atoms_fcc, 'P0-fcc-cij-c11', dict(energy=-15.7, forces=forces_fcc)),
                         '00001': (atoms_hcp, 'P0-hcp-stretch-bulk', dict(energy=-7.8, forces=get_forces(2))),
                         '00002': (get_long_cell(), 'P0-fcc-seqsf-FCC_111_ESF', None),
                         '00003': (get_long_cell(), 'P0-fcc-seqsf-FCC_111_TB',
                                   dict(energy=-7.81, forces=get_forces(2), if_ediff=False))}
            dir_ydir = write_round(path_root, dict_case, lexcluded=['00002'])
            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                path_dft = pei.collect_prop_dft(iround=0)
            self.assertIn('00002: excluded after select', stdout.getvalue())
            self.assertIn('00003: ediff_not_reached', stdout.getvalue())
            text = path_dft.read_text()
            self.assertEqual(text.count('begin set=train'), 2)
            self.assertIn('tag=P0-fcc-cij-c11 | frame=1/1 | file=dft_prop0_00000', text)
            lrecord = [json.loads(line) for line in
                       (dir_ydir / 'input.data.dft.manifest.jsonl').read_text().splitlines()]
            self.assertEqual([(record['dft_id'], record['natoms_dft'], record['natoms_data'], record['rep'])
                              for record in lrecord], [('00000', 4, 48, [3, 4, 1]), ('00001', 2, 48, [4, 6, 1])])
            self.assertEqual(lrecord[0]['sources'][0]['model'], 'c11')
            # 扩胞：E×M、力按 ASE repeat 的块序平铺。
            lblock = text.split('begin')[1:]
            self.assertAlmostEqual(float(lblock[0].split('energy')[-1].split()[0]), -15.7 * 12, places=6)
            arr_force = np.array([[float(value) for value in line.split()[7:10]]
                                  for line in lblock[0].splitlines() if line.startswith('atom')])
            np.testing.assert_allclose(arr_force, np.tile(forces_fcc, (12, 1)), atol=1e-6)
            # 补算好 00003 后阀门不跳过（条数不等），收齐后再调才跳过。
            write_outcar(dir_ydir / '00003' / 'OUTCAR', get_long_cell(), -7.81, get_forces(2))
            pei.collect_prop_dft(iround=0)
            self.assertEqual(path_dft.read_text().count('begin set=train'), 3)
            self.assertIn('tag=P0-fcc-seqsf-FCC_111_TB', path_dft.read_text())
            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                pei.collect_prop_dft(iround=0)
            self.assertIn('skip collect_prop_dft', stdout.getvalue())

    def test_collect_poscar_sha_mismatch_raises(self):
        with tempfile.TemporaryDirectory() as tmp:
            path_root, pei = get_root(tmp)
            atoms = bulk('Au', 'fcc', a=4.08, cubic=True)
            dir_ydir = write_round(path_root, {'00000': (atoms, 'P0-fcc-stretch-bulk',
                                                         dict(energy=-15.7, forces=get_forces(4)))})
            atoms.positions[0] += 0.001
            write(str(dir_ydir / '00000' / 'POSCAR'), atoms, format='vasp', direct=True, vasp5=True)
            with self.assertRaises(ValueError):
                pei.collect_prop_dft(iround=0)
            self.assertFalse((dir_ydir / 'input.data.dft').exists())

    def test_append_prop_to_train(self):
        with tempfile.TemporaryDirectory() as tmp:
            path_root, pei = get_root(tmp)
            atoms = bulk('Au', 'fcc', a=4.08, cubic=True)
            dict_case = {'00000': (atoms, 'P0-fcc-stretch-bulk', dict(energy=-15.7, forces=get_forces(4))),
                         '00001': (atoms, 'P0-fcc-cij-c44', dict(energy=-15.6, forces=get_forces(4),
                                                                   if_finished=False))}
            write_round(path_root, dict_case)
            path_src = path_root / 'stage2_input.data'
            path_src.write_text('begin set=train\ncomment | tag=A12-1\nend\n'
                                'begin set=test\ncomment | tag=A31-1\nend\n', encoding='utf-8')
            pei.collect_prop_dft(iround=0)
            with self.assertRaises(ValueError):
                pei.append_prop_to_train(liround=[0], path_data_ini_src=path_src)
            write_outcar(path_root / 'data_dft/0/y_dir/00001/OUTCAR', atoms, -15.6, get_forces(4))
            pei.collect_prop_dft(iround=0)
            with self.assertRaises(ValueError):
                pei.append_prop_to_train(liround=[0], path_data_ini_src=path_src, sha256_ini_src='0' * 64)
            path_train = pei.append_prop_to_train(liround=[0], path_data_ini_src=path_src,
                                                  sha256_ini_src=get_file_hash(path_src))
            self.assertEqual(pei._count_sample_types(path_train), {'train': 3, 'test': 1, 'unset': 0})
            self.assertTrue(path_train.read_text().startswith(path_src.read_text()))
            dict_source = json.loads((path_root / 'data/train/input.data.ini.source.json').read_text())
            self.assertEqual(dict_source['sha256'], get_file_hash(path_src))
            # 重跑幂等：从快照重建，不会把 P0 写两遍。
            sha_train = get_file_hash(path_train)
            pei.append_prop_to_train(liround=[0])
            self.assertEqual(get_file_hash(path_train), sha_train)
            # 快照一经建立，来源变了就拒绝。
            path_src.write_text('begin set=train\nend\n', encoding='utf-8')
            with self.assertRaises(ValueError):
                pei.append_prop_to_train(liround=[0], path_data_ini_src=path_src)


if __name__ == '__main__':
    unittest.main()
