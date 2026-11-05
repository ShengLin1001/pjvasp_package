"""-calc lammps for the y_full_relax workflows (E_in_1_2_bulk, NEB); never calls LAMMPS or Slurm.

A fake ``lmp`` on PATH stands in for LAMMPS where the runner is exercised.
"""

import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np
from ase.io import read

from mymetal.build.bulk.stacking import create_stacking_model
from mymetal.build.workflow.E_in_1_2_bulk import SEED_FILE, generate_E_in_1_2_bulk_dirs
from mymetal.build.workflow.lammps_input import (generate_lmp_neb_input,
                                                 generate_lmp_relax_input)
from mymetal.io.lammps import read_lmp_dump_atoms, write_lmp_data, write_lmp_neb_coords
from mymetal.io.vasp import my_read_vasp, my_write_vasp
from mymetal.post.lmp_univ import post_lmp_cases, read_lmp_neb_log
from myvasp import vasp_func as vf

PATH_REPO = Path(__file__).resolve().parents[1]
PATH_RUNNER = PATH_REPO / 'lmp_utils/lmp_universal/pei_lmp_univ_sbatch'
A0 = 2.8597
DICT_POT = {'pair_style': 'eam/alloy', 'pair_coeff': 'Au', 'lmass': ['mass 1 196.9666']}


def make_reference(path_root: Path, calc: str) -> Path:
    """y_full_relax with an ideal 6-layer HCP POSCAR (scale factor A0) and the calc files."""
    path_src = path_root / 'y_full_relax'
    path_src.mkdir(parents=True)
    atoms = create_stacking_model('ABABAB', a=A0, c=A0 * 1.63, lat='hcp')
    my_write_vasp(str(path_src / 'POSCAR'), atoms, lattice_scale_factor=A0, label='ref')
    if calc == 'lammps':
        generate_lmp_relax_input(path_src / 'lmp.in', box_relax='z', **DICT_POT)
    else:
        for name in ['INCAR', 'KPOINTS', 'POTCAR']:
            (path_src / name).write_text(name + '\n')
    return path_src


def assert_same_sites(atoms_a, atoms_b, atol=1e-8):
    """Positions equal modulo the lattice (minimum image)."""
    cell = np.array(atoms_b.get_cell())
    frac = np.linalg.solve(cell.T, (atoms_a.get_positions() - atoms_b.get_positions()).T).T
    frac -= np.round(frac)
    np.testing.assert_allclose(frac @ cell, 0.0, atol=atol)


class TestEIn12Bulk(unittest.TestCase):
    def test_grid_pair_and_cell(self):
        with tempfile.TemporaryDirectory() as tmp:
            path_root = Path(tmp)
            make_reference(path_root, 'vasp')
            lcase = generate_E_in_1_2_bulk_dirs(path_root, [2.6, 3.1], [2.7, 3.0], force=True)
            self.assertEqual(len(lcase), 4)
            lcase = generate_E_in_1_2_bulk_dirs(path_root, [2.6, 3.1], [2.7, 3.0],
                                                if_pair=True, force=True)
            self.assertEqual([p.name for p in lcase],
                             ['a1-2.60000000-a2-2.70000000', 'a1-3.10000000-a2-3.00000000'])
            atoms, a0 = my_read_vasp(str(lcase[0] / 'POSCAR'))
            self.assertAlmostEqual(a0, A0)
            cell = np.array(atoms.get_cell())
            np.testing.assert_allclose(cell[0], [2.6, 0, 0], atol=1e-10)
            np.testing.assert_allclose(cell[1][:2], [-1.3, 2.7 * np.sqrt(3) / 2], atol=1e-10)
            atoms_ref, _ = my_read_vasp(str(path_root / 'y_full_relax/POSCAR'))
            self.assertAlmostEqual(atoms.get_volume(), atoms_ref.get_volume(), places=6)
            for name in ['INCAR', 'KPOINTS', 'POTCAR']:
                self.assertTrue((lcase[0] / name).is_file())

    def test_pair_length_mismatch_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            path_root = Path(tmp)
            make_reference(path_root, 'vasp')
            with self.assertRaises(SystemExit):
                generate_E_in_1_2_bulk_dirs(path_root, [2.6, 2.7], [2.6], if_pair=True)

    def test_lammps_files_and_perturb(self):
        with tempfile.TemporaryDirectory() as tmp:
            path_root = Path(tmp)
            make_reference(path_root, 'lammps')
            lcase = generate_E_in_1_2_bulk_dirs(path_root, [2.8], [2.8], calc='lammps',
                                                perturb=0.01, outdir='y_pert', force=True)
            path_case = lcase[0]
            self.assertTrue((path_case / 'lmp.in').is_file())
            self.assertFalse((path_case / 'INCAR').exists())
            atoms_data = read(path_case / 'data.lmp', format='lammps-data')
            atoms_pos, _ = my_read_vasp(str(path_case / 'POSCAR'))
            assert_same_sites(atoms_data, atoms_pos)   # data.lmp wraps into the box
            self.assertIn('seed', (path_root / 'y_pert' / SEED_FILE).read_text())
            atoms_flat = generate_unperturbed(path_root)
            cell = np.array(atoms_flat.get_cell())
            frac = np.linalg.solve(cell.T, (atoms_pos.get_positions()
                                            - atoms_flat.get_positions()).T).T
            delta = (frac - np.round(frac)) @ cell     # POSCAR wraps a site kicked past 0
            self.assertTrue(np.all(np.abs(delta) <= 0.01 + 1e-12))
            self.assertGreater(np.abs(delta).max(), 0.0)

    def test_lammps_needs_lmp_in(self):
        with tempfile.TemporaryDirectory() as tmp:
            path_root = Path(tmp)
            make_reference(path_root, 'vasp')
            with self.assertRaises(SystemExit):
                generate_E_in_1_2_bulk_dirs(path_root, [2.8], [2.8], calc='lammps')


def generate_unperturbed(path_root: Path):
    """The same single case without perturbation, for the displacement check."""
    lcase = generate_E_in_1_2_bulk_dirs(path_root, [2.8], [2.8], calc='lammps',
                                        outdir='y_flat', force=True)
    atoms, _ = my_read_vasp(str(lcase[0] / 'POSCAR'))
    return atoms


class TestProfileIndex(unittest.TestCase):
    def test_sparse_grid_keeps_every_line_dense_grid_keeps_the_minimum(self):
        from mymetal.universal.plot.workflow import get_profile_index
        self.assertEqual(get_profile_index(8, 3), list(range(8)))
        lidx = get_profile_index(53, 27)
        self.assertIn(27, lidx)
        self.assertEqual((lidx[0], lidx[-1]), (0, 52))
        self.assertLessEqual(len(lidx), 10)


class TestLammpsInput(unittest.TestCase):
    def test_relax_render(self):
        with tempfile.TemporaryDirectory() as tmp:
            text = generate_lmp_relax_input(Path(tmp) / 'lmp.in', box_relax='z',
                                            **DICT_POT).read_text()
            llive = [line for line in text.splitlines() if not line.lstrip().startswith('#')]
            self.assertEqual(sum(line.startswith('read_data') for line in llive), 1)
            self.assertIn('fix br all box/relax z 0.0 fixedpoint 0 0 0', llive)
            self.assertIn('mass 1 196.9666', llive)
            self.assertIn('pair_coeff * * Au', ' '.join(llive))
            self.assertFalse(any('_template' in line for line in llive))
            text = generate_lmp_relax_input(Path(tmp) / 'lmp.in', box_relax='none',
                                            **DICT_POT).read_text()
            self.assertNotIn('box/relax', ' '.join(l for l in text.splitlines()
                                                   if not l.lstrip().startswith('#')))
            with self.assertRaises(SystemExit):
                generate_lmp_relax_input(Path(tmp) / 'lmp.in', box_relax='xyz', **DICT_POT)

    def test_neb_render(self):
        with tempfile.TemporaryDirectory() as tmp:
            text = generate_lmp_neb_input(Path(tmp) / 'lmp.in.neb', spring=5.0, ftol=0.01,
                                          n1=100, n2=200, nevery=10, **DICT_POT).read_text()
            self.assertIn('fix             fneb all neb 5 parallel neigh', text)
            self.assertIn('neb             0.0 0.01 100 200 10 each coords.${i}', text)
            self.assertIn('atom_modify     map array', text)


class TestLammpsIoPost(unittest.TestCase):
    def test_neb_coords_match_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            atoms = create_stacking_model('ABCABC', a=4.08, lat='fcc')
            write_lmp_data(Path(tmp) / 'data.lmp', atoms)
            write_lmp_neb_coords(Path(tmp) / 'coords.1', atoms)
            atoms_data = read(Path(tmp) / 'data.lmp', format='lammps-data')
            llines = (Path(tmp) / 'coords.1').read_text().splitlines()
            self.assertEqual(int(llines[0]), len(atoms))
            lpos = np.array([[float(x) for x in line.split()[1:]] for line in llines[1:]])
            np.testing.assert_allclose(lpos, atoms_data.get_positions(), atol=1e-8)

    def test_dump_origin_removed(self):
        with tempfile.TemporaryDirectory() as tmp:
            path_dump = Path(tmp) / 'y_final.lammpstrj'
            path_dump.write_text(
                'ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\n'
                'ITEM: BOX BOUNDS pp pp pp\n0 3\n0 3\n-0.5 5.5\n'
                'ITEM: ATOMS id type x y z\n1 1 0 0 -0.5\n2 1 1.5 1.5 2.5\n')
            atoms_ref = create_stacking_model('AB', a=3.0, c=6.0, lat='hcp')
            atoms = read_lmp_dump_atoms(path_dump, atoms_ref)
            np.testing.assert_allclose(atoms.get_positions()[:, 2], [0.0, 3.0], atol=1e-12)

    def test_post_cases_vasp_format(self):
        with tempfile.TemporaryDirectory() as tmp:
            path_work = Path(tmp)
            path_case = path_work / 'y_dir' / 'a1-2.80000000-a2-2.80000000'
            path_case.mkdir(parents=True)
            atoms = create_stacking_model('AB', a=2.8, c=4.6, lat='hcp')
            my_write_vasp(str(path_case / 'POSCAR'), atoms, lattice_scale_factor=2.8, label='x')
            (path_case / 'y_lmp_energy.txt').write_text(
                '-7.5 0.0000 in kB 1.0 2.0 3.0 4.0 5.0 6.0\n')
            (path_case / 'y_lmp_done.txt').write_text('pei_lmp_done\n')
            llines = ['ITEM: TIMESTEP', '0', 'ITEM: NUMBER OF ATOMS', '2',
                      'ITEM: BOX BOUNDS xy xz yz pp pp pp']
            cell = atoms.get_cell()
            llines += ['%.10f %.10f %.10f' % (min(0, cell[1][0]), cell[0][0], cell[1][0]),
                       '0 %.10f 0' % cell[1][1], '0 %.10f 0' % cell[2][2],
                       'ITEM: ATOMS id type x y z']
            for i, pos in enumerate(atoms.get_positions(), start=1):
                llines.append('%d 1 %.10f %.10f %.10f' % (i, pos[0], pos[1], pos[2]))
            (path_case / 'y_final.lammpstrj').write_text('\n'.join(llines) + '\n')
            path_unfinished = path_work / 'y_dir' / 'a1-3.00000000-a2-3.00000000'
            path_unfinished.mkdir()
            ldone, lunfinished = post_lmp_cases(path_work, sorted((path_work / 'y_dir').iterdir()))
            self.assertEqual(ldone, [path_case.name])
            self.assertEqual(lunfinished, [path_unfinished.name])
            jobn, etot, _, pres = vf.vasp_read_post_data(str(path_work / 'y_post_data.txt'))
            self.assertEqual(jobn, [path_case.name])
            self.assertAlmostEqual(etot[0], -7.5)
            np.testing.assert_allclose(pres[0], [1, 2, 3, 4, 5, 6])
            _, a0 = my_read_vasp(str(path_case / 'CONTCAR'))
            self.assertAlmostEqual(a0, 2.8)

    def test_neb_log_last_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            header = ('    Step     MaxReplicaForce MaxAtomForce GradV0 GradV1 GradVc EBF EBR RDT '
                      'RD1 PE1 RD2 PE2 RD3 PE3')
            text = '\n'.join([header, '0 1 1 0 0 0 0 0 1 0 -9 0.5 -8 1 -9.5',
                              'Climbing replica = 2', header,
                              '100 0.5 0.2 0 0 0 0 0 1 0 -9 0.5 -8.5 1 -9.5',
                              '200 0.005 0.002 0 0 0 0 0 2 0 -9 0.5 -8.4 1 -9.5',
                              'Loop time of 1.0'])
            (Path(tmp) / 'log.lammps').write_text(text + '\n')
            dict_neb = read_lmp_neb_log(Path(tmp) / 'log.lammps')
            self.assertEqual(dict_neb['step'], 200)
            self.assertEqual(dict_neb['lpe'], [-9.0, -8.4, -9.5])
            self.assertEqual(dict_neb['lrd'], [0.0, 0.5, 1.0])
            self.assertEqual(dict_neb['rdt'], 2.0)


class TestLammpsNebPost(unittest.TestCase):
    def test_force_b_and_vtst_files(self):
        from mymetal.post.lmp_univ import post_lmp_neb
        with tempfile.TemporaryDirectory() as tmp:
            path_ydir = Path(tmp)
            atoms = create_stacking_model('AB', a=2.9, c=4.7, lat='hcp')
            cell = atoms.get_cell()
            lshift = [0.0, 0.5, 1.0]                 # atom 2 slides along x
            lpe = [-7.0, -6.5, -6.9]
            lforce = [[0, 0, 0], [0.3, 0.2, 0], [0, 0, 0]]
            for i, (shift, force) in enumerate(zip(lshift, lforce)):
                path_image = path_ydir / ('%02d' % i)
                path_image.mkdir()
                my_write_vasp(str(path_image / 'POSCAR'), atoms, lattice_scale_factor=2.9)
                lpos = atoms.get_positions().copy()
                lpos[1, 0] += shift
                llines = ['ITEM: TIMESTEP', '0', 'ITEM: NUMBER OF ATOMS', '2',
                          'ITEM: BOX BOUNDS xy xz yz pp pp pp',
                          '%.10f %.10f %.10f' % (cell[1][0], cell[0][0], cell[1][0]),
                          '0 %.10f 0' % cell[1][1], '0 %.10f 0' % cell[2][2],
                          'ITEM: ATOMS id type x y z fx fy fz',
                          '1 1 %.10f %.10f %.10f 0 0 0' % tuple(lpos[0]),
                          '2 1 %.10f %.10f %.10f %g %g %g' % (tuple(lpos[1]) + tuple(force))]
                (path_ydir / ('y_final.%d.lammpstrj' % i)).write_text('\n'.join(llines) + '\n')
            header = ('    Step     MaxReplicaForce MaxAtomForce GradV0 GradV1 GradVc EBF EBR '
                      'RDT RD1 PE1 RD2 PE2 RD3 PE3')
            row = '10 0.005 0.002 0 0 0 0 0 1.0 0 %g 0.5 %g 1 %g' % tuple(lpe)
            (path_ydir / 'log.lammps').write_text(header + '\n' + row + '\n')
            (path_ydir / 'lmp.in').write_text('neb 0.0 0.01 100 100 10 each coords.${i}\n')
            (path_ydir / 'y_lmp_done.txt').write_text('pei_lmp_done\n')
            dict_neb = post_lmp_neb(path_ydir)
            self.assertTrue(dict_neb['converged'])
            self.assertAlmostEqual(dict_neb['barrier_eV'], 0.5)
            np.testing.assert_allclose(dict_neb['ldist'], [0, 0.5, 1.0], atol=1e-8)
            # image 1 is the maximum: tangent mixes to +x, so force_b = +0.3
            self.assertAlmostEqual(dict_neb['lforce_b'][1], 0.3, places=8)
            dat = np.loadtxt(path_ydir / 'neb.dat')
            np.testing.assert_allclose(dat[:, 2], [0, 0.5, 0.1], atol=1e-8)
            dat = np.loadtxt(path_ydir / 'nebef.dat')
            np.testing.assert_allclose(dat[1, 1], 0.2, atol=1e-8)    # perpendicular part
            np.testing.assert_allclose(dat[:, 2], lpe, atol=1e-8)


class TestLammpsRunner(unittest.TestCase):
    def run_runner(self, path_case, body):
        path_bin = path_case.parent / 'bin'
        path_bin.mkdir(exist_ok=True)
        (path_bin / 'lmp').write_text('#!/bin/bash\n' + body + '\n')
        (path_bin / 'lmp').chmod(0o755)
        env = dict(os.environ, MY_LAUNCHER='env', PATH=str(path_bin) + ':' + os.environ['PATH'])
        return subprocess.run(['bash', str(PATH_RUNNER), str(path_case)], env=env,
                              capture_output=True, text=True)

    def test_exit_codes(self):
        with tempfile.TemporaryDirectory() as tmp:
            path_case = Path(tmp) / 'case'
            path_case.mkdir()
            self.assertEqual(self.run_runner(path_case, 'exit 0').returncode, 1)   # no lmp.in
            (path_case / 'lmp.in').write_text('# fake\n')
            self.assertEqual(self.run_runner(path_case, 'exit 0').returncode, 1)   # no marker
            self.assertEqual(self.run_runner(path_case, 'echo done > y_lmp_done.txt').returncode, 0)
            self.assertEqual(self.run_runner(path_case, 'exit 3').returncode, 10)  # skipped


@unittest.skipUnless(shutil.which('nebmake.pl') and shutil.which('pei_vasp_univ_neb_nebmovie'),
                     'VTST nebmake.pl / nebmovie not on PATH')
class TestNebLammps(unittest.TestCase):
    def make_endpoints(self, path_root, c=2.9 * 1.63):
        for name, seq in [('ini', 'ABCABC'), ('fin', 'ABABAB')]:
            (path_root / 'y_full_relax_neb' / name).mkdir(parents=True, exist_ok=True)
            atoms = create_stacking_model(seq, a=2.9, c=c, lat='hcp')
            my_write_vasp(str(path_root / 'y_full_relax_neb' / name / 'CONTCAR'), atoms,
                          lattice_scale_factor=2.9, label=seq)
        generate_lmp_neb_input(path_root / 'y_full_relax_neb/ini/lmp.in.neb', **DICT_POT)

    def test_init_images_reinterpolated(self):
        from mymetal.build.workflow.neb import generate_neb_dirs
        with tempfile.TemporaryDirectory() as tmp:
            path_root = Path(tmp)
            self.make_endpoints(path_root)
            path_src = generate_neb_dirs(path_root, 3, calc='lammps', outdir='y_src', force=True)
            for name in ['01', '02', '03']:   # a bent source chain, as a relaxed NEB would be
                atoms, a0 = my_read_vasp(str(path_src / name / 'POSCAR'))
                atoms.positions[:, 0] += 0.1 * int(name)
                my_write_vasp(str(path_src / name / 'CONTCAR'), atoms, lattice_scale_factor=a0)
            # carried onto a longer c and refined from 3 to 7 images: t = k * 4 / 8
            path_long = path_root / 'long'
            self.make_endpoints(path_long, c=2.9 * 1.70)
            path_new = generate_neb_dirs(path_long, 7, calc='lammps', outdir='y_init',
                                         path_init_images=path_src, force=True)
            atoms_new = {k: my_read_vasp(str(path_new / ('%02d' % k) / 'POSCAR'))[0]
                         for k in range(9)}
            self.assertAlmostEqual(atoms_new[4].get_cell()[2][2], 6 * 2.9 * 1.70 / 2, places=6)
            for k, name in [(2, '01'), (4, '02'), (6, '03')]:   # integer t: source image itself
                atoms_src, _ = my_read_vasp(str(path_src / name / 'CONTCAR'))
                frac = atoms_new[k].get_scaled_positions() - atoms_src.get_scaled_positions()
                np.testing.assert_allclose(frac - np.round(frac), 0.0, atol=1e-8)
            frac_mid = 0.5 * (atoms_new[2].get_scaled_positions()
                              + atoms_new[4].get_scaled_positions())
            frac = atoms_new[3].get_scaled_positions() - frac_mid
            np.testing.assert_allclose(frac - np.round(frac), 0.0, atol=1e-8)
            self.assertEqual(len(list(path_new.glob('coords.*'))), 8)

    def test_lammps_images(self):
        from mymetal.build.workflow.neb import generate_neb_dirs
        with tempfile.TemporaryDirectory() as tmp:
            path_root = Path(tmp)
            self.make_endpoints(path_root)
            path_ydir = generate_neb_dirs(path_root, 3, calc='lammps', force=True)
            self.assertTrue((path_ydir / 'lmp.in').is_file())
            self.assertTrue((path_ydir / 'data.lmp').is_file())
            self.assertEqual(sorted(p.name for p in path_ydir.glob('coords.*')),
                             ['coords.1', 'coords.2', 'coords.3', 'coords.4'])
            llines = (path_ydir / 'coords.4').read_text().splitlines()
            lpos = np.array([[float(x) for x in line.split()[1:]] for line in llines[1:]])
            atoms_fin, _ = my_read_vasp(str(path_ydir / '04' / 'POSCAR'))
            np.testing.assert_allclose(lpos % 1e3, atoms_fin.get_positions() % 1e3, atol=1e-6)



class TestThermalExpansion(unittest.TestCase):
    """pei_vasp_run_thermal_expansion generation + mymetal.post.thermal_expansion on fake output."""

    @staticmethod
    def make_reference(path_root: Path, cubic: bool = True) -> Path:
        from ase.build import bulk
        from mymetal.build.workflow.lammps_input import generate_lmp_npt_input
        path_src = path_root / 'y_full_relax'
        path_src.mkdir(parents=True)
        atoms = bulk('Au', 'fcc', a=4.08, cubic=cubic)
        my_write_vasp(str(path_src / 'POSCAR'), atoms, lattice_scale_factor=1.0, label='ref')
        generate_lmp_npt_input(path_src / 'lmp.in.npt', n_equil=2000, n_prod=10000,
                               **DICT_POT)
        return path_src

    def test_npt_input(self):
        from mymetal.build.workflow.lammps_input import generate_lmp_npt_input
        with tempfile.TemporaryDirectory() as tmp:
            path_in = generate_lmp_npt_input(Path(tmp) / 'lmp.in.npt', **DICT_POT)
            llive = [line for line in path_in.read_text().splitlines()
                     if not line.lstrip().startswith('#')]
            self.assertFalse(any('_template' in line for line in llive))
            text = '\n'.join(llive)
            self.assertIn('run             20000', text)
            self.assertIn('run             100000', text)
            self.assertIn('ave/time 100 10 1000', text)
            self.assertIn('couple xy', text)
            with self.assertRaises(SystemExit):     # partial last block
                generate_lmp_npt_input(Path(tmp) / 'bad', n_prod=1500, **DICT_POT)

    def test_dirs(self):
        from mymetal.build.workflow.thermal_expansion import generate_thermal_expansion_dirs
        from mymetal.post.thermal_expansion import read_case_meta
        with tempfile.TemporaryDirectory() as tmp:
            path_root = Path(tmp)
            self.make_reference(path_root)
            lcase = generate_thermal_expansion_dirs(path_root, [300, 100], [2, 2, 3],
                                                    force=True)
            self.assertEqual([p.name for p in lcase], ['T-0100', 'T-0300'])
            llines = (lcase[1] / 'lmp.in').read_text().splitlines()
            self.assertEqual(llines[1].split(), ['variable', 'temp', 'equal', '300'])
            self.assertEqual(llines[2].split()[:3], ['variable', 'seed', 'equal'])
            self.assertTrue((lcase[1] / 'data.lmp').is_file())
            dict_meta = read_case_meta(lcase[1])
            self.assertEqual(dict_meta, {'T': 300, 'size': [2, 2, 3], 'natoms': 48})
            for bad in ([0], [150.5]):
                with self.assertRaises(SystemExit):
                    generate_thermal_expansion_dirs(path_root, bad, [2, 2, 2], force=True)

    def test_tilted_cell_refused(self):
        from mymetal.build.workflow.thermal_expansion import generate_thermal_expansion_dirs
        with tempfile.TemporaryDirectory() as tmp:
            path_root = Path(tmp)
            self.make_reference(path_root, cubic=False)       # fcc primitive: tilted
            with self.assertRaises(SystemExit):
                generate_thermal_expansion_dirs(path_root, [300], [2, 2, 2], force=True)

    def test_post(self):
        from mymetal.build.workflow.thermal_expansion import generate_thermal_expansion_dirs
        from mymetal.post.thermal_expansion import post_thermal_expansion
        with tempfile.TemporaryDirectory() as tmp:
            path_root = Path(tmp)
            self.make_reference(path_root)
            lcase = generate_thermal_expansion_dirs(path_root, [300, 600], [2, 2, 3],
                                                    force=True)
            path_case = lcase[0]
            # 20 production rows: lx alternates 8.2 +- 0.01 (a1 = lx / 2), lz = 3 * 4.1
            lrow = []
            for i in range(20):
                lx = 8.2 + (0.01 if i % 2 else -0.01)
                lrow.append('%d 300 0 %g %g 12.3 0 0 0 %g -147.0 0.05 0.95 0.01'
                            % (1000 * (i + 3), lx, lx, lx * lx * 12.3))
            (path_case / 'y_lmp_md_prod.txt').write_text(
                '# Time-averaged data\n# TimeStep v_md_t ...\n' + '\n'.join(lrow) + '\n')
            (path_case / 'y_lmp_relax0.txt').write_text('8.16 8.16 12.24 0 0 0 815.0 -150.0\n')
            (path_case / 'log.lammps').write_text(
                '### NNP EW SUMMARY ### TS:        100 EW          7 EWPERSTEP 7.0e-02\n'
                'run             2000\n'
                '### NNP EW SUMMARY ### TS:       2000 EW          5 EWPERSTEP 5.0e-02\n'
                'run             20000\n'
                '### NNP EW SUMMARY ### TS:       3000 EW          2 EWPERSTEP 2.0e-02\n'
                '### NNP EW SUMMARY ### TS:       4000 EW          1 EWPERSTEP 1.0e-02\n')
            (path_case / 'y_lmp_done.txt').write_text('pei_lmp_done\n')
            path_work = path_root / 'y_thermal_expansion'
            df, lunfinished = post_thermal_expansion(path_work, if_plot=False)
            self.assertEqual(lunfinished, ['T-0600'])
            self.assertEqual(df['T_set'].tolist(), [0, 300])
            row0, row = df.iloc[0], df.iloc[1]
            self.assertAlmostEqual(row0['a1'], 4.08)
            self.assertAlmostEqual(row0['a3'], 4.08)
            self.assertAlmostEqual(row0['E_atom'], -150.0 / 48)
            self.assertAlmostEqual(row['a1'], 4.1)
            self.assertAlmostEqual(row['a3'], 4.1)
            self.assertAlmostEqual(row['a1_sem'], 0.0)          # each block = one +- pair
            self.assertAlmostEqual(row['E_atom'], -147.0 / 48)
            self.assertEqual(row['ew_prod'], 3)                  # after the last run only
            self.assertTrue((path_work / 'p_post_thermal_expansion.txt').is_file())



class TestRefOverlay(unittest.TestCase):
    """-ref of pei_vasp_plot_E_in_1_2_bulk.py / pei_vasp_plot_neb.py on synthetic trees."""

    def test_E_in_ref(self):
        import pandas as pd
        from mymetal.post.E_in_1_2_bulk import (get_E_in_1_2_bulk_ref, my_write_E_in_1_2_bulk,
                                                write_E_in_1_2_bulk_ref)
        from mymetal.universal.plot.workflow import (my_plot_E_in_1_2_bulk,
                                                     my_plot_E_in_1_2_bulk_ref)
        def get_df(la, offset):
            lrow = [('a1-%.8f-a2-%.8f' % (a1, a2), '%.8f' % a1, '%.8f' % a2,
                     offset + (a1 - 2.9) ** 2 + 2 * (a2 - 2.9) ** 2, 1.6 + a1 - a2)
                    for a1 in la for a2 in la]
            return pd.DataFrame(lrow, columns=['jobn', 'a1', 'a2', 'Etot_eV_per_atom', 'lca'])
        df = get_df([2.8, 2.9, 3.0], -3.0)
        df.loc[0, 'Etot_eV_per_atom'] += 0.001        # (2.8, 2.8) off by 1 meV/atom
        df_ref = get_df([2.8, 2.9], -4.0)             # coarser grid, other absolute zero
        with tempfile.TemporaryDirectory() as tmp:
            ref_txt = str(Path(tmp) / 'ref.txt')
            my_write_E_in_1_2_bulk(df_ref.assign(Etot_eV=df_ref['Etot_eV_per_atom']), 1,
                                   (2.9, 2.9, -4.0), filename=ref_txt)
            df_cmp, lenergy_ref, lca_ref = get_E_in_1_2_bulk_ref(df, ref_txt)
            self.assertEqual(len(df_cmp), 4)
            np.testing.assert_allclose(df_cmp['dde'], [1.0, 0, 0, 0], atol=1e-6)
            self.assertEqual(int(np.isnan(lenergy_ref).sum()), 5)
            i = df.index[(df.a1 == '2.90000000') & (df.a2 == '2.90000000')][0]
            self.assertAlmostEqual(lenergy_ref[i], -3.0)   # anchored on this tree's scale
            path_txt = Path(tmp) / 'cmp.txt'
            write_E_in_1_2_bulk_ref(df_cmp, ref_txt, str(path_txt))
            self.assertEqual(len(pd.read_csv(path_txt, sep=r'\s+', comment='#')), 4)
            my_plot_E_in_1_2_bulk_ref(df_cmp, if_save=True, savefile=str(Path(tmp) / 'p.pdf'))
            my_plot_E_in_1_2_bulk(df.a1.astype(float).tolist(), df.a2.astype(float).tolist(),
                                  df.Etot_eV_per_atom.tolist(), df.lca.tolist(),
                                  str(Path(tmp) / 'c.pdf'), str(Path(tmp) / 'f.pdf'),
                                  lenergy_ref=lenergy_ref, lca_ref=lca_ref)
            self.assertTrue((Path(tmp) / 'f.pdf').is_file())

    def test_neb_ref(self):
        import pandas as pd
        from mymetal.post.neb import get_neb_ref, my_spline_neb
        from mymetal.universal.plot.workflow import my_plot_neb
        with tempfile.TemporaryDirectory() as tmp:
            path_ref = Path(tmp)
            (path_ref / 'y_dir/00').mkdir(parents=True)
            atoms = create_stacking_model('ABABAB', a=A0, c=A0 * 1.63, lat='hcp')
            my_write_vasp(str(path_ref / 'y_dir/00/POSCAR'), atoms, label='ref')
            (path_ref / 'p_neb.dat').write_text(
                'id dist_cum(Å) energy_rel(eV) force_b(eV/Å) image\n'
                '0 0.0 0.0 0.0 0\n1 1.0 0.6 0.0 1\n2 2.0 0.0 0.0 2\n', encoding='utf-8')
            (path_ref / 'p_nebef.dat').write_text(
                'id force_max(eV/Å) energy(eV) energy_rel(eV)\n'
                '0 0.0 -10.0 0.0\n1 0.0 -9.4 0.6\n2 0.0 -10.0 0.0\n', encoding='utf-8')
            dict_ref = get_neb_ref(str(path_ref))
            self.assertEqual(dict_ref['natoms'], len(atoms))
            nebdf = dict_ref['nebdf']
            nebefdf = dict_ref['nebefdf']
            spline_df = pd.DataFrame({'dist_cum(Å)': [0.0, 2.0], 'energy_rel(eV)': [0.0, 0.0],
                                      'force_b(eV/Å)': [0.0, 0.0]})
            extsdf = pd.DataFrame({'image': [1.0], 'energy_rel(eV)': [0.6]})
            fig, axes = my_plot_neb(nebdf, nebefdf, spline_df, extsdf,
                                    my_spline_neb(nebdf, nebefdf), len(atoms), dict_ref=dict_ref)
            llabel = [t.get_text() for t in axes[0, 0].get_legend().get_texts()]
            self.assertEqual(llabel, ['this, $E_b$ = %.2f' % (600 / len(atoms)),
                                      'ref, $E_b$ = %.2f' % (600 / len(atoms))])


if __name__ == '__main__':
    unittest.main()
