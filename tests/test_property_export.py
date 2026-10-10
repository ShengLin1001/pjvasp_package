"""Property export regression tests: no real solvers or scheduler calls."""

import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from ase import Atoms
from ase.io import write

from mymetal.post.property import (check_poscar_roundtrip, generate_candidate,
                                   get_property_frame_rules, get_stacking_letters,
                                   get_periodic_fault_count, read_property_data,
                                   read_property_dump, read_rule_frames)
from mymetal.build.bulk.stacking import create_stacking_model, get_stacking_sequence


class TestPropertyExport(unittest.TestCase):
    def write_dump(self, path_dump, lframe):
        """Synthetic metal-unit dump with bounds computed by the LAMMPS formula."""
        ltext = []
        for cell, origin, positions in lframe:
            lx, ly, lz = np.diag(cell)
            xy, xz, yz = cell[1, 0], cell[2, 0], cell[2, 1]
            xlo, ylo, zlo = origin
            ltext += ['ITEM: TIMESTEP', '0', 'ITEM: NUMBER OF ATOMS', str(len(positions)),
                      'ITEM: BOX BOUNDS xy xz yz pp pp pp',
                      '%s %s %s' % (xlo + min(0., xy, xz, xy + xz),
                                    xlo + lx + max(0., xy, xz, xy + xz), xy),
                      '%s %s %s' % (ylo + min(0., yz), ylo + ly + max(0., yz), xz),
                      '%s %s %s' % (zlo, zlo + lz, yz), 'ITEM: ATOMS id type x y z']
            # Deliberately write reversed atom-ID order.
            for index in reversed(range(len(positions))):
                pos = positions[index] + origin
                ltext.append('%d 1 %.12f %.12f %.12f' % (index + 1, *pos))
        path_dump.write_text('\n'.join(ltext) + '\n', encoding='utf-8')

    def test_opposite_tilt_and_origin(self):
        with tempfile.TemporaryDirectory() as tmp:
            path_root = Path(tmp)
            cell = np.array([[4., 0., 0.], [-2., 5., 0.], [3., -1., 8.]])
            origin = np.array([-3., 2., -4.])
            positions = np.array([[.1, .2, .3], [2., 3., 4.]])
            self.write_dump(path_root / 'movie', [(cell, origin, positions)])
            atoms = read_property_dump(path_root / 'movie', ['Au'])[0]
            np.testing.assert_allclose(atoms.cell, cell)
            np.testing.assert_allclose(atoms.positions, positions)
            self.assertEqual(atoms.info['box_origin_A'], list(origin))
            self.assertEqual(atoms.get_chemical_symbols(), ['Au', 'Au'])
            self.assertEqual(atoms.info['atom_ids'], [1, 2])

    def test_data_origin_and_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            path_data = Path(tmp) / 'data'
            path_data.write_text('LAMMPS\n\n2 atoms\n1 atom types\n\n'
                                 '-3 1 xlo xhi\n2 7 ylo yhi\n-4 4 zlo zhi\n'
                                 '-2 3 -1 xy xz yz\n\nAtoms # atomic\n\n'
                                 '2 1 -1 5 0\n1 1 -2.9 2.2 -3.7\n', encoding='utf-8')
            atoms = read_property_data(path_data, ['Au'])
            np.testing.assert_allclose(atoms.positions, [[.1, .2, .3], [2., 3., 4.]])
            self.assertEqual(atoms.info['box_origin_A'], [-3., 2., -4.])
            self.assertEqual(atoms.info['atom_ids'], [1, 2])

    def test_explicit_element_mapping(self):
        with self.assertRaises(ValueError):
            read_property_dump(Path('/unused'))

    def test_rules_complete_and_not_last_frame(self):
        lrule = get_property_frame_rules()
        self.assertEqual(len(lrule), 31)
        self.assertEqual(len([rule for rule in lrule if rule['group'] == 'cij']), 10)
        lsurface = [rule for rule in lrule if rule['group'] == 'surface']
        self.assertEqual(len(lsurface), 6)
        self.assertTrue(all(rule['states'][1] == 'unrelaxed' for rule in lsurface))
        lesf = [rule for rule in lrule if rule['model'] == 'FCC_111_ESF']
        self.assertNotEqual(len(lesf[0]['sequence_reference']), len(lesf[0]['sequence_fault']))
        self.assertEqual(lesf[0]['nfault'], 2)

    def test_poscar_geometry_and_idempotence(self):
        with tempfile.TemporaryDirectory() as tmp:
            path_root = Path(tmp)
            atoms = Atoms('Au2', positions=[[.1, .2, .3], [2., 3., 4.]],
                          cell=[[4., 0., 0.], [-2., 5., 0.], [3., -1., 8.]], pbc=True)
            dict_meta = dict(source_sha256='same')
            generate_candidate(atoms, path_root, 'test', dict_meta.copy())
            generate_candidate(atoms, path_root, 'test', dict_meta.copy())
            path_poscar = path_root / 'y_dir/test/POSCAR'
            self.assertLess(check_poscar_roundtrip(atoms, path_poscar)['cell_error_A'], 1e-9)
            with self.assertRaises(ValueError):
                generate_candidate(atoms, path_root, 'test', dict(source_sha256='changed'))

    def test_surface_wrap_not_relaxation(self):
        with tempfile.TemporaryDirectory() as tmp:
            path_root = Path(tmp)
            cell_ref = np.diag([4., 5., 8.])
            cell_slab = np.diag([4., 5., 28.])
            positions = np.array([[.1, .2, .3], [2., 3., 4.]])
            positions_wrapped = positions.copy()
            positions_wrapped[0] += cell_ref[0] + cell_ref[1]
            self.write_dump(path_root / 'movie', [(cell_ref, np.zeros(3), positions),
                                                  (cell_slab, np.zeros(3), positions_wrapped),
                                                  (cell_slab, np.zeros(3), positions)])
            (path_root / 'y_post_data.txt').write_text('header\n0 0\n1 1\n2 .5\n', encoding='utf-8')
            dict_rule = dict(group='surface', source='movie', format='lammps-dump-text',
                             states=['reference', 'unrelaxed', 'relaxed'])
            self.assertEqual(len(read_rule_frames(path_root, dict_rule, ['Au'])), 3)
            positions_wrapped[0, 2] += .1
            self.write_dump(path_root / 'movie', [(cell_ref, np.zeros(3), positions),
                                                  (cell_slab, np.zeros(3), positions_wrapped),
                                                  (cell_slab, np.zeros(3), positions)])
            with self.assertRaises(ValueError):
                read_rule_frames(path_root, dict_rule, ['Au'])

    def test_stacking_sequences_and_periodic_faults(self):
        dict_local_layers = {'HCP_basal_I1': 2, 'HCP_basal_I3': 6,
                             'FCC_111_ESF': 4, 'FCC_111_TB': 2}
        for model, count in dict_local_layers.items():
            phase, sequence_ref, sequence_fault = get_stacking_sequence(model)
            atoms_ref = create_stacking_model(sequence_ref, a=2.86 if phase == 'hcp' else 4.08,
                                               c=4.8, lat=phase)
            atoms_fault = create_stacking_model(sequence_fault, a=2.86 if phase == 'hcp' else 4.08,
                                                 c=4.8, lat=phase)
            self.assertEqual(get_stacking_letters(atoms_ref), sequence_ref)
            self.assertEqual(get_stacking_letters(atoms_fault), sequence_fault)
            self.assertEqual(get_periodic_fault_count(sequence_ref, phase), 0)
            self.assertEqual(get_periodic_fault_count(sequence_fault, phase), count)

    def test_symmetric_strain_selection(self):
        with tempfile.TemporaryDirectory() as tmp:
            path_root = Path(tmp)
            lstrain = np.linspace(-.003, .003, 101)
            self.write_dump(path_root / 'movie', [(np.diag([4. * (1. + strain), 5., 8.]),
                                                   np.zeros(3), np.array([[.1, .2, .3]]))
                                                  for strain in lstrain])
            (path_root / 'y_post_data.txt').write_text('header\n' + ''.join('%.12f 0\n' % (1. + strain)
                                                                          for strain in lstrain), encoding='utf-8')
            path_equil = path_root / 'y_stretch/fcc/dump/data.full_relax'
            path_equil.parent.mkdir(parents=True)
            write(str(path_equil), Atoms('Au', positions=[[.1, .2, .3]], cell=np.diag([4., 5., 8.]), pbc=True),
                  format='lammps-data', atom_style='atomic')
            dict_rule = dict(group='cij', phase='fcc', source='movie', format='lammps-dump-text', strain_points=5)
            lselected = read_rule_frames(path_root, dict_rule, ['Au'])
            self.assertEqual([row[1] for row in lselected], [0, 25, 50, 75, 100])
            np.testing.assert_allclose([row[4] for row in lselected], [-.003, -.0015, 0., .0015, .003])


if __name__ == '__main__':
    unittest.main()
