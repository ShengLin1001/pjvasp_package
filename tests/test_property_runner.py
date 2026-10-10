"""Safety guards in the property runner; never call LAMMPS or Slurm."""

from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

PATH_RUNNER = Path(__file__).resolve().parents[1] / 'lmp_utils/lmp_universal/pei_lmp_run_properties'


class TestPropertyRunner(unittest.TestCase):
    def run_guard(self, path_root):
        return subprocess.run(['bash', str(PATH_RUNNER), 'eam', 'unused', sys.executable,
                               'mass 1 196.97', str(path_root / 'template')],
                              cwd=path_root, capture_output=True, text=True)

    def test_legacy_output_not_deleted(self):
        with tempfile.TemporaryDirectory() as tmp:
            path_root = Path(tmp)
            path_file = path_root / 'y_stretch/fcc/keep.txt'
            path_file.parent.mkdir(parents=True)
            path_file.write_text('important', encoding='utf-8')
            result = self.run_guard(path_root)
            self.assertEqual(result.returncode, 1)
            self.assertIn('未登记', result.stdout)
            self.assertEqual(path_file.read_text(encoding='utf-8'), 'important')

    def test_changed_identity_not_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            path_root = Path(tmp)
            (path_root / '.pei_properties_signature').write_text('different\n', encoding='utf-8')
            result = self.run_guard(path_root)
            self.assertEqual(result.returncode, 1)
            self.assertIn('已改变', result.stdout)
            self.assertEqual((path_root / '.pei_properties_signature').read_text(encoding='utf-8'), 'different\n')

    def test_no_recursive_cleanup(self):
        text = PATH_RUNNER.read_text(encoding='utf-8')
        self.assertNotIn('rm -rf', text)
        self.assertNotIn('rmtree', text)
        self.assertIn('flock -n', text)
        self.assertIn('pei_slurm_univ_launch_retry', text)


if __name__ == '__main__':
    unittest.main()
