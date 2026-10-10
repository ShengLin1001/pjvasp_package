"""PeiN2p2Property epoch-selection tests; no scheduler or solver calls."""

import tempfile
import unittest
from pathlib import Path

from mymetal.ml.n2p2.workflow_property import PeiN2p2Property


class TestSelectExportEpochs(unittest.TestCase):

    @staticmethod
    def make_epoch_dirs(path_root, ldone_epoch, lincomplete_epoch):
        """Create epoch dirs with the 12 non-planar result markers (empty=incomplete)."""
        dir_scan = path_root / 'y_post' / '001' / 'properties' / 'y_epoch_scan' / 'y_dir'
        lmarkers = ([Path('y_stretch') / lat / 'p_post_stretch.txt' for lat in ['fcc', 'bcc', 'hcp']]
                    + [Path('y_Cij_energy') / lat / 'y_post_cij_energy.txt' for lat in ['fcc', 'bcc', 'hcp']]
                    + [Path('y_gsfe') / phase / t / 'y_post_gsfe.txt'
                       for phase, ltype in [('fcc', ['FCC_111', 'FCC_100']),
                                            ('hcp', ['HCP_basal', 'HCP_prism1w',
                                                     'HCP_pyr1w', 'HCP_pyr2'])]
                       for t in ltype])
        for epoch in ldone_epoch + lincomplete_epoch:
            for path_marker in lmarkers:
                path_marker = dir_scan / epoch / path_marker
                path_marker.parent.mkdir(parents=True, exist_ok=True)
                if epoch in ldone_epoch:
                    path_marker.write_text('result', encoding='utf-8')
        return dir_scan

    def test_done_epochs_selected_and_incomplete_skipped(self):
        with tempfile.TemporaryDirectory() as tmp:
            dir_scan = self.make_epoch_dirs(Path(tmp), ['000300', '000350'], ['000250'])
            ldone, lskip = PeiN2p2Property._select_export_epochs(dir_scan, if_planar=False)
            self.assertEqual(ldone, ['000300', '000350'])
            self.assertEqual(lskip, ['000250'])

    def test_named_epochs_filtered_and_missing_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            dir_scan = self.make_epoch_dirs(Path(tmp), ['000300', '000350'], [])
            ldone, _ = PeiN2p2Property._select_export_epochs(dir_scan, lepoch=[350], if_planar=False)
            self.assertEqual(ldone, ['000350'])
            with self.assertRaises(FileNotFoundError):
                PeiN2p2Property._select_export_epochs(dir_scan, lepoch=[500], if_planar=False)

    def test_missing_scan_root_rejected(self):
        with self.assertRaises(FileNotFoundError):
            PeiN2p2Property._select_export_epochs(Path('/unused/none'), if_planar=False)


if __name__ == '__main__':
    unittest.main()
