import unittest
from pathlib import Path

import numpy as np

from mymetal.ml.n2p2.workflow import PeiN2p2


class TestN2p2WorkflowUniformSelection(unittest.TestCase):

    def test_full_numeric_directory_name_controls_sort(self):
        ldir = [
            Path("+0.0200"),
            Path("-0.0100"),
            Path("+0.0000"),
            Path("-0.0200"),
            Path("+0.0100"),
        ]

        lselected = PeiN2p2._select_n_uniform_dirs(
            ldir,
            select_n_uniform_point=None,
        )

        self.assertEqual(
            [path_dir.name for path_dir in lselected],
            [
                "-0.0200",
                "-0.0100",
                "+0.0000",
                "+0.0100",
                "+0.0200",
            ],
        )

    def test_uniform_selection_includes_endpoints_and_midpoint(self):
        ldir = [Path(f"{strain:+.4f}")
                for strain in np.linspace(-0.12, 0.12, 97)]

        lselected = PeiN2p2._select_n_uniform_dirs(
            ldir,
            select_n_uniform_point=21,
        )

        self.assertEqual(len(lselected), 21)
        self.assertEqual(lselected[0].name, "-0.1200")
        self.assertEqual(lselected[10].name, "+0.0000")
        self.assertEqual(lselected[-1].name, "+0.1200")

    def test_data_tag_list_defaults_to_all_points(self):
        dict_config = PeiN2p2._get_data_tag_config({
            'A11-1': ['y_stretch'],
        })

        self.assertEqual(dict_config['A11-1']['lsubdir'], ['y_stretch'])
        self.assertIsNone(dict_config['A11-1']['select_n_uniform_point'])

    def test_data_tag_config_provides_uniform_point_count(self):
        dict_config = PeiN2p2._get_data_tag_config({
            'A11-2': {
                'lsubdir': ['y_hoec_energy/y_hoec_energy_M01'],
                'select_n_uniform_point': 21,
            },
        })

        self.assertEqual(dict_config['A11-2']['select_n_uniform_point'], 21)

    def test_selection_keeps_short_numeric_directory_list(self):
        ldir = [Path(name) for name in ["10", "2", "1"]]

        lselected = PeiN2p2._select_n_uniform_dirs(
            ldir,
            select_n_uniform_point=21,
        )

        self.assertEqual([path_dir.name for path_dir in lselected], ["1", "2", "10"])

    def test_selection_keeps_lexical_fallback_for_non_numeric_names(self):
        ldir = [Path(name) for name in ["y_full_relax", "reference"]]

        lselected = PeiN2p2._select_n_uniform_dirs(
            ldir,
            select_n_uniform_point=None,
        )

        self.assertEqual(
            [path_dir.name for path_dir in lselected],
            ["reference", "y_full_relax"],
        )

    def test_selection_rejects_non_positive_count(self):
        with self.assertRaises(ValueError):
            PeiN2p2._select_n_uniform_dirs([], select_n_uniform_point=0)


if __name__ == "__main__":
    unittest.main()
