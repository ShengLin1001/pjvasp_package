"""Regression tests for the reviewed HOEC solve and small-shear controls."""

import io
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

from mymetal.build.workflow.hoec import generate_hoec_dirs
from mymetal.calculate.calmechanics.hoec import (
    MODES_HEX,
    get_hoec_modes,
    get_mode_strain_lists,
    get_model,
)
from mymetal.post.hoec_energy import select_solve_modes


class TestHoecReviewFixes(unittest.TestCase):
    """Keep the three Claude review regressions closed."""

    def test_default_hex_foec_solve_includes_new_pure_normal_modes(self):
        model = get_model("hex")
        lsolve = select_solve_modes(
            model,
            4,
            list(MODES_HEX),
            dict_modes=MODES_HEX,
        )

        self.assertEqual(len(lsolve), len(model.names(4)))
        self.assertTrue({"M21", "M22", "M23"}.issubset(lsolve))
        arr_system = model.system(4, [MODES_HEX[name] for name in lsolve])
        self.assertEqual(np.linalg.matrix_rank(arr_system, tol=1e-7), len(model.names(4)))
        lpure = [name for name, d_dir in MODES_HEX.items() if not any(d_dir[3:])]
        self.assertTrue(set(lpure).issubset(lsolve))
        arr_pure = model.system(4, [MODES_HEX[name] for name in lpure])
        arr_solve_pure = model.system(
            4,
            [MODES_HEX[name] for name in lsolve if name in lpure],
        )
        self.assertEqual(
            np.linalg.matrix_rank(arr_solve_pure, tol=1e-7),
            np.linalg.matrix_rank(arr_pure, tol=1e-7),
        )

    def test_non_positive_small_shear_controls_fail_before_preparation(self):
        path_root = Path(__file__).resolve().parents[1]
        for dict_bad in (
            {"shear_scale": 0.0},
            {"shear_scale": -0.5},
            {"shear_cap": 0.0},
            {"shear_cap": -0.03},
        ):
            with self.subTest(**dict_bad):
                with patch("mymetal.build.workflow.hoec.prepare_hoec_reference") as mock_prepare:
                    with redirect_stdout(io.StringIO()):
                        with self.assertRaises(SystemExit):
                            generate_hoec_dirs(
                                path_root=path_root,
                                srcdir="mymetal",
                                symmetry="hex",
                                small_shear=True,
                                **dict_bad,
                            )
                mock_prepare.assert_not_called()

    def test_zero_shear_cap_is_not_silently_skipped(self):
        with self.assertRaisesRegex(ValueError, "shear_cap"):
            get_mode_strain_lists(
                "hex",
                dict_modes={"M09": MODES_HEX["M09"]},
                shear_cap=0.0,
                lcap_modes=["M09"],
            )

    def test_explicit_empty_small_shear_targets_change_no_directions(self):
        dict_modes = get_hoec_modes(
            "hex",
            small_shear=True,
            lsmall_shear_modes=[],
        )

        self.assertEqual(dict_modes, MODES_HEX)

    def test_generator_passes_explicit_empty_small_shear_targets_unchanged(self):
        path_root = Path(__file__).resolve().parents[1]
        atoms_ref = MagicMock()
        atoms_ref.get_cell.return_value = np.eye(3)
        atoms_ref.__len__.return_value = 2

        with patch("mymetal.build.workflow.hoec.prepare_hoec_reference",
                   return_value=path_root / "CONTCAR"):
            with patch("mymetal.build.workflow.hoec.read", return_value=atoms_ref):
                with patch("mymetal.build.workflow.hoec.get_hoec_modes",
                           side_effect=[MODES_HEX, RuntimeError("mode resolution complete")]) \
                        as mock_modes:
                    with redirect_stdout(io.StringIO()):
                        with self.assertRaisesRegex(RuntimeError, "mode resolution complete"):
                            generate_hoec_dirs(
                                path_root=path_root,
                                srcdir="mymetal",
                                symmetry="hex",
                                small_shear=True,
                                small_shear_modes=[],
                            )

        self.assertEqual(mock_modes.call_args_list[1].args[3], [])


if __name__ == "__main__":
    unittest.main()
