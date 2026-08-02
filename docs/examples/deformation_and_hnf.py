#!/usr/bin/env python3
"""Deformation gradient and Hermite Normal Form, side by side.

This tutorial is VASP-free and fully deterministic. It demonstrates two
linear-algebra helpers used across the project:

* :func:`mymetal.calculate.calmechanics.deformation.cal_deform_matrix` -- given
  an initial cell matrix and a deformed cell matrix (both as
  row-vector ``np.ndarray`` of shape ``(3, 3)``), returns the deformation
  gradient ``F = deformed_cell @ inv(initial_cell)``.

* :func:`mymetal.calculate.calmath.matrix.hermite_normal_form` -- given an
  integer ``np.ndarray``, returns its Hermite Normal Form (HNF) using pure
  integer arithmetic.

The figure has three panels in a single row:

* Panel 1 -- FCC Cu conventional cell drawn before (dashed) and after (solid)
  a uniaxial 5% stretch along x, with the deformation gradient annotated.
* Panel 2 -- a 3x3 integer supercell matrix ``M_c`` rendered as a coloured
  grid with the integer entries annotated.
* Panel 3 -- the Hermite Normal Form ``HNF(M_c)`` rendered the same way.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ase.build import bulk

from mymetal.calculate.calmath.matrix import hermite_normal_form
from mymetal.calculate.calmechanics.deformation import cal_deform_matrix
from mymetal.universal.plot.general import general_set_all_rcParams


# Lattice constant (Å). Textbook ambient value, hard-coded so the script is
# deterministic and free of network access.
A_FCC_CU = 3.61

# Three example integer matrices representing supercell transformations.
MATRICES = [
    np.array([[2, 1, 0], [0, 2, 1], [0, 0, 3]], dtype=np.int64),
    np.array([[1, 2, 3], [0, 1, 4], [0, 0, 1]], dtype=np.int64),
    np.array([[3, 1, 2], [1, 2, 1], [2, 1, 1]], dtype=np.int64),
]
MATRIX_NAMES = ["M_a (upper-triangular)", "M_b (unit upper-triangular)",
                "M_c (general dense)"]


# --------------------------------------------------------------------------- #
# Build step
# --------------------------------------------------------------------------- #
def build_deformation_pair():
    """Return (initial_Atoms, deformed_Atoms, initial_cell, deformed_cell, F).

    The deformed cell is produced by stretching the FCC Cu conventional cell
    by 5% along x, which makes the expected deformation gradient
    ``F = diag(1.05, 1.0, 1.0)``.
    """
    atoms_init = bulk("Cu", "fcc", a=A_FCC_CU, cubic=True)
    atoms_defm = atoms_init.copy()
    # Uniaxial 5% stretch along x: scale the first lattice vector.
    atoms_defm.cell[0] = atoms_defm.cell[0] * 1.05

    initial_cell = atoms_init.cell.array.copy()
    deformed_cell = atoms_defm.cell.array.copy()
    F = cal_deform_matrix(initial_cell, deformed_cell)
    return atoms_init, atoms_defm, initial_cell, deformed_cell, F


def build_hnf_table():
    """Return list of dicts with original matrix, HNF, and a few diagnostics."""
    rows = []
    for name, mat in zip(MATRIX_NAMES, MATRICES):
        hnf = hermite_normal_form(mat)
        rows.append({
            "name": name,
            "matrix": mat,
            "hnf": hnf,
            "det": int(round(np.linalg.det(mat))),
        })
    return rows


# --------------------------------------------------------------------------- #
# Render step
# --------------------------------------------------------------------------- #
def _draw_cell_box(ax, cell, color, linestyle, label):
    """Draw the projection of a 3x3 cell onto the xy-plane as a parallelogram."""
    a = cell[0, :2]
    b = cell[1, :2]
    quad = np.array([np.zeros(2), a, a + b, b, np.zeros(2)])
    ax.plot(quad[:, 0], quad[:, 1], color=color, linestyle=linestyle,
            linewidth=2.0, label=label)
    # Corner markers for clarity.
    ax.scatter(quad[:-1, 0], quad[:-1, 1], color=color, s=30, zorder=5)


def render_figure(path_image: Path, defm_data, hnf_rows) -> None:
    """Render the three-panel figure described in the module docstring."""
    path_image.parent.mkdir(parents=True, exist_ok=True)

    general_set_all_rcParams(
        backend="Agg",
        axes_grid=False,
        figure_subp=(1, 3),
        figure_one_figsize=(4.2, 4.0),
        figure_left=0.50,
        figure_top=0.85,
        figure_wspace=0.55,
        figure_hspace=0.40,
        font_family=["DejaVu Sans"],
        fontsize=10,
        axes_titlepad=8,
    )
    fig, axes = plt.subplots(1, 3)

    # ---- Panel 1: deformation schematic (xy-projection) ------------------ #
    ax = axes[0]
    _, _, initial_cell, deformed_cell, F = defm_data
    _draw_cell_box(ax, initial_cell, color="#1f77b4", linestyle="--",
                   label="initial")
    _draw_cell_box(ax, deformed_cell, color="#d62728", linestyle="-",
                   label="deformed (x +5%)")
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_title("FCC Cu: uniaxial stretch along x")
    ax.set_xlabel("x (Å)")
    ax.set_ylabel("y (Å)")
    ax.legend(loc="upper right", fontsize=9, framealpha=0.9)
    # Annotate the deformation gradient.
    f00 = float(F[0, 0])
    f11 = float(F[1, 1])
    f22 = float(F[2, 2])
    ax.text(
        0.02, 0.02,
        f"F = diag({f00:.3f}, {f11:.3f}, {f22:.3f})",
        transform=ax.transAxes, fontsize=9, va="bottom", ha="left",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.6", alpha=0.9),
    )

    # ---- Panels 2 and 3: matrix (c) original and HNF --------------------- #
    # Use the third (general dense) matrix for the visual comparison.
    matrix_c = hnf_rows[2]["matrix"]
    hnf_c = hnf_rows[2]["hnf"]
    _render_matrix_heatmap(axes[1], matrix_c, "M_c (original)")
    _render_matrix_heatmap(axes[2], hnf_c, "HNF(M_c)")

    fig.savefig(path_image, bbox_inches="tight", dpi=150)
    plt.close(fig)

    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError("image was not created: " + str(path_image))
    image_rgb = plt.imread(path_image)[..., :3]
    if float(np.mean(np.abs(image_rgb - 1.0))) < 0.002:
        raise AssertionError("image is effectively blank: " + str(path_image))


def _render_matrix_heatmap(ax, matrix, title):
    """Render a small integer matrix as a coloured grid with annotations."""
    m = np.asarray(matrix)
    vmax = max(abs(int(m.min())), abs(int(m.max())), 1)
    im = ax.imshow(m, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(m.shape[1]))
    ax.set_yticks(range(m.shape[0]))
    ax.set_title(title)
    # Annotate each cell with its integer value.
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            val = int(m[i, j])
            text_color = "white" if abs(val) > vmax * 0.55 else "black"
            ax.text(j, i, str(val), ha="center", va="center",
                    color=text_color, fontsize=11, fontweight="bold")
    ax.set_xlabel("column")
    ax.set_ylabel("row")


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def run_example(path_output: Path) -> dict:
    """Build, render, print, and assert; return a summary dict."""
    path_output = path_output.resolve()
    path_output.mkdir(parents=True, exist_ok=True)

    defm_data = build_deformation_pair()
    hnf_rows = build_hnf_table()

    # ---- Print deformation summary --------------------------------------- #
    _, _, initial_cell, deformed_cell, F = defm_data
    print("== cal_deform_matrix ==")
    print("initial_cell (Å):")
    print(np.array2string(initial_cell, precision=4, suppress_small=True))
    print("deformed_cell (Å):")
    print(np.array2string(deformed_cell, precision=4, suppress_small=True))
    print("deformation gradient F = deformed_cell @ inv(initial_cell):")
    print(np.array2string(np.asarray(F), precision=4, suppress_small=True))
    print()

    # ---- Print HNF summary ----------------------------------------------- #
    print("== hermite_normal_form ==")
    print(f"{'matrix':<28} {'det':>6}   original -> HNF")
    for row in hnf_rows:
        mat_str = np.array2string(row["matrix"]).replace("\n", " ")
        hnf_str = np.array2string(row["hnf"]).replace("\n", " ")
        print(f"{row['name']:<28} {row['det']:>6}   {mat_str}  ->  {hnf_str}")
    print()

    # ---- Render ---------------------------------------------------------- #
    path_image = path_output / "deformation_and_hnf.png"
    render_figure(path_image, defm_data, hnf_rows)
    print("wrote: " + str(path_image))

    # ---- Deterministic assertions ---------------------------------------- #
    # Uniaxial 5% stretch along x of a cubic cell -> F = diag(1.05, 1, 1).
    assert abs(float(F[0, 0]) - 1.05) < 1e-9, f"F[0,0] != 1.05, got {F[0, 0]}"
    assert abs(float(F[1, 1]) - 1.00) < 1e-9
    assert abs(float(F[2, 2]) - 1.00) < 1e-9
    # Off-diagonal entries must remain zero for a pure axis-aligned stretch.
    off_diag_mask = ~np.eye(3, dtype=bool)
    assert np.allclose(np.asarray(F)[off_diag_mask], 0.0, atol=1e-9)

    # HNF of the unit upper-triangular matrix (M_b) must be the identity.
    assert np.array_equal(hnf_rows[1]["hnf"],
                          np.eye(3, dtype=np.int64)), \
        "HNF of unit upper-triangular matrix must be the identity"

    # HNF of the upper-triangular M_a with positive diagonal must equal M_a.
    assert np.array_equal(hnf_rows[0]["hnf"], MATRICES[0]), \
        "HNF of upper-triangular matrix with positive diagonal must equal itself"

    # HNF results must be integer-valued (the function uses int64 arithmetic).
    for row in hnf_rows:
        assert row["hnf"].dtype == np.int64, \
            f"HNF of {row['name']} must be int64"
        assert np.all(np.asarray(row["hnf"]) == np.round(np.asarray(row["hnf"]))), \
            f"HNF of {row['name']} must contain only integers"

    return {
        "image": str(path_image),
        "F00": float(F[0, 0]),
        "hnf_rows": hnf_rows,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Deformation gradient and Hermite Normal Form example.")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("deformation-and-hnf-output"),
        help="Directory for the comparison PNG.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    summary = run_example(parse_args().output)
    print("\nOK: F[0,0] = %.6f, image = %s"
          % (summary["F00"], summary["image"]))
