#!/usr/bin/env python3
"""Build and compare FCC / BCC / HCP / diamond conventional cells.

This tutorial is VASP-free. It demonstrates how ``ase.build.bulk`` plus the
project helpers in :mod:`mymetal.build.bulk.create` produce the four most
common crystal structures used as bulk references:

* FCC Cu (cubic, ``a = 3.61 A``)
* BCC Fe (cubic, ``a = 2.87 A``)
* HCP Mg (primitive, ``a = 3.21 A``, ``c = 5.21 A``)
* Diamond Si (primitive, ``a = 5.43 A``)

Each cell is built, repeated to a comparable number of atoms, written to a
``POSCAR``, read back, and rendered in a 2x4 comparison figure reused by the
documentation page.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ase import Atoms
from ase.build import bulk
from ase.visualize.plot import plot_atoms

from mymetal.io.vasp import my_read_vasp, my_write_vasp
from mymetal.universal.plot.general import general_set_all_rcParams


# Lattice constants (Å). These are textbook ambient values, kept here so the
# script is fully deterministic and free of network access.
A_FCC_CU = 3.61
A_BCC_FE = 2.87
A_HCP_MG = 3.21
C_HCP_MG = 5.21
A_DIA_SI = 5.43


def build_bulk_cells() -> list[Atoms]:
    """Return four conventional bulk cells in display order."""
    cells: list[Atoms] = [
        bulk("Cu", "fcc", a=A_FCC_CU, cubic=True),
        bulk("Fe", "bcc", a=A_BCC_FE, cubic=True),
        bulk("Mg", "hcp", a=A_HCP_MG, covera=C_HCP_MG / A_HCP_MG),
        bulk("Si", "diamond", a=A_DIA_SI, cubic=True),
    ]
    return cells


def cell_label(idx: int) -> str:
    labels = [
        "FCC Cu (cubic)",
        "BCC Fe (cubic)",
        "HCP Mg (primitive)",
        "Diamond Si (cubic)",
    ]
    return labels[idx]


def repeat_for_display(cell: Atoms) -> Atoms:
    """Repeat each conventional cell to a comparable visual size."""
    if cell.get_chemical_formula().startswith("Cu"):
        return cell.repeat((2, 2, 2))
    if cell.get_chemical_formula().startswith("Fe"):
        return cell.repeat((2, 2, 2))
    if cell.get_chemical_formula().startswith("Mg"):
        return cell.repeat((3, 3, 2))
    if cell.get_chemical_formula().startswith("Si"):
        return cell.repeat((1, 1, 1))
    return cell


def summarize_cell(idx: int, cell: Atoms) -> dict[str, object]:
    """Return a deterministic summary used by the table and the assertions."""
    lengths = np.linalg.norm(cell.cell.array, axis=1)
    angles = cell.cell.angles()
    return {
        "label": cell_label(idx),
        "formula": cell.get_chemical_formula(),
        "atoms": len(cell),
        "a_A": round(float(lengths[0]), 4),
        "b_A": round(float(lengths[1]), 4),
        "c_A": round(float(lengths[2]), 4),
        "alpha_deg": round(float(angles[0]), 4),
        "beta_deg": round(float(angles[1]), 4),
        "gamma_deg": round(float(angles[2]), 4),
        "volume_A3": round(float(cell.get_volume()), 4),
        "pbc": cell.pbc.tolist(),
    }


def write_and_verify(cell: Atoms, path_poscar: Path) -> None:
    """Write a VASP file and verify a round trip through the project I/O."""
    my_write_vasp(path_poscar, cell, label="bulk structure")
    cell_read, scale = my_read_vasp(path_poscar)
    assert scale == 1.0
    assert cell_read.get_chemical_formula() == cell.get_chemical_formula()
    assert len(cell_read) == len(cell)
    assert np.array_equal(cell_read.pbc, cell.pbc)
    assert np.allclose(cell_read.cell.array, cell.cell.array, atol=1e-10)
    assert np.allclose(cell_read.positions, cell.positions, atol=1e-10)


def render_comparison(cells: list[Atoms], path_image: Path) -> None:
    """Render a 2x4 comparison figure (top: 3D oblique, bottom: side view)."""
    path_image.parent.mkdir(parents=True, exist_ok=True)

    general_set_all_rcParams(
        backend="Agg",
        axes_grid=False,
        figure_subp=(2, 4),
        figure_one_figsize=(3.6, 4.4),
        figure_left=0.50,
        figure_top=0.50,
        figure_wspace=0.55,
        figure_hspace=0.45,
        font_family=["DejaVu Sans"],
        fontsize=10,
        axes_titlepad=6,
    )
    fig, axes = plt.subplots(2, 4)
    for col, cell in enumerate(cells):
        display = repeat_for_display(cell)
        plot_atoms(
            display,
            axes[0, col],
            rotation="10x,20y,0z",
            show_unit_cell=1,
            radii=0.6,
        )
        axes[0, col].set_title(cell_label(col) + "\n3D view")
        axes[0, col].set_axis_off()

        plot_atoms(
            display,
            axes[1, col],
            rotation="90x,0y,0z",
            show_unit_cell=1,
            radii=0.6,
        )
        axes[1, col].set_title("side view")
        axes[1, col].set_axis_off()
    fig.savefig(path_image, bbox_inches="tight")
    plt.close(fig)

    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError("comparison image was not created: " + str(path_image))
    image_rgb = plt.imread(path_image)[..., :3]
    if float(np.mean(np.abs(image_rgb - 1.0))) < 0.002:
        raise AssertionError("comparison image is effectively blank: " + str(path_image))


def run_example(path_output: Path) -> list[dict[str, object]]:
    """Build, write, render, and return the per-structure summary rows."""
    path_output = path_output.resolve()
    path_output.mkdir(parents=True, exist_ok=True)

    cells = build_bulk_cells()
    rows: list[dict[str, object]] = []
    for idx, cell in enumerate(cells):
        rows.append(summarize_cell(idx, cell))
        slug = cell_label(idx).split()[0]
        write_and_verify(cell, path_output / f"POSCAR_{slug}")

    path_image = path_output / "bulk_structures.png"
    render_comparison(cells, path_image)

    print("label                      formula  atoms  a_A     b_A     c_A     "
          "alpha   beta    gamma   volume_A3  pbc")
    for row in rows:
        print(
            f"{row['label']:<26}  {row['formula']:<7}  "
            f"{row['atoms']:<5}  "
            f"{row['a_A']:<6.4f}  "
            f"{row['b_A']:<6.4f}  "
            f"{row['c_A']:<6.4f}  "
            f"{row['alpha_deg']:<6.2f}  "
            f"{row['beta_deg']:<6.2f}  "
            f"{row['gamma_deg']:<6.2f}  "
            f"{row['volume_A3']:<9.4f}  "
            f"{row['pbc']}"
        )
    print("wrote: " + str(path_image))
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build and compare FCC / BCC / HCP / diamond bulk cells.")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("bulk-structures-output"),
        help="Directory for POSCAR files and the comparison PNG.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    rows = run_example(parse_args().output)
    formulas = [row["formula"] for row in rows]
    assert all(row["atoms"] > 0 for row in rows)
    assert all(row["pbc"] == [True, True, True] for row in rows)
    # Cubic cells must have 90 degree angles.
    assert rows[0]["alpha_deg"] == 90.0 and rows[0]["gamma_deg"] == 90.0
    assert rows[1]["alpha_deg"] == 90.0 and rows[1]["gamma_deg"] == 90.0
    assert rows[3]["alpha_deg"] == 90.0 and rows[3]["gamma_deg"] == 90.0
    # HCP primitive cell must have 120 degree gamma.
    assert abs(rows[2]["gamma_deg"] - 120.0) < 1e-3
