#!/usr/bin/env python3
"""Build, compare, write, and render HCP Mg low-index surfaces.

The script is deterministic and VASP-free. It builds four HCP Mg cells with
the project helpers in :mod:`mymetal.build.bulk.create`:

* :func:`create_hcp_basal` -> Mg(0001) basal plane
* :func:`create_hcp_prism1` (``mode='wide'``) -> Mg(10-10) prism I, wide spacing
* :func:`create_hcp_prism1` (``mode='narrow'``) -> Mg(10-10) prism I, narrow spacing
* :func:`create_hcp_prism2` -> Mg(10-11) prism II

It then repeats each cell to a comparable thickness, prints a summary table,
writes one ``POSCAR`` per surface, runs an ASE round-trip check, and renders a
2x4 comparison figure reused by the documentation page.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ase import Atoms
from ase.visualize.plot import plot_atoms

from mymetal.build.bulk.create import (
    create_hcp_basal,
    create_hcp_prism1,
    create_hcp_prism2,
)
from mymetal.io.vasp import my_read_vasp, my_write_vasp
from mymetal.universal.plot.general import general_set_all_rcParams


# Mg HCP lattice constants (Å). a is the basal edge, c the axial height.
A_HCP_MG = 3.21
C_HCP_MG = 5.21
# Build enough repeats so every surface has a comparable slab thickness.
# These values are picked so the resulting c-axis length is in the 25-30 Å
# range before vacuum is added; the cells stay small enough for a CI smoke test.
SIZE_BASAL = (1, 1, 8)
SIZE_PRISM1_WIDE = (1, 1, 5)
SIZE_PRISM1_NARROW = (1, 1, 5)
SIZE_PRISM2 = (1, 1, 6)
VACUUM = 15.0


def build_hcp_cells() -> list[Atoms]:
    """Return the four deterministic HCP Mg cells in display order."""
    cells: list[Atoms] = [
        create_hcp_basal(a=A_HCP_MG, c=C_HCP_MG, size=SIZE_BASAL, symbol="Mg"),
        create_hcp_prism1(
            a=A_HCP_MG, c=C_HCP_MG, size=SIZE_PRISM1_WIDE,
            symbol="Mg", mode="wide",
        ),
        create_hcp_prism1(
            a=A_HCP_MG, c=C_HCP_MG, size=SIZE_PRISM1_NARROW,
            symbol="Mg", mode="narrow",
        ),
        create_hcp_prism2(a=A_HCP_MG, c=C_HCP_MG, size=SIZE_PRISM2, symbol="Mg"),
    ]
    for cell in cells:
        # Add vacuum along z so the side view reads as a slab, not a bulk rod.
        cell.center(vacuum=VACUUM, axis=2)
    return cells


def cell_label(idx: int) -> str:
    labels = [
        "Mg(0001) basal",
        "Mg(10-10) prism I (wide)",
        "Mg(10-10) prism I (narrow)",
        "Mg(10-11) prism II",
    ]
    return labels[idx]


def summarize_cell(idx: int, cell: Atoms) -> dict[str, object]:
    """Return a deterministic summary used by the table and the assertions."""
    lengths = np.linalg.norm(cell.cell.array, axis=1)
    return {
        "label": cell_label(idx),
        "formula": cell.get_chemical_formula(),
        "atoms": len(cell),
        "a_A": round(float(lengths[0]), 4),
        "b_A": round(float(lengths[1]), 4),
        "c_total_A": round(float(lengths[2]), 4),
        "pbc": cell.pbc.tolist(),
    }


def write_and_verify(cell: Atoms, path_poscar: Path) -> None:
    """Write a VASP file and verify a round trip through the project I/O."""
    my_write_vasp(path_poscar, cell, label="Mg slab")
    cell_read, scale = my_read_vasp(path_poscar)
    assert scale == 1.0
    assert cell_read.get_chemical_formula() == cell.get_chemical_formula()
    assert len(cell_read) == len(cell)
    assert np.array_equal(cell_read.pbc, cell.pbc)
    assert np.allclose(cell_read.cell.array, cell.cell.array, atol=1e-10)
    assert np.allclose(cell_read.positions, cell.positions, atol=1e-10)


def render_comparison(cells: list[Atoms], path_image: Path) -> None:
    """Render a 2x4 comparison figure (top row: side view, bottom row: top view)."""
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
        plot_atoms(
            cell.repeat((2, 1, 1)),
            axes[0, col],
            rotation="90x,0y,0z",
            show_unit_cell=0,
            radii=0.6,
        )
        axes[0, col].set_title(cell_label(col) + "\nside view")
        axes[0, col].set_axis_off()

        plot_atoms(
            cell.repeat((2, 2, 1)),
            axes[1, col],
            rotation="0x,0y,0z",
            show_unit_cell=1,
            radii=0.6,
        )
        axes[1, col].set_title("top view")
        axes[1, col].set_axis_off()
    fig.savefig(path_image, bbox_inches="tight")
    plt.close(fig)

    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError("comparison image was not created: " + str(path_image))
    image_rgb = plt.imread(path_image)[..., :3]
    if float(np.mean(np.abs(image_rgb - 1.0))) < 0.002:
        raise AssertionError("comparison image is effectively blank: " + str(path_image))


def run_example(path_output: Path) -> list[dict[str, object]]:
    """Build, write, render, and return the per-surface summary rows."""
    path_output = path_output.resolve()
    path_output.mkdir(parents=True, exist_ok=True)

    cells = build_hcp_cells()
    rows: list[dict[str, object]] = []
    for idx, cell in enumerate(cells):
        rows.append(summarize_cell(idx, cell))
        slug = cell_label(idx).split()[0].replace("(", "_").replace(")", "").replace("-", "")
        write_and_verify(cell, path_output / f"POSCAR_{slug}")

    path_image = path_output / "hcp_surfaces.png"
    render_comparison(cells, path_image)

    print("label                          formula  atoms  a_A      b_A      c_total_A  pbc")
    for row in rows:
        print(
            f"{row['label']:<30}  {row['formula']:<7}  "
            f"{row['atoms']:<5}  "
            f"{row['a_A']:<7.4f}  "
            f"{row['b_A']:<7.4f}  "
            f"{row['c_total_A']:<9.4f}  "
            f"{row['pbc']}"
        )
    print("wrote: " + str(path_image))
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build and compare HCP Mg low-index surfaces.")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("hcp-surfaces-output"),
        help="Directory for POSCAR files and the comparison PNG.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    rows = run_example(parse_args().output)
    formulas = [row["formula"] for row in rows]
    assert all(f.startswith("Mg") for f in formulas), formulas
    assert all(row["atoms"] > 0 for row in rows)
    assert all(row["pbc"] == [True, True, True] for row in rows)
