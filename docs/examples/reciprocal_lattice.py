#!/usr/bin/env python3
"""Reciprocal lattice vectors and RK-based k-point mesh sizing for bulk cells.

This tutorial is VASP-free. It complements ``kpoints_sampling.py`` (which
covers Monkhorst-Pack vs Gamma-centered positions and an RK scan on a slab)
by focusing on two things the other tutorial does not show:

1. The two equivalent implementations of the reciprocal lattice in
   :mod:`mymetal.calculate.calqm.kpoints` --
   :func:`cal_reciprocal_matrix` (cross product) and
   :func:`cal_reciprocal_matrix2` (matrix inversion) -- compared side by
   side on three bulk cells so the agreement is visible and asserted.
2. How :func:`get_size_by_distance` maps a single ``RK`` product onto very
   different k-point meshes for FCC, BCC, and HCP conventional cells,
   driven by the spread of their reciprocal lattice vector lengths.

The script renders a two-panel figure reused by the documentation page:

* Left:  grouped bar chart of ``|b1|``, ``|b2|``, ``|b3|`` for FCC Cu,
  BCC Fe, and HCP Mg.
* Right: grouped bar chart of the old-VASP k-point mesh ``(nx, ny, nz)``
  returned by ``get_size_by_distance`` at ``rk = 80`` for the same three
  cells.
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

from mymetal.calculate.calqm.kpoints import (
    cal_reciprocal_matrix,
    cal_reciprocal_matrix2,
    get_size_by_distance,
)
from mymetal.universal.plot.general import general_set_all_rcParams


# Lattice constants (A). Textbook ambient values, kept inline so the script
# is fully deterministic and needs no network access.
A_FCC_CU = 3.61
A_BCC_FE = 2.87
A_HCP_MG = 3.21
C_HCP_MG = 5.21

# Single RK product used for the mesh-sizing comparison. 80 is a typical
# "medium-density" bulk value and makes the differences between the three
# lattices easy to read off the bar chart.
RK_TARGET = 80


def build_bulk_cells() -> list[Atoms]:
    """Return FCC Cu, BCC Fe, HCP Mg conventional cells in display order."""
    return [
        bulk("Cu", "fcc", a=A_FCC_CU, cubic=True),
        bulk("Fe", "bcc", a=A_BCC_FE, cubic=True),
        bulk("Mg", "hcp", a=A_HCP_MG, covera=C_HCP_MG / A_HCP_MG),
    ]


def cell_labels() -> list[str]:
    """Display labels in the same order as :func:`build_bulk_cells`."""
    return ["FCC Cu", "BCC Fe", "HCP Mg"]


def reciprocal_pair(cell_matrix: np.ndarray) -> dict[str, np.ndarray]:
    """Return reciprocal vectors from cross product and matrix inversion.

    Both methods use ``scale = 2*pi`` so the units are ``2*pi/A``, matching
    what VASP reports. The cross-product method returns three separate
    vectors; we stack them into a ``(3, 3)`` array whose rows are
    ``(b1, b2, b3)`` so it can be compared element-wise with the
    matrix-inversion result.
    """
    b1, b2, b3 = cal_reciprocal_matrix(cell_matrix, scale=2.0 * np.pi)
    cross = np.vstack([b1, b2, b3])
    inv = cal_reciprocal_matrix2(cell_matrix, scale=2.0 * np.pi)
    return {"cross": cross, "inv": inv}


def summarize_cell(idx: int, atoms: Atoms) -> dict[str, object]:
    """Return a deterministic summary row for one bulk cell."""
    cell_matrix = np.asarray(atoms.get_cell())
    recip = reciprocal_pair(cell_matrix)
    b_lengths = np.linalg.norm(recip["cross"], axis=1)
    old_k, new_k = get_size_by_distance(atoms=atoms, rk=RK_TARGET)
    old_k = np.asarray(old_k)
    new_k = np.asarray(new_k)
    return {
        "label": cell_labels()[idx],
        "formula": atoms.get_chemical_formula(),
        "atoms": len(atoms),
        "cell_matrix": cell_matrix,
        "recip_cross": recip["cross"],
        "recip_inv": recip["inv"],
        "b_lengths": b_lengths,
        "old_mesh": old_k,
        "new_mesh": new_k,
    }


def render_figure(rows: list[dict[str, object]], path_image: Path) -> None:
    """Render the two-panel comparison figure.

    Panel A (left):  grouped bars of ``|b1|``, ``|b2|``, ``|b3|`` per cell.
    Panel B (right): grouped bars of ``(nx, ny, nz)`` from
    ``get_size_by_distance`` at ``rk = RK_TARGET`` per cell.
    """
    path_image.parent.mkdir(parents=True, exist_ok=True)

    general_set_all_rcParams(
        backend="Agg",
        axes_grid=True,
        figure_subp=(1, 2),
        figure_one_figsize=(5.6, 4.4),
        figure_left=0.60,
        figure_top=0.55,
        figure_wspace=1.0,
        figure_hspace=0.6,
        font_family=["DejaVu Sans"],
        fontsize=11,
        axes_titlepad=6,
    )
    fig, axes = plt.subplots(1, 2)

    labels = [row["label"] for row in rows]
    x = np.arange(len(labels))
    bar_w = 0.25
    colors = ["#1f77b4", "#d62728", "#2ca02c"]

    # Panel A: reciprocal lattice vector magnitudes.
    ax = axes[0]
    for i in range(3):
        vals = np.array([row["b_lengths"][i] for row in rows])
        ax.bar(
            x + (i - 1) * bar_w,
            vals,
            bar_w,
            color=colors[i],
            label=rf"$|b_{{{i+1}}}|$",
        )
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel(r"reciprocal vector length ($2\pi/\mathrm{\AA}$)")
    ax.set_title("Reciprocal lattice vectors")
    ax.legend(fontsize=9)

    # Panel B: k-point mesh sizes from get_size_by_distance.
    ax = axes[1]
    for i in range(3):
        vals = np.array([row["old_mesh"][i] for row in rows])
        ax.bar(
            x + (i - 1) * bar_w,
            vals,
            bar_w,
            color=colors[i],
            label=rf"$n_{{{['x','y','z'][i]}}}$",
        )
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel(rf"k-point mesh at $R_k = {RK_TARGET}$")
    ax.set_title("RK-based mesh sizing (old VASP)")
    ax.legend(fontsize=9)

    fig.savefig(path_image, bbox_inches="tight")
    plt.close(fig)

    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError("reciprocal image was not created: " + str(path_image))
    image_rgb = plt.imread(path_image)[..., :3]
    if float(np.mean(np.abs(image_rgb - 1.0))) < 0.002:
        raise AssertionError("reciprocal image is effectively blank: " + str(path_image))


def run_example(path_output: Path) -> list[dict[str, object]]:
    """Build cells, compute reciprocal vectors and meshes, render, return rows."""
    path_output = path_output.resolve()
    path_output.mkdir(parents=True, exist_ok=True)

    cells = build_bulk_cells()
    rows: list[dict[str, object]] = []
    for idx, atoms in enumerate(cells):
        rows.append(summarize_cell(idx, atoms))

    path_image = path_output / "reciprocal_lattice.png"
    render_figure(rows, path_image)

    # Summary table.
    print("label    formula  |b1|     |b2|     |b3|     old_mesh        new_mesh")
    for row in rows:
        bl = row["b_lengths"]
        print(
            f"{row['label']:<8} {row['formula']:<7}  "
            f"{bl[0]:<6.4f}  {bl[1]:<6.4f}  {bl[2]:<6.4f}  "
            f"{row['old_mesh'].tolist()}  {row['new_mesh'].tolist()}"
        )
    print("wrote: " + str(path_image))
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Reciprocal lattice vectors and RK-based k-point mesh sizing."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("reciprocal-output"),
        help="Directory for the reciprocal-lattice PNG.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    result = run_example(parse_args().output)
    # Both reciprocal implementations must agree to 1e-10 for every cell.
    for row in result:
        assert np.allclose(row["recip_cross"], row["recip_inv"], atol=1e-10), (
            f"reciprocal mismatch for {row['label']}"
        )
        # K-point meshes must be positive integers in every direction.
        assert np.all(row["old_mesh"] >= 1), f"old mesh non-positive for {row['label']}"
        assert np.all(row["new_mesh"] >= 1), f"new mesh non-positive for {row['label']}"
        assert np.all(np.mod(row["old_mesh"], 1) == 0), (
            f"old mesh non-integer for {row['label']}"
        )
        # Reciprocal vector lengths must be positive and finite.
        assert np.all(np.isfinite(row["b_lengths"])) and np.all(row["b_lengths"] > 0), (
            f"bad reciprocal lengths for {row['label']}"
        )
