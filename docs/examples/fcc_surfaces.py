#!/usr/bin/env python3
"""Build, compare, write, and render FCC Cu(100)/(110)/(111) slabs.

The script is intentionally deterministic and VASP-free: it builds three
12-layer Cu slabs with :func:`mymetal.build.film.stretch.generate_film`,
prints a comparable summary table, writes one ``POSCAR`` per surface, runs
an ASE round-trip check, and renders a 3-column comparison figure that is
reused by the documentation page.
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

from mymetal.build.film.stretch import generate_film
from mymetal.io.vasp import my_read_vasp, my_write_vasp
from mymetal.universal.plot.general import general_set_all_rcParams


# Cu conventional lattice constant (Å). Using the same a_fcc for every
# surface makes the layer-by-layer comparison meaningful.
A_FCC_CU = 3.61
NUM_LAYERS = 12
VACUUM = 15.0

# Ordered (label, miller) list used by every helper below.
FCC_SURFACES: list[tuple[str, tuple[int, int, int]]] = [
    ("Cu(100)", (1, 0, 0)),
    ("Cu(110)", (1, 1, 0)),
    ("Cu(111)", (1, 1, 1)),
]


def build_fcc_slabs() -> list[Atoms]:
    """Return the three deterministic slabs in :data:`FCC_SURFACES` order."""
    slabs: list[Atoms] = []
    for _label, miller in FCC_SURFACES:
        slab = generate_film(
            symbols="Cu",
            structure="fcc",
            num_layers=NUM_LAYERS,
            my_vacuum=VACUUM,
            slice_plane=miller,
            a_fcc=A_FCC_CU,
        )
        slabs.append(slab)
    return slabs


def summarize_slab(label: str, slab: Atoms) -> dict[str, object]:
    """Return a deterministic summary used by the table and the assertions."""
    cell = slab.cell.array
    lengths = np.linalg.norm(cell, axis=1)
    return {
        "label": label,
        "formula": slab.get_chemical_formula(),
        "atoms": len(slab),
        "a_in_plane_A": round(float(lengths[0]), 4),
        "b_in_plane_A": round(float(lengths[1]), 4),
        "c_total_A": round(float(lengths[2]), 4),
        "pbc": slab.pbc.tolist(),
    }


def write_and_verify(slab: Atoms, path_poscar: Path) -> None:
    """Write a VASP file and verify a round trip through the project I/O."""
    my_write_vasp(path_poscar, slab, label="Cu slab")
    slab_read, scale = my_read_vasp(path_poscar)
    assert scale == 1.0
    assert slab_read.get_chemical_formula() == slab.get_chemical_formula()
    assert len(slab_read) == len(slab)
    assert np.array_equal(slab_read.pbc, slab.pbc)
    assert np.allclose(slab_read.cell.array, slab.cell.array, atol=1e-10)
    assert np.allclose(slab_read.positions, slab.positions, atol=1e-10)


def render_comparison(slabs: list[Atoms], path_image: Path) -> None:
    """Render a 3-column (side view / top view / cell box) comparison figure."""
    path_image.parent.mkdir(parents=True, exist_ok=True)

    general_set_all_rcParams(
        backend="Agg",
        axes_grid=False,
        figure_subp=(2, 3),
        figure_one_figsize=(4.4, 4.6),
        figure_left=0.55,
        figure_top=0.55,
        figure_wspace=0.65,
        figure_hspace=0.50,
        font_family=["DejaVu Sans"],
        fontsize=11,
        axes_titlepad=6,
    )
    fig, axes = plt.subplots(2, 3)
    for col, (label, _miller) in enumerate(FCC_SURFACES):
        slab = slabs[col]
        plot_atoms(
            slab.repeat((3, 1, 1)),
            axes[0, col],
            rotation="90x,0y,0z",
            show_unit_cell=0,
            radii=0.7,
        )
        axes[0, col].set_title(label + "\nside view (3x in-plane repeat)")
        axes[0, col].set_axis_off()

        plot_atoms(
            slab.repeat((2, 2, 1)),
            axes[1, col],
            rotation="0x,0y,0z",
            show_unit_cell=1,
            radii=0.7,
        )
        axes[1, col].set_title("top view (2x2 repeat)")
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

    slabs = build_fcc_slabs()
    rows: list[dict[str, object]] = []
    for col, (label, _miller) in enumerate(FCC_SURFACES):
        slab = slabs[col]
        rows.append(summarize_slab(label, slab))
        write_and_verify(slab, path_output / f"POSCAR_{label.replace('(', '_').replace(')', '').replace(' ', '')}")

    path_image = path_output / "fcc_surfaces.png"
    render_comparison(slabs, path_image)

    print("label     formula  atoms  a_in_plane_A  b_in_plane_A  c_total_A  pbc")
    for row in rows:
        print(
            f"{row['label']:<8}  {row['formula']:<7}  "
            f"{row['atoms']:<5}  "
            f"{row['a_in_plane_A']:<12.4f}  "
            f"{row['b_in_plane_A']:<12.4f}  "
            f"{row['c_total_A']:<9.4f}  "
            f"{row['pbc']}"
        )
    print("wrote: " + str(path_image))
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build and compare FCC Cu(100)/(110)/(111) slabs.")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("fcc-surfaces-output"),
        help="Directory for POSCAR files and the comparison PNG.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    rows = run_example(parse_args().output)
    formulas = [row["formula"] for row in rows]
    assert all(f == "Cu12" for f in formulas), formulas
    assert all(row["atoms"] == NUM_LAYERS for row in rows)
    assert all(row["pbc"] == [True, True, True] for row in rows)
