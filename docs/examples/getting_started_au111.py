#!/usr/bin/env python3
"""Build, write, render, and verify a 12-layer Au(111) slab."""

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


def build_au111_slab() -> Atoms:
    """Return the deterministic structure used throughout Getting Started."""
    return generate_film(
        symbols="Au",
        structure="fcc",
        num_layers=12,
        my_vacuum=20.0,
        slice_plane=(1, 1, 1),
    )


def check_nonblank_image(path_image: Path) -> None:
    """Fail when a generated image is missing or effectively blank."""
    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError("structure image was not created: " + str(path_image))
    image_rgb = plt.imread(path_image)[..., :3]
    if float(np.mean(np.abs(image_rgb - 1.0))) < 0.002:
        raise AssertionError("structure image is effectively blank: " + str(path_image))


def render_slab(
        slab: Atoms,
        path_image: Path,
        title: str = "12-layer Au(111) slab") -> None:
    """Render side and top views with the periodic cell and vacuum visible."""
    path_image.parent.mkdir(parents=True, exist_ok=True)

    # A two-panel layout makes both the long vacuum direction and the compact
    # in-plane lattice readable; DejaVu Sans is available on CI and Windows.
    general_set_all_rcParams(
        backend="Agg",
        axes_grid=False,
        figure_subp=(1, 2),
        figure_one_figsize=(4.8, 4.4),
        figure_left=0.55,
        figure_top=0.45,
        figure_wspace=1.15,
        figure_hspace=1.0,
        font_family=["DejaVu Sans"],
        fontsize=12,
        axes_titlepad=6,
    )
    fig, axes = plt.subplots(1, 2)
    plot_atoms(
        slab.repeat((6, 1, 1)),
        axes[0],
        rotation="90x,0y,0z",
        show_unit_cell=0,
        radii=0.7,
    )
    xmin, xmax = axes[0].get_xlim()
    ymin, ymax = axes[0].get_ylim()
    vacuum = 20.0
    axes[0].set_ylim(ymin - vacuum, ymax + vacuum)
    axes[0].add_patch(plt.Rectangle(
        (xmin, ymin - vacuum),
        xmax - xmin,
        ymax - ymin + 2 * vacuum,
        fill=False,
        linestyle="--",
        linewidth=1.2,
        color="black",
    ))
    axes[0].text(
        (xmin + xmax) / 2,
        ymax + vacuum / 2,
        "vacuum",
        ha="center",
        va="center",
    )
    plot_atoms(
        slab.repeat((2, 2, 1)),
        axes[1],
        rotation="0x,0y,0z",
        show_unit_cell=1,
        radii=0.7,
    )
    axes[0].set_title(title + "\nSide view: 20 Å vacuum on each side")
    axes[1].set_title("Top view: Au(111)")
    for ax in axes:
        ax.set_axis_off()
    fig.savefig(path_image, bbox_inches="tight")
    plt.close(fig)
    check_nonblank_image(path_image)


def write_and_verify(slab: Atoms, path_poscar: Path) -> None:
    """Write POSCAR and verify a round trip through the project I/O helpers."""
    my_write_vasp(path_poscar, slab, label="Au(111) 12-layer slab")
    slab_read, lattice_scale_factor = my_read_vasp(path_poscar)

    assert lattice_scale_factor == 1.0
    assert slab_read.get_chemical_formula() == "Au12"
    assert len(slab_read) == 12
    assert np.array_equal(slab_read.pbc, slab.pbc)
    assert np.allclose(slab_read.cell.array, slab.cell.array, atol=1e-10)
    assert np.allclose(slab_read.positions, slab.positions, atol=1e-10)


def run_example(path_output: Path) -> None:
    """Run the complete no-VASP example and print its observable results."""
    path_output = path_output.resolve()
    path_output.mkdir(parents=True, exist_ok=True)
    path_poscar = path_output / "POSCAR"
    path_image = path_output / "au111_slab.png"

    slab = build_au111_slab()
    assert slab.get_chemical_formula() == "Au12"
    assert len(slab) == 12
    assert slab.pbc.tolist() == [True, True, True]

    write_and_verify(slab, path_poscar)
    render_slab(slab, path_image)

    print("formula: " + slab.get_chemical_formula())
    print("atoms: " + str(len(slab)))
    print("cell (angstrom):")
    for vector in slab.cell.array:
        print("  [" + ", ".join(f"{value:.6f}" for value in vector) + "]")
    print("pbc: " + str(slab.pbc.tolist()))
    print("round-trip: ok")
    print("image-check: nonblank")
    print("wrote: " + str(path_poscar))
    print("wrote: " + str(path_image))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build and verify the documentation Au(111) slab.")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("au111-output"),
        help="Directory for POSCAR and au111_slab.png.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    run_example(parse_args().output)
