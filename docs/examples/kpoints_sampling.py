#!/usr/bin/env python3
"""Generate and visualize k-point grids for VASP calculations.

This tutorial is VASP-free. It demonstrates three commonly-used k-point
schemes provided by :mod:`mymetal.calculate.calqm.kpoints`:

1. :func:`get_kpoints_by_size` -- explicit Monkhorst-Pack vs Gamma-centered
   grids for the same ``size`` and ``offset``.
2. :func:`get_size_by_distance` -- automatic mesh from a target ``RK``
   product (old VASP scheme) or ``KSPACING`` (new VASP scheme), using the
   reciprocal cell of an ``ase.Atoms`` object.
3. :func:`cal_reciprocal_matrix` and :func:`cal_reciprocal_matrix2` --
   real-space cell -> reciprocal lattice vectors, compared on the same
   structure so the two equivalent implementations agree.

The script renders a three-panel figure reused by the documentation page:

* Left: Monkhorst-Pack vs Gamma-centered k-points in the (kx, ky) plane.
* Middle: k-point count per direction as a function of RK.
* Right: RK -> k-point count for a slab-like Cu(111) cell, comparing the
  old (rounding) and new (ceil) VASP recipes.
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

from mymetal.build.film.stretch import generate_film
from mymetal.calculate.calqm.kpoints import (
    cal_reciprocal_matrix,
    cal_reciprocal_matrix2,
    get_kpoints_by_size,
    get_size_by_distance,
)
from mymetal.universal.plot.general import general_set_all_rcParams


# Lattice constant (Å) and grid size used as a teaching baseline. The
# numbers come from the Cu(111) tutorial so the k-point density matches
# what a real surface calculation would use.
A_FCC_CU = 3.61
GRID_SIZE = (6, 6, 1)
GRID_OFFSET = (0.5, 0.5, 0.5)
SLAB_VACUUM = 20.0
SLAB_LAYERS = 4


def build_slab() -> Atoms:
    """Return a 4-layer Cu(111) slab used as the k-point reference cell."""
    return generate_film(
        symbols="Cu",
        structure="fcc",
        num_layers=SLAB_LAYERS,
        my_vacuum=SLAB_VACUUM,
        slice_plane=(1, 1, 1),
        a_fcc=A_FCC_CU,
    )


def get_mp_and_gamma(size: tuple[int, int, int] = GRID_SIZE,
                     offset: tuple[float, float, float] = GRID_OFFSET):
    """Return Monkhorst-Pack and Gamma-centered k-point arrays."""
    mpk, gk = get_kpoints_by_size(size, offset)
    return np.asarray(mpk), np.asarray(gk)


def rk_scan(atoms: Atoms, rk_list: list[int]) -> dict[str, np.ndarray]:
    """Return old and new VASP k-point meshes for each RK in ``rk_list``."""
    old_list = []
    new_list = []
    for rk in rk_list:
        old_k, new_k = get_size_by_distance(atoms=atoms, rk=rk)
        old_list.append(np.asarray(old_k))
        new_list.append(np.asarray(new_k))
    return {
        "rk": np.asarray(rk_list, dtype=int),
        "old": np.asarray(old_list),
        "new": np.asarray(new_list),
    }


def compare_reciprocal(cell_matrix: np.ndarray) -> dict[str, np.ndarray]:
    """Return reciprocal cell from cross product and matrix inversion."""
    b1, b2, b3 = cal_reciprocal_matrix(cell_matrix, scale=2.0 * np.pi)
    cross = np.vstack([b1, b2, b3])
    inv = cal_reciprocal_matrix2(cell_matrix, scale=2.0 * np.pi)
    return {"cross": cross, "inv": inv}


def render_figure(
    mpk: np.ndarray,
    gk: np.ndarray,
    rk_result: dict[str, np.ndarray],
    path_image: Path,
) -> None:
    """Render the three-panel k-point figure."""
    path_image.parent.mkdir(parents=True, exist_ok=True)

    general_set_all_rcParams(
        backend="Agg",
        axes_grid=True,
        figure_subp=(1, 3),
        figure_one_figsize=(4.8, 4.2),
        figure_left=0.65,
        figure_top=0.55,
        figure_wspace=1.1,
        figure_hspace=0.6,
        font_family=["DejaVu Sans"],
        fontsize=11,
        axes_titlepad=6,
    )
    fig, axes = plt.subplots(1, 3)

    # Left: Monkhorst-Pack vs Gamma-centered.
    ax = axes[0]
    ax.plot(mpk[:, 0], mpk[:, 1], "o", color="#1f77b4",
            label="Monkhorst-Pack", markersize=8)
    ax.plot(gk[:, 0], gk[:, 1], "x", color="#d62728",
            label="Gamma-centered", markersize=10, mew=2)
    ax.axhline(0.0, color="gray", linestyle=":", linewidth=0.8)
    ax.axvline(0.0, color="gray", linestyle=":", linewidth=0.8)
    ax.set_xlabel(r"$k_x$ ($2\pi/a$)")
    ax.set_ylabel(r"$k_y$ ($2\pi/a$)")
    ax.set_title(f"MP vs Gamma-centered\nsize = {GRID_SIZE}")
    ax.legend(loc="upper right", fontsize=9)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-0.6, 0.6)
    ax.set_ylim(-0.6, 0.6)

    # Middle: per-direction k-point count vs RK for the slab.
    ax = axes[1]
    rk = rk_result["rk"]
    old = rk_result["old"]
    new = rk_result["new"]
    ax.plot(rk, old[:, 0], "o-", color="#1f77b4", label=r"$n_x$ (old, round)")
    ax.plot(rk, new[:, 0], "x--", color="#1f77b4", label=r"$n_x$ (new, ceil)",
            mew=2)
    ax.plot(rk, old[:, 2], "s-", color="#d62728", label=r"$n_z$ (old, round)")
    ax.plot(rk, new[:, 2], "v--", color="#d62728", label=r"$n_z$ (new, ceil)")
    ax.set_xlabel(r"$R_k$ product")
    ax.set_ylabel("k-point count per direction")
    ax.set_title("RK -> mesh for Cu(111) slab")
    ax.legend(loc="upper left", fontsize=9)

    # Right: ratio between in-plane and out-of-plane density.
    ax = axes[2]
    ratio_old = old[:, 0] / np.maximum(old[:, 2], 1)
    ratio_new = new[:, 0] / np.maximum(new[:, 2], 1)
    ax.plot(rk, ratio_old, "o-", color="#2ca02c", label=r"$n_x / n_z$ (old)")
    ax.plot(rk, ratio_new, "x--", color="#9467bd", label=r"$n_x / n_z$ (new)",
            mew=2)
    ax.set_xlabel(r"$R_k$ product")
    ax.set_ylabel(r"$n_x / n_z$")
    ax.set_title("In-plane vs out-of-plane density")
    ax.legend(loc="lower right", fontsize=9)

    fig.savefig(path_image, bbox_inches="tight")
    plt.close(fig)

    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError("kpoints image was not created: " + str(path_image))


def run_example(path_output: Path):
    """Build the slab, generate k-points, render the figure, return results."""
    path_output = path_output.resolve()
    path_output.mkdir(parents=True, exist_ok=True)

    slab = build_slab()
    mpk, gk = get_mp_and_gamma()
    rk_list = [20, 40, 60, 80, 100, 120]
    rk_result = rk_scan(slab, rk_list)

    cell_matrix = np.asarray(slab.get_cell())
    recip = compare_reciprocal(cell_matrix)

    path_image = path_output / "kpoints_sampling.png"
    render_figure(mpk, gk, rk_result, path_image)

    print(f"slab formula: {slab.get_chemical_formula()}")
    print(f"slab atoms: {len(slab)}")
    print(f"slab cell lengths (A): "
          f"{np.linalg.norm(slab.cell.array, axis=1).round(4).tolist()}")
    print(f"MP k-points (first 5): {mpk[:5].tolist()}")
    print(f"Gamma k-points (first 5): {gk[:5].tolist()}")
    print(f"MP count: {len(mpk)}, Gamma count: {len(gk)}")
    print(f"RK list: {rk_list}")
    print(f"old VASP mesh per RK: {rk_result['old'].tolist()}")
    print(f"new VASP mesh per RK: {rk_result['new'].tolist()}")
    print("reciprocal (cross-product) b1, b2, b3:")
    for vec in recip["cross"]:
        print("  " + ", ".join(f"{v:.6f}" for v in vec))
    print("reciprocal (matrix-inverse) rows:")
    for vec in recip["inv"]:
        print("  " + ", ".join(f"{v:.6f}" for v in vec))
    print("wrote: " + str(path_image))
    return {
        "mpk": mpk,
        "gk": gk,
        "rk_result": rk_result,
        "recip": recip,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate and visualize k-point grids for a Cu(111) slab.")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("kpoints-output"),
        help="Directory for the k-points PNG.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    result = run_example(parse_args().output)
    # Both grids must have the same total number of k-points.
    assert len(result["mpk"]) == len(result["gk"])
    # Gamma-centered grid must differ from MP grid when the size is even.
    assert not np.allclose(result["mpk"], result["gk"])
    # The two reciprocal-lattice implementations must agree to 1e-10.
    assert np.allclose(result["recip"]["cross"], result["recip"]["inv"],
                       atol=1e-10)
    # In-plane k-point count must grow monotonically with RK.
    old_nx = result["rk_result"]["old"][:, 0]
    assert np.all(np.diff(old_nx) >= 0)
    # Out-of-plane count must stay small because of the 20 Å vacuum.
    assert np.all(result["rk_result"]["old"][:, 2] <= 3)
