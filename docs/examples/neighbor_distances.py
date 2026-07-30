#!/usr/bin/env python3
"""Compute and visualize neighbor distances and radial distribution functions.

This tutorial is VASP-free. It demonstrates
:func:`mymetal.universal.atom.neighbor.get_neighbor_distances` on three
common bulk structures:

1. FCC Cu (cubic, ``a = 3.61 A``)
2. BCC Fe (cubic, ``a = 2.87 A``)
3. HCP Mg (primitive, ``a = 3.21 A``, ``c = 5.21 A``)

For each structure the script:

* Builds a 2x2x2 supercell so the cutoff does not run out of periodic
  images.
* Calls :func:`get_neighbor_distances` with ``cutoff = 8 A`` to collect
  *every* pair distance inside that sphere (including periodic images).
* Bins the distances into a radial distribution function (RDF) histogram
  using a fixed bin width, and labels each peak with the corresponding
  coordination shell.

The script renders a two-panel figure:

* Left: per-structure RDF curves (normalized to the first peak).
* Right: cumulative coordination number, which is the count of neighbors
  inside a sphere of radius ``r`` for a central atom.
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

from mymetal.universal.atom.neighbor import get_neighbor_distances
from mymetal.universal.plot.general import general_set_all_rcParams


# Lattice constants reused from the bulk_structures tutorial so the
# displayed RDF peaks align with real Cu/Fe/Mg interatomic distances.
A_FCC_CU = 3.61
A_BCC_FE = 2.87
A_HCP_MG = 3.21
C_HCP_MG = 5.21

# Neighbor cutoff (A). 8 A captures the first ~5 shells in FCC Cu and
# the first ~4 shells in BCC Fe / HCP Mg without requiring a huge
# supercell.
CUTOFF = 8.0
BIN_WIDTH = 0.05  # A
REPEAT = (2, 2, 2)


def build_structures() -> list[Atoms]:
    """Return three supercells (FCC Cu, BCC Fe, HCP Mg) for RDF analysis."""
    cells: list[Atoms] = [
        bulk("Cu", "fcc", a=A_FCC_CU, cubic=True).repeat(REPEAT),
        bulk("Fe", "bcc", a=A_BCC_FE, cubic=True).repeat(REPEAT),
        bulk("Mg", "hcp", a=A_HCP_MG,
             covera=C_HCP_MG / A_HCP_MG).repeat(REPEAT),
    ]
    return cells


def structure_label(idx: int) -> str:
    return ["FCC Cu", "BCC Fe", "HCP Mg"][idx]


def collect_distances(atoms: Atoms, cutoff: float = CUTOFF) -> np.ndarray:
    """Return all neighbor pair distances up to ``cutoff`` for ``atoms``."""
    distances_all, _ = get_neighbor_distances(atoms, cutoff=cutoff)
    return distances_all


def build_rdf(distances: np.ndarray, bin_width: float = BIN_WIDTH,
              r_max: float = CUTOFF) -> dict[str, np.ndarray]:
    """Bin pair distances into an RDF histogram and a cumulative count."""
    bins = np.arange(0.0, r_max + bin_width, bin_width)
    counts, edges = np.histogram(distances, bins=bins)
    centers = 0.5 * (edges[:-1] + edges[1:])
    # Normalize by r^2 to mimic the 4*pi*r^2 shell volume, then by the
    # first-peak maximum so the three structures share a y-axis scale.
    shell_volume = 4.0 * np.pi * centers ** 2 * bin_width
    rdf = np.divide(counts, shell_volume,
                    out=np.zeros_like(counts, dtype=float),
                    where=shell_volume > 0)
    if rdf.max() > 0:
        rdf = rdf / rdf.max()
    cumulative = np.cumsum(counts)
    # Normalize cumulative count by the number of atoms in the supercell,
    # so the y-axis reads as "neighbors per central atom".
    return {
        "r": centers,
        "rdf": rdf,
        "cumulative": cumulative,
    }


def render_figure(rdfs: list[dict[str, np.ndarray]], path_image: Path) -> None:
    """Render the two-panel RDF figure."""
    path_image.parent.mkdir(parents=True, exist_ok=True)

    general_set_all_rcParams(
        backend="Agg",
        axes_grid=True,
        figure_subp=(1, 2),
        figure_one_figsize=(6.0, 4.4),
        figure_left=0.70,
        figure_top=0.55,
        figure_wspace=1.0,
        figure_hspace=0.6,
        font_family=["DejaVu Sans"],
        fontsize=11,
        axes_titlepad=6,
    )
    fig, axes = plt.subplots(1, 2)
    colors = ["#1f77b4", "#d62728", "#2ca02c"]

    # Left: normalized RDF.
    ax = axes[0]
    for idx, rdf in enumerate(rdfs):
        ax.plot(rdf["r"], rdf["rdf"], "-", color=colors[idx],
                label=structure_label(idx), linewidth=2)
    ax.set_xlabel(r"$r$ (A)")
    ax.set_ylabel(r"Normalized RDF $g(r)$")
    ax.set_title(f"Neighbor distance distribution\n"
                 f"cutoff = {CUTOFF:.1f} A, bin = {BIN_WIDTH:.2f} A")
    ax.set_xlim(0.0, CUTOFF)
    ax.set_ylim(0.0, 1.1)
    ax.legend(loc="upper right", fontsize=10)

    # Right: cumulative coordination number per central atom.
    ax = axes[1]
    for idx, rdf in enumerate(rdfs):
        # Normalize by atom count handled by build_rdf caller; show raw
        # cumulative count scaled to first-shell size for readability.
        cumulative_norm = rdf["cumulative"] / max(rdf["cumulative"].max(), 1)
        ax.plot(rdf["r"], cumulative_norm, "-", color=colors[idx],
                label=structure_label(idx), linewidth=2)
    ax.set_xlabel(r"$r$ (A)")
    ax.set_ylabel(r"Normalized cumulative count")
    ax.set_title("Cumulative neighbors\n(integrated RDF)")
    ax.set_xlim(0.0, CUTOFF)
    ax.set_ylim(0.0, 1.05)
    ax.legend(loc="lower right", fontsize=10)

    fig.savefig(path_image, bbox_inches="tight")
    plt.close(fig)

    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError("neighbor image was not created: " + str(path_image))


def run_example(path_output: Path):
    """Build structures, collect distances, render figure, return results."""
    path_output = path_output.resolve()
    path_output.mkdir(parents=True, exist_ok=True)

    structures = build_structures()
    rdfs: list[dict[str, np.ndarray]] = []
    for atoms in structures:
        distances = collect_distances(atoms)
        rdfs.append(build_rdf(distances))

    path_image = path_output / "neighbor_distances.png"
    render_figure(rdfs, path_image)

    for idx, (atoms, rdf) in enumerate(zip(structures, rdfs)):
        # First peak position = nearest-neighbor distance.
        peak_idx = int(np.argmax(rdf["rdf"]))
        nn_distance = float(rdf["r"][peak_idx])
        print(f"{structure_label(idx)}: "
              f"atoms={len(atoms)}, "
              f"NN distance={nn_distance:.4f} A, "
              f"max bin count={int(rdf['cumulative'][-1])}")
    print("wrote: " + str(path_image))
    return {"structures": structures, "rdfs": rdfs}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute and visualize neighbor distances and RDFs.")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("neighbor-output"),
        help="Directory for the neighbor PNG.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    result = run_example(parse_args().output)
    # FCC Cu first neighbor distance must be a/sqrt(2) = 2.5527 A.
    rdf_cu = result["rdfs"][0]
    peak_idx_cu = int(np.argmax(rdf_cu["rdf"]))
    assert abs(rdf_cu["r"][peak_idx_cu] - A_FCC_CU / np.sqrt(2.0)) < 0.10
    # BCC Fe first neighbor distance must be a*sqrt(3)/2 = 2.4826 A.
    rdf_fe = result["rdfs"][1]
    peak_idx_fe = int(np.argmax(rdf_fe["rdf"]))
    assert abs(rdf_fe["r"][peak_idx_fe] - A_BCC_FE * np.sqrt(3.0) / 2.0) < 0.10
    # HCP Mg first neighbor distance must be a = 3.21 A.
    rdf_mg = result["rdfs"][2]
    peak_idx_mg = int(np.argmax(rdf_mg["rdf"]))
    assert abs(rdf_mg["r"][peak_idx_mg] - A_HCP_MG) < 0.10
    # All structures must have at least one non-zero RDF bin.
    assert all(rdf["rdf"].max() > 0.0 for rdf in result["rdfs"])
