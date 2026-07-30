#!/usr/bin/env python3
"""Compute deformation and strain matrices for simple cell deformations.

This tutorial is VASP-free. It demonstrates the mechanics helpers in
:mod:`mymetal.calculate.calmechanics`:

* :func:`cal_deform_matrix` -- the deformation gradient ``F`` between an
  initial cell and a deformed cell, computed as ``F = B * A^-1`` where
  columns of ``A`` and ``B`` are the initial and deformed lattice
  vectors.
* :func:`cal_strain_matrix` -- the Lagrangian (``E = 1/2 (F^T F - I)``)
  and Eulerian (``e = 1/2 (I - (F F^T)^{-1})``) strain tensors from ``F``.
* :func:`cal_principal_and_shear_strain` -- the principal strains
  (eigenvalues), their directions (eigenvectors), and the deviatoric
  shear strain matrix.

The script applies three textbook deformations to a cubic Cu reference
cell, computes the strain tensor for each, and renders a 1x3 panel
figure that visualizes the resulting symmetric strain matrix as a
heatmap with annotations.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from unittest.mock import MagicMock

# ``hetbuilder`` is an optional dependency of
# :mod:`mymetal.calculate.calmechanics.strain` (used only by the
# heterostructure-aware helpers). This tutorial only uses the pure-numpy
# strain-matrix helpers, so we mock the missing import to keep the example
# runnable in a minimal installation.
for _name in ("hetbuilder", "hetbuilder.algorithm", "hetbuilder.plotting"):
    sys.modules.setdefault(_name, MagicMock())

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ase import Atoms
from ase.build import bulk

from mymetal.calculate.calmechanics.deformation import cal_deform_matrix
from mymetal.calculate.calmechanics.strain import (
    cal_principal_and_shear_strain,
    cal_strain_matrix,
)
from mymetal.universal.plot.general import general_set_all_rcParams


A_FCC_CU = 3.61


def build_reference_cell() -> np.ndarray:
    """Return the column-vector lattice matrix of cubic Cu (3x3)."""
    atoms = bulk("Cu", "fcc", a=A_FCC_CU, cubic=True)
    return np.asarray(atoms.cell.array).T


def make_uniaxial_x_strain(cell: np.ndarray, strain: float) -> np.ndarray:
    """Return a cell stretched by ``strain`` along x only."""
    deformed = cell.copy()
    deformed[:, 0] *= (1.0 + strain)
    return deformed


def make_biaxial_xy_strain(cell: np.ndarray, strain: float) -> np.ndarray:
    """Return a cell stretched by ``strain`` along x and y simultaneously."""
    deformed = cell.copy()
    deformed[:, 0] *= (1.0 + strain)
    deformed[:, 1] *= (1.0 + strain)
    return deformed


def make_simple_shear_xy(cell: np.ndarray, gamma: float) -> np.ndarray:
    """Return a cell with simple shear ``gamma`` on the x-y plane."""
    deformed = cell.copy()
    # F = I + gamma * e_x (cross) e_y  -> add gamma to the (0,1) entry
    # in the column-vector lattice matrix means cell[0, 1] += gamma * cell[1, 1]
    # but to keep the deformation textbook-clean we apply F directly to
    # the lattice vectors: B = F * A.
    F = np.eye(3)
    F[0, 1] = gamma
    deformed = F @ cell
    return deformed


def compute_strain_pack(initial: np.ndarray, deformed: np.ndarray) -> dict[str, np.ndarray]:
    """Return the deformation gradient, Lagrangian/Euler strain, principal strains."""
    F = cal_deform_matrix(initial, deformed)
    E_lagrangian, E_eulerian = cal_strain_matrix(F)
    principal_matrix, shear_matrix, eigvals, eigvecs = (
        cal_principal_and_shear_strain(E_lagrangian)
    )
    return {
        "F": F,
        "E_lagrangian": E_lagrangian,
        "E_eulerian": E_eulerian,
        "principal_matrix": principal_matrix,
        "shear_matrix": shear_matrix,
        "eigvals": eigvals,
        "eigvecs": eigvecs,
    }


def render_figure(
    packs: list[dict[str, np.ndarray]],
    titles: list[str],
    path_image: Path,
) -> None:
    """Render the 1x3 strain-matrix heatmap figure."""
    path_image.parent.mkdir(parents=True, exist_ok=True)

    general_set_all_rcParams(
        backend="Agg",
        axes_grid=True,
        figure_subp=(1, 3),
        figure_one_figsize=(4.6, 4.4),
        figure_left=0.60,
        figure_top=0.55,
        figure_wspace=1.05,
        figure_hspace=0.6,
        font_family=["DejaVu Sans"],
        fontsize=11,
        axes_titlepad=10,
    )
    fig, axes = plt.subplots(1, 3)

    vmax = max(np.abs(pack["E_lagrangian"]).max() for pack in packs)
    vmax = float(max(vmax, 1e-6))

    for ax, pack, title in zip(axes, packs, titles):
        matrix = pack["E_lagrangian"]
        # Display the matrix as a heatmap with explicit +0.000 annotations.
        im = ax.imshow(matrix, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_xticks([0, 1, 2])
        ax.set_yticks([0, 1, 2])
        ax.set_xticklabels(["x", "y", "z"])
        ax.set_yticklabels(["x", "y", "z"])
        ax.set_title(title)
        for ii in range(3):
            for jj in range(3):
                value = matrix[ii, jj]
                color = "white" if abs(value) > 0.5 * vmax else "black"
                ax.text(jj, ii, f"{value:+.3f}",
                        ha="center", va="center",
                        color=color, fontsize=10)
        cb = fig.colorbar(im, ax=ax, orientation="vertical",
                          pad=0.10, shrink=0.85)
        cb.set_label("Lagrangian strain")

    fig.suptitle("Lagrangian strain tensors for three reference deformations",
                 y=0.98, fontsize=12)
    fig.savefig(path_image, bbox_inches="tight")
    plt.close(fig)

    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError("strain image was not created: " + str(path_image))


def run_example(path_output: Path):
    """Build the reference cell, deform it, compute strain, render figure."""
    path_output = path_output.resolve()
    path_output.mkdir(parents=True, exist_ok=True)

    cell_ref = build_reference_cell()

    strain_amp = 0.05  # 5% strain, well inside the small-strain regime
    shear_amp = 0.10   # 10% simple shear for visibility

    deformed_uniaxial = make_uniaxial_x_strain(cell_ref, strain_amp)
    deformed_biaxial = make_biaxial_xy_strain(cell_ref, strain_amp)
    deformed_shear = make_simple_shear_xy(cell_ref, shear_amp)

    pack_uniaxial = compute_strain_pack(cell_ref, deformed_uniaxial)
    pack_biaxial = compute_strain_pack(cell_ref, deformed_biaxial)
    pack_shear = compute_strain_pack(cell_ref, deformed_shear)

    packs = [pack_uniaxial, pack_biaxial, pack_shear]
    titles = [
        f"Uniaxial x\nstrain = {strain_amp:+.2f}",
        f"Biaxial xy\nstrain = {strain_amp:+.2f}",
        f"Simple shear xy\ngamma = {shear_amp:+.2f}",
    ]

    path_image = path_output / "strain_deformation.png"
    render_figure(packs, titles, path_image)

    print(f"reference cell (A):")
    for vec in cell_ref.T:
        print("  [" + ", ".join(f"{v:.6f}" for v in vec) + "]")
    for title, pack in zip(titles, packs):
        print(f"\n{title}")
        print("  deformation gradient F:")
        for vec in pack["F"]:
            print("    [" + ", ".join(f"{v:+.6f}" for v in vec) + "]")
        print("  Lagrangian strain E:")
        for vec in pack["E_lagrangian"]:
            print("    [" + ", ".join(f"{v:+.6f}" for v in vec) + "]")
        print("  principal strains (eigenvalues of E):")
        print("    " + str(pack["eigvals"]))
    print("\nwrote: " + str(path_image))
    return {"packs": packs, "titles": titles}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute and visualize deformation / strain matrices.")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("strain-output"),
        help="Directory for the strain PNG.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    result = run_example(parse_args().output)
    pack_uniaxial, pack_biaxial, pack_shear = result["packs"]

    # Uniaxial x strain: E[0, 0] must equal 0.5 * ((1+s)^2 - 1) for s = 0.05.
    s = 0.05
    expected_xx = 0.5 * ((1.0 + s) ** 2 - 1.0)
    assert abs(pack_uniaxial["E_lagrangian"][0, 0] - expected_xx) < 1e-10
    # Off-diagonal entries must vanish for pure uniaxial x strain.
    for ii, jj in [(0, 1), (0, 2), (1, 2)]:
        assert abs(pack_uniaxial["E_lagrangian"][ii, jj]) < 1e-10
    # Biaxial xy strain: E[0, 0] == E[1, 1] and equals expected_xx.
    assert abs(pack_biaxial["E_lagrangian"][0, 0]
               - pack_biaxial["E_lagrangian"][1, 1]) < 1e-12
    assert abs(pack_biaxial["E_lagrangian"][0, 0] - expected_xx) < 1e-10
    # Simple shear: E[0, 1] = E[1, 0] = gamma / 2 (small-strain convention
    # on the Lagrangian tensor); E[0, 0] = gamma^2 / 2 to second order.
    g = 0.10
    assert abs(pack_shear["E_lagrangian"][0, 1] - g / 2.0) < 1e-6
    assert abs(pack_shear["E_lagrangian"][1, 0] - g / 2.0) < 1e-6
    # All strain matrices must be symmetric.
    for pack in result["packs"]:
        assert np.allclose(pack["E_lagrangian"],
                           pack["E_lagrangian"].T, atol=1e-12)
