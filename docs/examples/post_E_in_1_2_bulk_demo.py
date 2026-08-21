#!/usr/bin/env python3
"""Demonstrate the E_in_1/2 bulk (biaxial) energy landscape.

This tutorial is VASP-free. It demonstrates :mod:`mymetal.post.E_in_1_2_bulk`,
which post-processes a grid of biaxially (a1, a2) strained bulk cells to map
the per-atom energy surface, locate the equilibrium (a1, a2) at the energy
minimum, and render a 2D contour plus profile plots, mirroring
:func:`mymetal.universal.plot.workflow.my_plot_E_in_1_2_bulk`.

The real workflow reads energies from ``y_post_data.txt`` scraped by
``pei_vasp_univ_post`` and the (a1, a2) labels from the job names; this demo
synthesises the energy surface analytically as a 2D quadratic well around a
known equilibrium (plus a reproducible small noise), runs the same grid +
equilibrium extraction the library uses, and draws a 2D contour with the
equilibrium marker and two profile cuts (E vs a1 at a2_eq, E vs a2 at a1_eq).
No external calculation is required.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mymetal.universal.plot.general import (
    general_modify_legend,
    general_set_all_rcParams,
)


# ---------------------------------------------------------------------------
# Reference "true" equilibrium parameters. These are *not* DFT data; they are
# the analytic inputs used to synthesise the 2D energy surface.
# ---------------------------------------------------------------------------
A1_EQ = 2.90        # equilibrium a1, angstrom (HCP Mg-like a)
A2_EQ = 2.90        # equilibrium a2, angstrom
E0_EV = -1.29       # equilibrium energy per atom, eV (HCP Mg-like)
CURV1 = 8.0         # curvature along a1 (eV/A^2)
CURV2 = 8.0         # curvature along a2 (eV/A^2)
CROSS = 1.5         # a1-a2 coupling (eV/A^2)

# Grid sampling around the equilibrium.
N_GRID = 7
SPAN = 0.10         # +/- 0.10 A around equilibrium

NOISE_SEED = 42
NOISE_SIGMA_EV = 1e-4   # eV


def build_synthetic_data() -> dict[str, object]:
    """Return a synthetic (a1, a2, energy-per-atom) grid for the biaxial well.

    The energy surface is built analytically as a 2D quadratic well

        E(a1, a2) = E0 + 0.5*C1*(a1 - a1_eq)^2 + 0.5*C2*(a2 - a2_eq)^2
                    + C12*(a1 - a1_eq)*(a2 - a2_eq)

    over a regular grid around the equilibrium, plus a reproducible small noise
    so the equilibrium extraction is exercised on non-exact data, mirroring
    real DFT. The returned dict carries ``la1``, ``la2``, ``lenergy``
    (per-atom energy in eV), and ``eq`` (a1_eq, a2_eq, E0).
    """
    rng = np.random.default_rng(NOISE_SEED)
    a1_vals = np.linspace(A1_EQ - SPAN, A1_EQ + SPAN, N_GRID)
    a2_vals = np.linspace(A2_EQ - SPAN, A2_EQ + SPAN, N_GRID)
    la1: list[float] = []
    la2: list[float] = []
    lenergy: list[float] = []
    for a2 in a2_vals:
        for a1 in a1_vals:
            e = (E0_EV
                 + 0.5 * CURV1 * (a1 - A1_EQ) ** 2
                 + 0.5 * CURV2 * (a2 - A2_EQ) ** 2
                 + CROSS * (a1 - A1_EQ) * (a2 - A2_EQ))
            e = e + rng.normal(0.0, NOISE_SIGMA_EV)
            la1.append(float(a1))
            la2.append(float(a2))
            la1[-1] = round(la1[-1], 6)
            la2[-1] = round(la2[-1], 6)
            lenergy.append(float(e))
    return {"la1": la1, "la2": la2, "lenergy": lenergy,
            "a1_vals": a1_vals, "a2_vals": a2_vals,
            "eq": (A1_EQ, A2_EQ, E0_EV)}


def find_equilibrium(data: dict[str, object]) -> tuple[float, float, float]:
    """Locate the (a1, a2) minimum of the energy surface.

    This reproduces the equilibrium-extraction step of
    :func:`mymetal.universal.plot.workflow.my_plot_E_in_1_2_bulk`: that function
    gridding finds the single ``energy_grid == 0`` (minimum) point. Here we
    find the minimum over the flattened grid.
    """
    la1 = data["la1"]
    la2 = data["la2"]
    lenergy = data["lenergy"]
    i_min = int(np.argmin(lenergy))
    return float(la1[i_min]), float(la2[i_min]), float(lenergy[i_min])


def render_figure(data: dict[str, object],
                  eq_fit: tuple[float, float, float],
                  path_image: Path) -> None:
    """Render a 2-panel figure: 2D contour + equilibrium, and profile cuts."""
    path_image.parent.mkdir(parents=True, exist_ok=True)

    general_set_all_rcParams(
        backend="Agg",
        axes_grid=True,
        figure_subp=(1, 2),
        figure_one_figsize=(6.4, 5.4),
        figure_left=0.75,
        figure_top=0.60,
        figure_wspace=0.90,
        figure_hspace=0.60,
        font_family=["DejaVu Sans"],
        fontsize=12,
        axes_linewidth=1.2,
        axes_titlepad=8,
        grid_linewidth=0.8,
        legend_linewidth=1.2,
        lines_markersize=7,
        lines_linewidth=1.6,
        lines_markeredgewidth=1.0,
        patch_linewidth=1.2,
        xtick_major_width=1.2,
        xtick_minor_width=0.8,
        ytick_major_width=1.2,
        ytick_minor_width=0.8,
    )

    fig, (ax_contour, ax_profile) = plt.subplots(1, 2)

    a1_vals = data["a1_vals"]
    a2_vals = data["a2_vals"]
    la1 = data["la1"]
    la2 = data["la2"]
    lenergy = data["lenergy"]
    eq_a1, eq_a2, eq_energy = eq_fit

    # ---- Panel A: 2D contour of E vs a1, a2 ----
    a1_grid, a2_grid = np.meshgrid(a1_vals, a2_vals)
    energy_grid = np.full(a1_grid.shape, np.nan)
    for i in range(a1_grid.shape[0]):
        for j in range(a1_grid.shape[1]):
            a1 = float(a1_grid[i, j])
            a2 = float(a2_grid[i, j])
            for k in range(len(la1)):
                if abs(la1[k] - a1) < 1e-9 and abs(la2[k] - a2) < 1e-9:
                    energy_grid[i, j] = (lenergy[k] - eq_energy) * 1000.0
                    break

    contourf = ax_contour.contourf(a1_grid, a2_grid, energy_grid, 30,
                                    cmap="coolwarm")
    ax_contour.scatter([eq_a1], [eq_a2], color="red", marker="o", s=120,
                       edgecolor="black", zorder=5, label="equilibrium")
    ax_contour.set_xlabel(r"$a_1$ (Angstrom)")
    ax_contour.set_ylabel(r"$a_2$ (Angstrom)")
    ax_contour.set_title("(A) E(a1, a2) contour")
    legend = ax_contour.legend(loc="upper right", fontsize=10,
                               framealpha=0.9)
    general_modify_legend(legend, linewidth=1.2)
    # equilibrium annotation
    ax_contour.text(0.05, 0.95,
                    "eq:\na1=%.4f\na2=%.4f\nE=%.6f eV"
                    % (eq_a1, eq_a2, eq_energy),
                    transform=ax_contour.transAxes, ha="left", va="top",
                    fontsize=9,
                    bbox=dict(boxstyle="round", facecolor="white",
                              alpha=0.5))
    cbar = fig.colorbar(contourf, ax=ax_contour, label=r"$E - E_0$ (meV/atom)")

    # ---- Panel B: profile cuts at the equilibrium a2 (vs a1) and a1 (vs a2) ----
    # E vs a1 at a2 ~ eq_a2
    i_eq_a2 = int(np.argmin(np.abs(a2_vals - eq_a2)))
    a1_profile = a1_grid[i_eq_a2, :]
    e_a1_profile = energy_grid[i_eq_a2, :]
    ax_profile.plot(a1_profile, e_a1_profile, marker="o", linestyle="-",
                    color="tab:blue", label="vs a1 (at a2_eq)")
    # E vs a2 at a1 ~ eq_a1
    j_eq_a1 = int(np.argmin(np.abs(a1_vals - eq_a1)))
    a2_profile = a2_grid[:, j_eq_a1]
    e_a2_profile = energy_grid[:, j_eq_a1]
    ax_profile.plot(a2_profile, e_a2_profile, marker="s", linestyle="--",
                    color="tab:green", label="vs a2 (at a1_eq)")
    ax_profile.axhline(0.0, color="gray", linestyle=":", linewidth=1.2)
    ax_profile.set_xlabel(r"lattice parameter (Angstrom)")
    ax_profile.set_ylabel(r"$E - E_0$ (meV/atom)")
    ax_profile.set_title("(B) Profiles at equilibrium")
    legend = ax_profile.legend(loc="upper center", fontsize=10,
                                framealpha=0.9)
    general_modify_legend(legend, linewidth=1.2)

    fig.savefig(path_image, bbox_inches="tight")
    plt.close(fig)

    # ---- non-blank guard ----
    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError("E_in_1_2_bulk image was not created: "
                              + str(path_image))
    image_rgb = plt.imread(path_image)[..., :3]
    if float(np.mean(np.abs(image_rgb - 1.0))) < 0.002:
        raise AssertionError("E_in_1_2_bulk image is effectively blank: "
                              + str(path_image))


def run_example(path_output: Path) -> dict[str, object]:
    """Generate synthetic biaxial data, find equilibrium, render the figure."""
    path_output = path_output.resolve()
    path_output.mkdir(parents=True, exist_ok=True)

    data = build_synthetic_data()
    eq_fit = find_equilibrium(data)

    path_image = path_output / "post_E_in_1_2_bulk_demo.png"
    render_figure(data, eq_fit, path_image)

    # ---- summary ----
    print("=" * 64)
    print("E_in_1/2 bulk (biaxial) energy landscape demo (synthetic data)")
    print("  grid = %dx%d   span = +/-%.2f A" % (N_GRID, N_GRID, SPAN))
    print("  noise sigma = %.1e eV  (seed=%d)"
          % (NOISE_SIGMA_EV, NOISE_SEED))
    print("-" * 64)
    print("  reference eq: a1=%.4f  a2=%.4f  E=%.6f eV"
          % (A1_EQ, A2_EQ, E0_EV))
    print("  fitted    eq: a1=%.4f  a2=%.4f  E=%.6f eV"
          % (eq_fit[0], eq_fit[1], eq_fit[2]))
    print("-" * 64)
    print("NOTE: data is synthetic (2D quadratic well + noise);")
    print("      plot mirrors my_plot_E_in_1_2_bulk, not real DFT.")
    print("wrote: " + str(path_image))
    print("=" * 64)

    return {"eq_fit": eq_fit, "path_image": path_image}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot synthetic E_in_1/2 bulk biaxial landscape (VASP-free).")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("post-E-in-1-2-bulk-output"),
        help="Directory for the E_in_1/2 bulk PNG.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    result = run_example(parse_args().output)
    eq_fit = result["eq_fit"]

    # Deterministic assertions: fitted equilibrium within grid spacing of inputs.
    # N_GRID points over +/- SPAN gives a step of 2*SPAN/(N_GRID-1).
    grid_step = 2.0 * SPAN / (N_GRID - 1)
    assert abs(eq_fit[0] - A1_EQ) < grid_step * 1.5, (
        "fitted a1 %.4f not within 1.5*step of %.4f"
        % (eq_fit[0], A1_EQ))
    assert abs(eq_fit[1] - A2_EQ) < grid_step * 1.5, (
        "fitted a2 %.4f not within 1.5*step of %.4f"
        % (eq_fit[1], A2_EQ))
    # Fitted energy must be the grid minimum (<= E0 + noise band).
    assert eq_fit[2] <= E0_EV + 5 * NOISE_SIGMA_EV, (
        "fitted energy %.6f eV not near E0 %.6f" % (eq_fit[2], E0_EV))
    # Image must exist and be non-blank (also checked inside render_figure).
    assert result["path_image"].is_file()
