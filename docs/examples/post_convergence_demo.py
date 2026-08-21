#!/usr/bin/env python3
"""Demonstrate ENCUT and KPOINTS convergence tests.

This tutorial is VASP-free. It demonstrates :mod:`mymetal.post.convergence`,
which post-processes a series of single-point energy calculations at varying
ENCUT and k-point grids and plots the per-atom energy and its delta from the
converged (highest-parameter) value, mirroring
:func:`mymetal.universal.plot.workflow.my_plot_convergence` with
``if_difference=True``.

The real workflow reads energies from ``y_post_data.txt`` scraped by
``pei_vasp_univ_post``; this demo synthesises the per-atom energies
analytically from a saturating-exponential form (energy approaches the
converged value as ENCUT/kpoint density rises), plus a reproducible small
noise, and reports the delta-from-converged in meV/atom. No external
calculation is required.
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
# Reference "true" converged energies (eV/atom). These are *not* DFT data;
# they are the analytic inputs used to synthesise the convergence curves.
# ---------------------------------------------------------------------------
E_CONV_ENCUT = -3.74000   # converged per-atom energy at high ENCUT, eV
E_CONV_KPTS = -3.74005     # converged per-atom energy at dense k-grid, eV

# ENCUT and kpoint sampling grids.
ENCUT_GRID = np.array([300, 350, 400, 450, 500, 550, 600])
KPOINTS_GRID = np.array([
    [4, 4, 4], [6, 6, 6], [8, 8, 8], [10, 10, 10],
    [12, 12, 12], [14, 14, 14], [16, 16, 16],
])

NOISE_SEED = 42
NOISE_SIGMA_EV = 1e-4   # tiny, well below convergence precision


def build_synthetic_encut() -> dict[str, np.ndarray]:
    """Return synthetic per-atom energies vs ENCUT (eV/atom).

    Energy approaches the converged value as a saturating exponential in
    ENCUT, plus a reproducible small noise. The returned dict carries
    ``encut`` and ``e`` (per-atom energy).
    """
    rng = np.random.default_rng(NOISE_SEED)
    encut = ENCUT_GRID.copy()
    # saturating exponential: E(encut) = E_conv - A * exp(-encut / tau)
    A_encut = 0.020
    tau_encut = 120.0
    e = E_CONV_ENCUT - A_encut * np.exp(-encut / tau_encut)
    e = e + rng.normal(0.0, NOISE_SIGMA_EV, size=encut.shape)
    return {"encut": encut, "e": e, "e_conv": E_CONV_ENCUT}


def build_synthetic_kpoints() -> dict[str, object]:
    """Return synthetic per-atom energies vs k-point grid (eV/atom).

    Energy approaches the converged value as a saturating exponential in the
    total k-point count, plus a reproducible small noise. The returned dict
    carries ``grid`` (list of [kx, ky, kz]) and ``e`` (per-atom energy).
    """
    rng = np.random.default_rng(NOISE_SEED + 1)
    grid = KPOINTS_GRID.copy()
    ktot = np.array([np.prod(g) for g in grid], dtype=float)
    A_kpts = 0.015
    tau_kpts = 400.0
    e = E_CONV_KPTS - A_kpts * np.exp(-ktot / tau_kpts)
    e = e + rng.normal(0.0, NOISE_SIGMA_EV, size=ktot.shape)
    return {"grid": grid, "e": e, "e_conv": E_CONV_KPTS}


def render_figure(encut_data: dict[str, np.ndarray],
                  kpoints_data: dict[str, object],
                  path_image: Path) -> None:
    """Render a 2-panel figure: delta-E vs ENCUT and delta-E vs kpoint grid."""
    path_image.parent.mkdir(parents=True, exist_ok=True)

    general_set_all_rcParams(
        backend="Agg",
        axes_grid=True,
        figure_subp=(1, 2),
        figure_one_figsize=(6.4, 5.0),
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
        lines_markersize=8,
        lines_linewidth=1.6,
        lines_markeredgewidth=1.0,
        patch_linewidth=1.2,
        xtick_major_width=1.2,
        xtick_minor_width=0.8,
        ytick_major_width=1.2,
        ytick_minor_width=0.8,
    )

    fig, (ax_encut, ax_kpts) = plt.subplots(1, 2)

    # ---- Panel A: delta-E vs ENCUT ----
    encut = encut_data["encut"]
    e_encut = encut_data["e"]
    e_conv_encut = encut_data["e_conv"]
    delta_encut = (e_encut - e_conv_encut) * 1000.0   # meV/atom
    ax_encut.plot(encut, delta_encut, marker="o", linestyle="-",
                  color="tab:blue", label="delta from converged")
    ax_encut.axhline(0.0, color="gray", linestyle=":", linewidth=1.2)
    ax_encut.set_xlabel("ENCUT (eV)")
    ax_encut.set_ylabel(r"$E - E_{\mathrm{conv}}$ (meV per atom)")
    ax_encut.set_title("(A) ENCUT convergence")
    # annotate non-converged points (|delta| > 1 meV/atom), mirroring if_mask
    for i, d in enumerate(delta_encut):
        if abs(d) > 1.0:
            ax_encut.annotate("%.2f" % d, (encut[i], d),
                               textcoords="offset points", xytext=(0, 8),
                               fontsize=9, color="tab:red")
    legend = ax_encut.legend(loc="upper right", fontsize=10, framealpha=0.9)
    general_modify_legend(legend, linewidth=1.2)

    # ---- Panel B: delta-E vs kpoint grid ----
    grid = kpoints_data["grid"]
    e_kpts = kpoints_data["e"]
    e_conv_kpts = kpoints_data["e_conv"]
    delta_kpts = (e_kpts - e_conv_kpts) * 1000.0   # meV/atom
    labels = [r"$%d\times%d\times%d$" % (int(g[0]), int(g[1]), int(g[2]))
              for g in grid]
    x = np.arange(len(grid))
    ax_kpts.plot(x, delta_kpts, marker="o", linestyle="-",
                 color="tab:green", label="delta from converged")
    ax_kpts.axhline(0.0, color="gray", linestyle=":", linewidth=1.2)
    ax_kpts.set_xticks(x)
    ax_kpts.set_xticklabels(labels, rotation=45, fontsize=10)
    ax_kpts.set_xlabel("k-point grid")
    ax_kpts.set_ylabel(r"$E - E_{\mathrm{conv}}$ (meV per atom)")
    ax_kpts.set_title("(B) KPOINTS convergence")
    for i, d in enumerate(delta_kpts):
        if abs(d) > 1.0:
            ax_kpts.annotate("%.2f" % d, (x[i], d),
                              textcoords="offset points", xytext=(0, 8),
                              fontsize=9, color="tab:red")
    legend = ax_kpts.legend(loc="upper right", fontsize=10, framealpha=0.9)
    general_modify_legend(legend, linewidth=1.2)

    fig.savefig(path_image, bbox_inches="tight")
    plt.close(fig)

    # ---- non-blank guard ----
    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError("convergence image was not created: " + str(path_image))
    image_rgb = plt.imread(path_image)[..., :3]
    if float(np.mean(np.abs(image_rgb - 1.0))) < 0.002:
        raise AssertionError("convergence image is effectively blank: " + str(path_image))


def run_example(path_output: Path) -> dict[str, object]:
    """Generate synthetic ENCUT/KPOINTS data, render the figure, return summary."""
    path_output = path_output.resolve()
    path_output.mkdir(parents=True, exist_ok=True)

    encut_data = build_synthetic_encut()
    kpoints_data = build_synthetic_kpoints()

    path_image = path_output / "post_convergence_demo.png"
    render_figure(encut_data, kpoints_data, path_image)

    # ---- summary ----
    print("=" * 64)
    print("ENCUT / KPOINTS convergence demo (synthetic data)")
    print("  noise sigma = %.1e eV  (seed=%d)" % (NOISE_SIGMA_EV, NOISE_SEED))
    print("-" * 64)
    print("ENCUT panel:")
    print("  ENCUT grid = " + str(ENCUT_GRID.tolist()))
    print("  delta (meV/atom) = " + str(
        [float("%.3f" % d) for d in
         (encut_data["e"] - encut_data["e_conv"]) * 1000.0]))
    print("KPOINTS panel:")
    print("  grid = " + str([g.tolist() for g in KPOINTS_GRID]))
    print("  delta (meV/atom) = " + str(
        [float("%.3f" % d) for d in
         (kpoints_data["e"] - kpoints_data["e_conv"]) * 1000.0]))
    print("-" * 64)
    print("NOTE: data is synthetic (built from saturating exponentials);")
    print("      delta mirrors my_plot_convergence(if_difference=True), not DFT.")
    print("wrote: " + str(path_image))
    print("=" * 64)

    return {"encut": encut_data, "kpoints": kpoints_data,
            "path_image": path_image}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot synthetic ENCUT/KPOINTS convergence (VASP-free).")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("post-convergence-output"),
        help="Directory for the convergence PNG.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    result = run_example(parse_args().output)

    # Deterministic assertions: the last (converged) point is within 0.1 meV/atom
    # of the converged value by construction (delta ~ noise).
    encut = result["encut"]
    kpts = result["kpoints"]
    delta_encut_last = abs(float((encut["e"][-1] - encut["e_conv"]) * 1000.0))
    delta_kpts_last = abs(float((kpts["e"][-1] - kpts["e_conv"]) * 1000.0))
    assert delta_encut_last < 0.5, (
        "last ENCUT delta %.3f meV/atom >= 0.5" % delta_encut_last)
    assert delta_kpts_last < 0.5, (
        "last KPOINTS delta %.3f meV/atom >= 0.5" % delta_kpts_last)
    # The first (coarsest) point must be more off than the last (converged).
    delta_encut_first = abs(float((encut["e"][0] - encut["e_conv"]) * 1000.0))
    assert delta_encut_first > delta_encut_last
    # Image must exist and be non-blank (also checked inside render_figure).
    assert result["path_image"].is_file()
