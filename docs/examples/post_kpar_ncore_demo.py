#!/usr/bin/env python3
"""Demonstrate the KPAR/NCORE parallel benchmark.

This tutorial is VASP-free. It demonstrates :mod:`mymetal.post.kpar_ncore`,
which post-processes a grid of VASP runs at varying KPAR and NCORE to measure
elapsed time and the relative energy delta ``E - E0`` in meV/atom, then plots
both vs NCORE grouped by KPAR, mirroring
:func:`mymetal.universal.plot.workflow.my_plot_kpar_ncore`.

The real workflow reads times from ``y_post_time.txt`` and energies from
``y_post_data.txt`` scraped by ``pei_vasp_univ_post`` via
:func:`mymetal.post.kpar_ncore.read_kpar_ncore_times` /
:func:`mymetal.post.kpar_ncore.read_kpar_ncore_energies`; this demo
synthesises the (KPAR, NCORE) -> time/energy table analytically (time decreases
with KPAR and has a NCORE optimum; energy delta is tiny and noise-like), runs
the same grouping and relative-energy computation the library uses, and draws
the same two-panel figure with a base-2 log NCORE axis. No external
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
# Benchmark grid: KPAR x NCORE. These mirror the KPAR_ORDER / EXPECTED_PAIRS
# table in mymetal.post.kpar_ncore, trimmed to a compact 3x4 subset that still
# shows the typical U-shape (time has a NCORE optimum; energy is flat).
# ---------------------------------------------------------------------------
KPAR_ORDER = [16, 8, 4]
NCORE_VALUES = [1, 2, 4, 8]
NATOMS = 8

# Reference "true" timing model (minutes): time = base / KPAR + a / NCORE + b*NCORE
# (the a/NCORE + b*NCORE term gives the classic NCORE optimum). These are *not*
# real HPC timings; they are analytic inputs so the plotted shape is controlled.
T_BASE_MIN = 40.0
T_A_MIN = 8.0     # serial-ish overhead that shrinks with NCORE
T_B_MIN = 0.15    # communication overhead that grows with NCORE

# Reference converged energy (eV). The energy delta is pure noise around it.
E0_EV = -8.20

NOISE_SEED = 42
NOISE_SIGMA_T_MIN = 0.05   # minutes
NOISE_SIGMA_E_EV = 1e-5    # eV


def build_synthetic_data() -> tuple[dict[tuple[int, int], float],
                                     dict[tuple[int, int], float]]:
    """Return synthetic (KPAR,NCORE) -> elapsed-minutes and total-energy tables.

    Elapsed time follows the analytic model ``T = T_base/KPAR + T_a/NCORE +
    T_b*NCORE`` (a classic U-shape with a NCORE optimum per KPAR), plus a
    reproducible small noise. Energy is the reference ``E0`` plus a tiny
    reproducible noise, so the delta-from-converged is noise-like (a few
    micro-eV/atom), exactly mirroring the real benchmark's flat energy panel.
    """
    rng = np.random.default_rng(NOISE_SEED)
    dict_time: dict[tuple[int, int], float] = {}
    dict_energy: dict[tuple[int, int], float] = {}
    for kpar in KPAR_ORDER:
        for ncore in NCORE_VALUES:
            t = T_BASE_MIN / kpar + T_A_MIN / ncore + T_B_MIN * ncore
            t = t + rng.normal(0.0, NOISE_SIGMA_T_MIN)
            dict_time[(kpar, ncore)] = float(t)
            e = E0_EV + rng.normal(0.0, NOISE_SIGMA_E_EV)
            dict_energy[(kpar, ncore)] = float(e)
    return dict_time, dict_energy


def get_delta_energies(
        dict_energy: dict[tuple[int, int], float],
        natoms: int) -> dict[tuple[int, int], float]:
    """Convert total energies to ``(E - E0)`` in meV/atom.

    This reproduces :func:`mymetal.post.kpar_ncore.get_delta_energies`: the
    minimum measured total energy is used as ``E0``.
    """
    energy0 = min(dict_energy.values())
    return {
        pair: (energy - energy0) * 1000.0 / natoms
        for pair, energy in dict_energy.items()
    }


def render_figure(dict_time: dict[tuple[int, int], float],
                  dict_delta: dict[tuple[int, int], float],
                  path_image: Path) -> None:
    """Render a 2-panel figure: time vs NCORE (by KPAR) and delta-E vs NCORE."""
    path_image.parent.mkdir(parents=True, exist_ok=True)

    general_set_all_rcParams(
        backend="Agg",
        axes_grid=True,
        figure_subp=(1, 2),
        figure_one_figsize=(6.0, 5.0),
        figure_left=0.75,
        figure_top=0.60,
        figure_wspace=0.80,
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

    fig, (ax_t, ax_e) = plt.subplots(1, 2)

    colors = ["tab:blue", "tab:orange", "tab:green"]

    # ---- Panel A: elapsed time vs NCORE, grouped by KPAR ----
    for index, kpar in enumerate(KPAR_ORDER):
        lpairs = sorted(
            [pair for pair in dict_time if pair[0] == kpar],
            key=lambda pair: pair[1])
        if not lpairs:
            continue
        lx = [pair[1] for pair in lpairs]
        ly = [dict_time[pair] for pair in lpairs]
        ax_t.plot(lx, ly, marker="o", color=colors[index % len(colors)],
                  label="KPAR=%d" % kpar)
    ax_t.set_xlabel("NCORE")
    ax_t.set_ylabel("Elapsed time (min)")
    ax_t.set_title("(A) Time vs NCORE (by KPAR)")
    ax_t.set_xticks(NCORE_VALUES)
    ax_t.set_xticklabels([str(n) for n in NCORE_VALUES])
    legend = ax_t.legend(loc="upper right", fontsize=10, framealpha=0.9)
    general_modify_legend(legend, linewidth=1.2)

    # ---- Panel B: delta energy vs NCORE, grouped by KPAR (log-x base-2) ----
    for index, kpar in enumerate(KPAR_ORDER):
        lpairs = sorted(
            [pair for pair in dict_delta if pair[0] == kpar],
            key=lambda pair: pair[1])
        if not lpairs:
            continue
        lx = [pair[1] for pair in lpairs]
        ly = [dict_delta[pair] for pair in lpairs]
        ax_e.plot(lx, ly, marker="o", color=colors[index % len(colors)],
                  label="KPAR=%d" % kpar)
    # base-2 log NCORE axis, mirroring my_plot_kpar_ncore
    ax_e.set_xscale("log", base=2)
    ax_e.set_xticks(NCORE_VALUES)
    ax_e.set_xticklabels([str(n) for n in NCORE_VALUES])
    ax_e.set_xlabel("NCORE")
    ax_e.set_ylabel(r"$\Delta E$ (meV/atom)")
    ax_e.set_title("(B) Energy delta vs NCORE")
    legend = ax_e.legend(loc="upper right", fontsize=10, framealpha=0.9)
    general_modify_legend(legend, linewidth=1.2)

    fig.savefig(path_image, bbox_inches="tight")
    plt.close(fig)

    # ---- non-blank guard ----
    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError("kpar-ncore image was not created: " + str(path_image))
    image_rgb = plt.imread(path_image)[..., :3]
    if float(np.mean(np.abs(image_rgb - 1.0))) < 0.002:
        raise AssertionError("kpar-ncore image is effectively blank: " + str(path_image))


def run_example(path_output: Path) -> dict[str, object]:
    """Generate synthetic KPAR/NCORE data, render the figure, return summary."""
    path_output = path_output.resolve()
    path_output.mkdir(parents=True, exist_ok=True)

    dict_time, dict_energy = build_synthetic_data()
    dict_delta = get_delta_energies(dict_energy, NATOMS)

    path_image = path_output / "post_kpar_ncore_demo.png"
    render_figure(dict_time, dict_delta, path_image)

    # ---- summary ----
    best_pair = min(dict_time, key=dict_time.get)
    e0_measured = min(dict_energy.values())
    print("=" * 64)
    print("KPAR/NCORE benchmark demo (synthetic data)")
    print("  KPAR order  = " + str(KPAR_ORDER))
    print("  NCORE grid  = " + str(NCORE_VALUES))
    print("  natoms      = %d" % NATOMS)
    print("  noise sigma = %.1e min / %.1e eV  (seed=%d)"
          % (NOISE_SIGMA_T_MIN, NOISE_SIGMA_E_EV, NOISE_SEED))
    print("-" * 64)
    print("  best (fastest) = KPAR=%d NCORE=%d  time=%.4f min"
          % (best_pair[0], best_pair[1], dict_time[best_pair]))
    print("  E0 (measured)  = %.10f eV" % e0_measured)
    print("  max deltaE     = %.6f meV/atom" % max(dict_delta.values()))
    print("-" * 64)
    print("NOTE: data is synthetic (analytic timing model + flat energy);")
    print("      plot mirrors my_plot_kpar_ncore, not real HPC timings.")
    print("wrote: " + str(path_image))
    print("=" * 64)

    return {"dict_time": dict_time, "dict_delta": dict_delta,
            "best_pair": best_pair, "path_image": path_image}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot synthetic KPAR/NCORE benchmark (VASP-free).")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("post-kpar-ncore-output"),
        help="Directory for the KPAR/NCORE PNG.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    result = run_example(parse_args().output)

    # Deterministic assertions: best pair has KPAR=16 (highest parallelism) and
    # a non-trivial NCORE (the analytic model has an interior optimum, not
    # NCORE=1, because of the T_a/NCORE + T_b*NCORE trade-off).
    best_pair = result["best_pair"]
    assert best_pair[0] == max(KPAR_ORDER), (
        "best KPAR %d != %d (highest parallelism should be fastest)"
        % (best_pair[0], max(KPAR_ORDER)))
    # All energy deltas are tiny (sub-meV/atom) by construction.
    assert max(result["dict_delta"].values()) < 1.0, (
        "energy delta %.6f meV/atom >= 1.0 (should be noise-like)"
        % max(result["dict_delta"].values()))
    # Image must exist and be non-blank (also checked inside render_figure).
    assert result["path_image"].is_file()
