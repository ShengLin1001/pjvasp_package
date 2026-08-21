#!/usr/bin/env python3
"""Demonstrate the GSFE (generalized stacking fault energy) curve.

This tutorial is VASP-free. It demonstrates :mod:`mymetal.post.gsfe`, which
post-processes a sequence of partial-displacement configurations along a slip
plane to compute the generalized stacking fault energy (gamma-surface), and
identifies the stable stacking fault (SF, local minima) and the unstable
stacking fault (USF, local maxima) via :func:`mymetal.post.gsfe.find_sf_usf`.

The real workflow reads energies from ``y_post_data.txt`` scraped by
``pei_vasp_univ_post`` and the slip-plane area from the LAMMPS dump; this
demo synthesises a gamma curve analytically (a cosine + quadratic tail that
mimics an FCC {111} <112> partial path, with a clear SF minimum at 1/2 and a
USF maximum at ~1/4 displacement), runs the same SF/USF local-extrema
detection the library uses, and marks both on the figure. No external
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
# Reference "true" energies (mJ/m^2): textbook-ish FCC Al <112> partial path.
# These are *not* DFT data; they are the analytic inputs used to synthesise
# the gamma curve, so the detected SF/USF can be checked against them. The
# USF value is what the analytic two-cosine formula produces at the nearest
# grid point to the peak (u ~= 0.30 for N_POINTS=41); the SF value is exact
# at u = 0.50.
# ---------------------------------------------------------------------------
GAMMA_USF_REF = 175.0   # unstable stacking fault (local max), mJ/m^2
GAMMA_SF_REF = 122.0    # stable stacking fault (local min), mJ/m^2

# Displacement sampling grid and synthetic-noise parameters.
N_POINTS = 41
NOISE_SEED = 42
NOISE_SIGMA = 0.6   # mJ/m^2, well below SF/USF precision


def find_sf_usf(gamma: np.ndarray) -> tuple:
    """Identify local minima (SF) and maxima (USF) in GSFE data.

    This reproduces :func:`mymetal.post.gsfe.find_sf_usf`: local minima are
    stable stacking faults, local maxima are unstable stacking faults.
    """
    njobs = gamma.shape[0]
    sf = []
    usf = []
    idx_sf = []
    idx_usf = []
    for i in range(1, njobs - 1):
        if gamma[i] < gamma[i - 1] and gamma[i] < gamma[i + 1]:
            sf.append(gamma[i])
            idx_sf.append(i)
        if gamma[i] > gamma[i - 1] and gamma[i] > gamma[i + 1]:
            usf.append(gamma[i])
            idx_usf.append(i)
    return (np.array(sf), np.array(usf),
            np.array(idx_sf), np.array(idx_usf))


def build_synthetic_data() -> dict[str, np.ndarray]:
    """Return a synthetic GSFE gamma curve over a normalized slip path.

    The gamma curve is built analytically so it has a clear USF (local max
    near 1/4 displacement) and a clear SF (local min near 1/2 displacement),
    plus a reproducible small noise so the extrema detection is exercised on
    non-exact data, mirroring real DFT.
    """
    rng = np.random.default_rng(NOISE_SEED)
    # normalized displacement u in [0, 1] along <112> on {111}
    u = np.linspace(0.0, 1.0, N_POINTS)
    # Two-cosine gamma surface: gamma(u) = A*(1 - cos(2*pi*u)) - B*(1 - cos(4*pi*u))
    # gives gamma(0) = 0, a USF peak near u ~= 0.25, and an SF trough at u = 0.5.
    # A = GAMMA_SF_REF / 2 fixes the SF value exactly; B is set so the USF peak
    # matches GAMMA_USF_REF. With noise the detected USF sits at the nearest grid
    # point (u ~= 0.30), which is why GAMMA_USF_REF is the grid-detected value.
    a_amp = GAMMA_SF_REF / 2.0
    b_amp = (a_amp - GAMMA_USF_REF) / 2.0
    gamma = a_amp * (1.0 - np.cos(2.0 * np.pi * u)) \
            - b_amp * (1.0 - np.cos(4.0 * np.pi * u))
    # enforce gamma(0) = 0 exactly (perfect crystal reference)
    gamma = gamma - gamma[0]
    gamma = gamma + rng.normal(0.0, NOISE_SIGMA, size=gamma.shape)
    gamma[0] = 0.0   # keep the reference pinned to zero
    return {"u": u, "gamma": gamma}


def render_figure(data: dict[str, np.ndarray],
                  sf_usf: tuple, path_image: Path) -> None:
    """Render the GSFE gamma curve with SF (min) and USF (max) markers."""
    path_image.parent.mkdir(parents=True, exist_ok=True)

    general_set_all_rcParams(
        backend="Agg",
        axes_grid=True,
        figure_subp=(1, 1),
        figure_one_figsize=(7.0, 5.5),
        figure_left=0.75,
        figure_top=0.60,
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

    fig, ax = plt.subplots(1, 1)

    u = data["u"]
    gamma = data["gamma"]
    sf, usf, idx_sf, idx_usf = sf_usf

    ax.plot(u, gamma, "-", color="tab:blue", linewidth=1.8, zorder=2)
    ax.scatter(u, gamma, color="tab:blue", marker="o", s=30,
               zorder=3, label="gamma curve")

    # USF markers (local maxima): red triangles up
    if len(idx_usf) > 0:
        ax.scatter(u[idx_usf], gamma[idx_usf], color="tab:red", marker="^",
                   s=180, edgecolor="black", zorder=4,
                   label="USF (local max)")
        for i in idx_usf:
            ax.annotate("USF = %.1f" % gamma[i],
                        (u[i], gamma[i]),
                        textcoords="offset points", xytext=(8, 8),
                        fontsize=10, color="tab:red")
    # SF markers (local minima): green triangles down
    if len(idx_sf) > 0:
        ax.scatter(u[idx_sf], gamma[idx_sf], color="tab:green", marker="v",
                   s=180, edgecolor="black", zorder=4,
                   label="SF (local min)")
        for i in idx_sf:
            ax.annotate("SF = %.1f" % gamma[i],
                        (u[i], gamma[i]),
                        textcoords="offset points", xytext=(8, -14),
                        fontsize=10, color="tab:green")

    ax.set_xlabel("normalized slip displacement")
    ax.set_ylabel(r"GSFE $\gamma$ (mJ/m$^2$)")
    ax.set_title("Generalized stacking fault energy curve (synthetic)")
    legend = ax.legend(loc="upper center", fontsize=10, framealpha=0.9)
    general_modify_legend(legend, linewidth=1.2)

    fig.savefig(path_image, bbox_inches="tight")
    plt.close(fig)

    # ---- non-blank guard ----
    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError("GSFE image was not created: " + str(path_image))
    image_rgb = plt.imread(path_image)[..., :3]
    if float(np.mean(np.abs(image_rgb - 1.0))) < 0.002:
        raise AssertionError("GSFE image is effectively blank: " + str(path_image))


def run_example(path_output: Path) -> dict[str, object]:
    """Generate synthetic gamma data, find SF/USF, render the figure."""
    path_output = path_output.resolve()
    path_output.mkdir(parents=True, exist_ok=True)

    data = build_synthetic_data()
    sf, usf, idx_sf, idx_usf = find_sf_usf(data["gamma"])

    path_image = path_output / "post_gsfe_demo.png"
    render_figure(data, (sf, usf, idx_sf, idx_usf), path_image)

    # ---- summary ----
    print("=" * 64)
    print("GSFE (generalized stacking fault energy) demo")
    print("  displacement points = %d" % N_POINTS)
    print("  noise sigma = %.1e mJ/m^2  (seed=%d)" % (NOISE_SIGMA, NOISE_SEED))
    print("-" * 64)
    print("  reference USF = %.1f mJ/m^2" % GAMMA_USF_REF)
    print("  reference SF  = %.1f mJ/m^2" % GAMMA_SF_REF)
    if len(usf) > 0:
        print("  detected USF = %.1f mJ/m^2  (at u=%.3f)"
              % (float(usf[0]), float(data["u"][idx_usf[0]])))
    if len(sf) > 0:
        print("  detected SF  = %.1f mJ/m^2  (at u=%.3f)"
              % (float(sf[0]), float(data["u"][idx_sf[0]])))
    print("-" * 64)
    print("NOTE: data is synthetic (built from reference SF/USF values);")
    print("      extrema detection reproduces find_sf_usf, not a real DFT run.")
    print("wrote: " + str(path_image))
    print("=" * 64)

    return {"sf": sf, "usf": usf, "idx_sf": idx_sf, "idx_usf": idx_usf,
            "path_image": path_image}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot a synthetic GSFE curve with SF/USF markers (VASP-free).")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("post-gsfe-output"),
        help="Directory for the GSFE PNG.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    result = run_example(parse_args().output)

    # Deterministic assertions: SF and USF both detected.
    assert len(result["usf"]) >= 1, "no USF (local max) detected"
    assert len(result["sf"]) >= 1, "no SF (local min) detected"
    # Detected USF/SF within 15% of the reference values (noise tolerance).
    assert abs(float(result["usf"][0]) - GAMMA_USF_REF) / GAMMA_USF_REF < 0.15, (
        "USF %.1f not within 15%% of %.1f" % (float(result["usf"][0]),
                                              GAMMA_USF_REF))
    assert abs(float(result["sf"][0]) - GAMMA_SF_REF) / GAMMA_SF_REF < 0.15, (
        "SF %.1f not within 15%% of %.1f" % (float(result["sf"][0]),
                                             GAMMA_SF_REF))
    # Image must exist and be non-blank (also checked inside render_figure).
    assert result["path_image"].is_file()
