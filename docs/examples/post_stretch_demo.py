#!/usr/bin/env python3
"""Demonstrate uniaxial stretch energy curve fitting and equilibrium extraction.

This tutorial is VASP-free. It demonstrates :mod:`mymetal.post.stretch`, which
post-processes a series of uniformly stretched cells to extract the
equilibrium lattice parameter (minimum energy) by fitting the per-atom energy
``E(stretch factor)`` to a quadratic

    E = a*f^2 + b*f + c

and reading off the vertex ``f0 = -b/(2a)``; the equilibrium row vectors are
``rvectors_ref * f0``. This is exactly the math of
:func:`mymetal.post.stretch.post_stretch` /
:func:`mymetal.universal.plot.workflow.my_plot_stretch`.

The real workflow reads energies from ``y_post_data.txt`` scraped by
``pei_vasp_univ_post`` and the stretch factors from the job names; this demo
synthesises the energy-stretch curve analytically from a known equilibrium
factor and noise, runs the same ``np.polyfit`` + vertex extraction, and
reports the recovered equilibrium lattice parameter. No external calculation
is required.
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
# Reference "true" parameters used to synthesise the stretch-energy curve.
# These are *not* DFT data; they are the analytic inputs so the fitted
# equilibrium can be checked against them.
# ---------------------------------------------------------------------------
A_REF = 4.05          # equilibrium lattice parameter, angstrom (FCC Al-like)
F0_REF = 1.0          # equilibrium stretch factor (stretch is relative to ref)
E0_REF_EV = -3.74     # equilibrium energy per atom, eV
CURVATURE = 4.5       # quadratic curvature a (eV per atom per unit stretch^2)

# Stretch-factor sampling grid and synthetic-noise parameters.
F_GRID = np.array([0.94, 0.96, 0.98, 1.00, 1.02, 1.04, 1.06])
NOISE_SEED = 42
NOISE_SIGMA_EV = 1.5e-3   # tiny, well below equilibrium-precision


def build_synthetic_data() -> dict[str, np.ndarray]:
    """Return a synthetic stretch-factor / per-atom-energy table.

    The energy is built analytically as ``E = a*(f - f0)^2 + E0`` from the
    reference curvature and equilibrium factor, plus a reproducible small
    noise so the quadratic fit is exercised on non-exact data, mirroring real
    DFT. The returned dict carries ``f`` (stretch factors), ``e`` (per-atom
    energies in eV), and ``lca`` (c/a ratios, held constant for a uniaxial a
    stretch to mirror the second panel of ``my_plot_stretch``).
    """
    rng = np.random.default_rng(NOISE_SEED)
    f = F_GRID.copy()
    e = CURVATURE * (f - F0_REF) ** 2 + E0_REF_EV
    e = e + rng.normal(0.0, NOISE_SIGMA_EV, size=f.shape)
    # c/a is constant for a pure a-axis stretch of a cubic cell (c = a0)
    lca = np.full_like(f, 1.0)
    return {"f": f, "e": e, "lca": lca}


def fit_stretch_equilibrium(
        data: dict[str, np.ndarray]) -> dict[str, float]:
    """Fit E(f) to a quadratic and extract the equilibrium factor and energy.

    This reproduces the fit + vertex extraction of
    :func:`mymetal.universal.plot.workflow.my_plot_stretch`: it calls
    ``np.polyfit(f, E/natoms, 2)`` (here energies are already per-atom), builds
    a ``np.poly1d``, and reads the vertex ``f0 = -b/(2a)``, ``E0 = p(f0)``.
    """
    f = np.asarray(data["f"], dtype=float)
    e = np.asarray(data["e"], dtype=float)
    coeffs = np.polyfit(f, e, 2)     # highest power first: [a, b, c]
    p = np.poly1d(coeffs)
    a, b, c = coeffs
    f0 = -b / (2.0 * a)
    e0 = float(p(f0))
    return {"coeffs": coeffs, "f0": float(f0), "e0": e0}


def render_figure(data: dict[str, np.ndarray],
                  fit: dict[str, float],
                  path_image: Path) -> None:
    """Render a 2-panel figure: energy-stretch + quadratic fit, and c/a ratio."""
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
        lines_linewidth=1.4,
        lines_markeredgewidth=1.0,
        patch_linewidth=1.2,
        xtick_major_width=1.2,
        xtick_minor_width=0.8,
        ytick_major_width=1.2,
        ytick_minor_width=0.8,
    )

    fig, (ax_e, ax_ca) = plt.subplots(1, 2)

    f = np.asarray(data["f"], dtype=float)
    e = np.asarray(data["e"], dtype=float)
    coeffs = fit["coeffs"]
    f0 = fit["f0"]
    e0 = fit["e0"]

    # ---- Panel A: per-atom energy vs stretch factor with quadratic fit ----
    f_fine = np.linspace(f.min(), f.max(), 200)
    e_fine = np.polyval(coeffs, f_fine)
    # shift to meV above equilibrium, mirroring my_plot_stretch
    e_mev = (e - e0) * 1000.0
    e_fit_mev = (e_fine - e0) * 1000.0

    ax_e.plot(f, e_mev, marker="o", linestyle="", color="tab:blue",
              label="synthetic data")
    ax_e.plot(f_fine, e_fit_mev, color="tab:blue", linestyle="--",
              linewidth=1.6, label="quadratic fit")
    ax_e.axvline(f0, color="tab:red", linestyle=":", linewidth=1.6,
                 label="equilibrium f0")
    ax_e.set_xlabel("stretch factor")
    ax_e.set_ylabel(r"$E - E_0$ (meV per atom)")
    ax_e.set_title("(A) Uniaxial stretch energy curve")
    textstr = ("f0 = %.4f\nE0 = %.6f eV/atom\na0 = %.4f A"
               % (f0, e0, A_REF * f0))
    ax_e.text(0.05, 0.95, textstr, transform=ax_e.transAxes,
              ha="left", va="top", fontsize=10,
              bbox=dict(boxstyle="round", facecolor="white", alpha=0.5))
    legend = ax_e.legend(loc="upper center", fontsize=9, framealpha=0.9)
    general_modify_legend(legend, linewidth=1.2)

    # ---- Panel B: c/a ratio vs stretch factor ----
    ax_ca.plot(f, data["lca"], marker="o", linestyle="", color="tab:green",
               label="c/a (constant for a-stretch)")
    ax_ca.axvline(f0, color="tab:red", linestyle=":", linewidth=1.6,
                  label="equilibrium f0")
    ax_ca.set_xlabel("stretch factor")
    ax_ca.set_ylabel(r"$c/a$", rotation=0)
    ax_ca.set_title("(B) c/a ratio")
    legend = ax_ca.legend(loc="upper right", fontsize=10, framealpha=0.9)
    general_modify_legend(legend, linewidth=1.2)
    ax_ca.set_ylim(0.8, 1.2)

    fig.savefig(path_image, bbox_inches="tight")
    plt.close(fig)

    # ---- non-blank guard ----
    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError("stretch image was not created: " + str(path_image))
    image_rgb = plt.imread(path_image)[..., :3]
    if float(np.mean(np.abs(image_rgb - 1.0))) < 0.002:
        raise AssertionError("stretch image is effectively blank: " + str(path_image))


def run_example(path_output: Path) -> dict[str, object]:
    """Generate synthetic data, fit the stretch curve, render the figure."""
    path_output = path_output.resolve()
    path_output.mkdir(parents=True, exist_ok=True)

    data = build_synthetic_data()
    fit = fit_stretch_equilibrium(data)

    path_image = path_output / "post_stretch_demo.png"
    render_figure(data, fit, path_image)

    # ---- summary ----
    a0_fit = A_REF * fit["f0"]
    print("=" * 64)
    print("Uniaxial stretch demo (synthetic Al-like data)")
    print("  stretch grid = " + str(F_GRID.tolist()))
    print("  noise sigma  = %.1e eV  (seed=%d)" % (NOISE_SIGMA_EV, NOISE_SEED))
    print("-" * 64)
    print("  reference f0 = %.6f   a0 = %.4f A" % (F0_REF, A_REF))
    print("  fitted    f0 = %.6f   a0 = %.4f A" % (fit["f0"], a0_fit))
    print("  fitted    E0 = %.6f eV/atom" % fit["e0"])
    print("-" * 64)
    print("NOTE: data is synthetic (built from reference parameters);")
    print("      fit reproduces my_plot_stretch polynomial math, not real DFT.")
    print("wrote: " + str(path_image))
    print("=" * 64)

    return {"fit": fit, "a0_fit": a0_fit, "path_image": path_image}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fit uniaxial stretch energy and extract equilibrium (VASP-free).")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("post-stretch-output"),
        help="Directory for the stretch PNG.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    result = run_example(parse_args().output)
    fit = result["fit"]

    # Deterministic assertions: fitted equilibrium within tolerance of inputs.
    assert abs(fit["f0"] - F0_REF) < 0.01, (
        "f0 %.6f not within 0.01 of %.6f" % (fit["f0"], F0_REF))
    assert abs(fit["e0"] - E0_REF_EV) < 0.01, (
        "E0 %.6f not within 0.01 eV of %.6f" % (fit["e0"], E0_REF_EV))
    # Image must exist and be non-blank (also checked inside render_figure).
    assert result["path_image"].is_file()
