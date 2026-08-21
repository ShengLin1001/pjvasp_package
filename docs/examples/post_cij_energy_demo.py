#!/usr/bin/env python3
"""Demonstrate the Cij energy-strain method on synthetic HCP-like data.

This tutorial is VASP-free. It demonstrates :mod:`mymetal.post.Cij_energy`,
which post-processes LAMMPS/vasp deformation runs to extract the five
independent 2nd-order elastic constants of an HCP crystal
``[C11, C12, C13, C33, C44]`` from the quadratic coefficients of the
strain-energy density ``U/V0`` of five deformation modes.

The real workflow reads energies from ``y_post_data.txt`` scraped by
``pei_vasp_univ_post`` and fits each mode with ``np.polyfit(e, ed, 2)`` via
:func:`mymetal.post.Cij_energy.fit_cij_energy`; this demo synthesises the
strain-energy curves analytically from textbook HCP Mg constants, applies the
same ``np.polyfit`` + coefficient-extraction math, and reports the solved
Cij. No external calculation is required.

Energy-strain method
--------------------
For small strains the elastic energy density obeys

    U / V0 = (1/2) * C_ij * eta_i * eta_j

(sum over Voigt indices i, j). Each deformation mode gives a different
linear combination of the Cij through the leading quadratic coefficient
``lp0`` of ``ed = lp0 * e^2``:

* c11 mode:  lp0 = C11 / 2            ->  C11 = 2 * lp0
* c12 mode:  lp0 = (C11 + C12) / 2    ->  C12 = 2 * lp0 - C11
* c13 mode:  lp0 = (C11 + C13) / 2    ->  C13 = 2 * lp0 - C11
* c33 mode:  lp0 = C33 / 2            ->  C33 = 2 * lp0
* c44 mode:  lp0 = C44 / 2             ->  C44 = 2 * lp0
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
# Reference "true" constants (GPa): textbook ambient HCP Mg. These are *not*
# DFT data; they are the analytic inputs used to synthesise the strain-energy
# curves, so the fitted Cij can be checked against them.
# ---------------------------------------------------------------------------
C11_REF = 59.7
C12_REF = 25.4
C13_REF = 21.5
C33_REF = 65.5
C44_REF = 16.4

# Strain grid and synthetic-noise parameters.
STRAIN_GRID = np.array([-0.02, -0.01, 0.0, 0.01, 0.02])
NOISE_SEED = 42
NOISE_SIGMA_GPA = 1e-5   # tiny, well below Cij precision


def _mode_quadratic_coeff(mode: str) -> float:
    """Return the analytic leading quadratic coefficient (GPa) per mode."""
    if mode == "c11":
        return 0.5 * C11_REF
    if mode == "c12":
        return 0.5 * (C11_REF + C12_REF)
    if mode == "c13":
        return 0.5 * (C11_REF + C13_REF)
    if mode == "c33":
        return 0.5 * C33_REF
    if mode == "c44":
        return 0.5 * C44_REF
    raise ValueError("unknown mode: " + mode)


def build_synthetic_data() -> dict[str, dict[str, np.ndarray]]:
    """Return synthetic strain-energy-density curves for the five HCP modes.

    Each entry has keys ``strain`` (dimensionless e) and ``ed`` (energy
    density in GPa, ``U/V0``) built analytically as ``ed = lp0 * e^2`` from
    the textbook HCP constants, plus a reproducible small noise so the fit
    is exercised on non-exact data, mirroring real DFT.
    """
    rng = np.random.default_rng(NOISE_SEED)
    e = STRAIN_GRID.copy()
    out: dict[str, dict[str, np.ndarray]] = {}
    for mode in ("c11", "c12", "c13", "c33", "c44"):
        lp0 = _mode_quadratic_coeff(mode)
        ed = lp0 * e ** 2
        ed = ed + rng.normal(0.0, NOISE_SIGMA_GPA, size=e.shape)
        out[mode] = {"strain": e.copy(), "ed": ed, "lp0": float(lp0)}
    return out


def fit_cij_from_modes(
        modes: dict[str, dict[str, np.ndarray]]) -> dict[str, float]:
    """Fit quadratic energy-density vs strain and extract the five HCP Cij.

    This reproduces the polynomial-fit step of
    :func:`mymetal.post.Cij_energy.fit_cij_energy`: that function calls
    ``np.polyfit(e, ed, 2)`` on each deformation mode and reads the leading
    coefficient ``lp0 = param[0]`` to build the Cij vector. Here we apply the
    same ``np.polyfit`` to each synthetic mode and invert the mode->Cij
    relations documented in the module docstring.
    """
    lp0: dict[str, float] = {}
    for mode in ("c11", "c12", "c13", "c33", "c44"):
        p = np.polyfit(modes[mode]["strain"], modes[mode]["ed"], 2)
        lp0[mode] = float(p[0])

    c11 = 2.0 * lp0["c11"]
    c12 = 2.0 * lp0["c12"] - c11
    c13 = 2.0 * lp0["c13"] - c11
    c33 = 2.0 * lp0["c33"]
    c44 = 2.0 * lp0["c44"]
    return {"C11": c11, "C12": c12, "C13": c13, "C33": c33, "C44": c44,
            "lp0": lp0}


def render_figure(modes: dict[str, dict[str, np.ndarray]],
                  fitted: dict[str, float],
                  path_image: Path) -> None:
    """Render a 2-panel figure: strain-energy curves + Cij bar comparison."""
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

    fig, (ax_curves, ax_bar) = plt.subplots(1, 2)

    # ---- Panel A: five strain-energy curves with quadratic fits ----
    mode_colors = {"c11": "tab:blue", "c12": "tab:orange",
                   "c13": "tab:green", "c33": "tab:red",
                   "c44": "tab:purple"}
    eta_fine = np.linspace(STRAIN_GRID.min(), STRAIN_GRID.max(), 200)
    for mode in ("c11", "c12", "c13", "c33", "c44"):
        m = modes[mode]
        ax_curves.scatter(m["strain"], m["ed"],
                          color=mode_colors[mode], marker="o",
                          label=mode, zorder=3)
        p = np.polyfit(m["strain"], m["ed"], 2)
        ax_curves.plot(eta_fine, np.polyval(p, eta_fine),
                       color=mode_colors[mode], linestyle="--",
                       linewidth=1.6, zorder=2)
    ax_curves.set_xlabel(r"strain $\eta$")
    ax_curves.set_ylabel(r"energy density $U/V_0$ (GPa)")
    ax_curves.set_title("(A) Five HCP modes (data + quadratic fit)")
    legend = ax_curves.legend(loc="upper center", fontsize=9,
                              framealpha=0.9, ncol=3)
    general_modify_legend(legend, linewidth=1.2)
    ax_curves.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))

    # ---- Panel B: bar chart input vs fitted C11, C12, C13, C33, C44 ----
    labels = ["C11", "C12", "C13", "C33", "C44"]
    input_vals = [C11_REF, C12_REF, C13_REF, C33_REF, C44_REF]
    fitted_vals = [fitted["C11"], fitted["C12"], fitted["C13"],
                   fitted["C33"], fitted["C44"]]
    x = np.arange(len(labels))
    width = 0.36
    ax_bar.bar(x - width / 2, input_vals, width,
               color="tab:gray", edgecolor="black", label="input")
    ax_bar.bar(x + width / 2, fitted_vals, width,
               color="tab:orange", edgecolor="black", label="fitted")
    for xi_i, (vi, vf) in enumerate(zip(input_vals, fitted_vals)):
        ax_bar.text(xi_i - width / 2, vi + 1.2, "%.1f" % vi,
                    ha="center", va="bottom", fontsize=9)
        ax_bar.text(xi_i + width / 2, vf + 1.2, "%.1f" % vf,
                    ha="center", va="bottom", fontsize=9)
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels(labels)
    ax_bar.set_ylabel("elastic constant (GPa)")
    ax_bar.set_title("(B) Input vs fitted HCP Cij (seed=42)")
    legend = ax_bar.legend(loc="upper right", fontsize=10, framealpha=0.9)
    general_modify_legend(legend, linewidth=1.2)
    ax_bar.set_ylim(0, max(max(input_vals), max(fitted_vals)) * 1.18)

    fig.savefig(path_image, bbox_inches="tight")
    plt.close(fig)

    # ---- non-blank guard ----
    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError("Cij image was not created: " + str(path_image))
    image_rgb = plt.imread(path_image)[..., :3]
    if float(np.mean(np.abs(image_rgb - 1.0))) < 0.002:
        raise AssertionError("Cij image is effectively blank: " + str(path_image))


def run_example(path_output: Path) -> dict[str, object]:
    """Generate synthetic data, fit Cij, render figure, return summary rows."""
    path_output = path_output.resolve()
    path_output.mkdir(parents=True, exist_ok=True)

    modes = build_synthetic_data()
    fitted = fit_cij_from_modes(modes)

    path_image = path_output / "post_cij_energy_demo.png"
    render_figure(modes, fitted, path_image)

    # ---- summary table ----
    input_map = {"C11": C11_REF, "C12": C12_REF, "C13": C13_REF,
                 "C33": C33_REF, "C44": C44_REF}
    print("=" * 64)
    print("Cij energy-strain fitting (synthetic HCP Mg-like data)")
    print("  strain grid = " + str(STRAIN_GRID.tolist()))
    print("  noise sigma = %.1e GPa  (seed=%d)" % (NOISE_SIGMA_GPA, NOISE_SEED))
    print("-" * 64)
    print("%10s %12s %13s %12s" % ("constant", "input(GPa)", "fitted(GPa)",
                                   "rel.err(%)"))
    rows: list[dict[str, object]] = []
    for key in ("C11", "C12", "C13", "C33", "C44"):
        vi = input_map[key]
        vf = fitted[key]
        rel = abs(vf - vi) / vi * 100.0
        print("%10s %12.3f %13.3f %12.4f" % (key, vi, vf, rel))
        rows.append({"constant": key, "input": vi, "fitted": vf,
                     "rel_err_pct": rel})
    print("-" * 64)
    print("NOTE: data is synthetic (generated from the input constants),")
    print("      not real DFT. Fit reproduces fit_cij_energy polynomial math.")
    print("wrote: " + str(path_image))
    print("=" * 64)

    return {"rows": rows, "fitted": fitted, "path_image": path_image}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fit HCP Cij from synthetic energy-strain data (VASP-free).")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("post-cij-energy-output"),
        help="Directory for the Cij fitting PNG.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    result = run_example(parse_args().output)
    rows = result["rows"]

    # Deterministic assertions: fitted constants within 5% of inputs.
    for row in rows:
        assert row["rel_err_pct"] < 5.0, (
            "%s relative error %.3f%% >= 5%%" % (row["constant"],
                                                 row["rel_err_pct"]))
    # Sanity: all five constants recovered.
    assert len(rows) == 5
    # Image must exist and be non-blank (also checked inside render_figure).
    assert result["path_image"].is_file()
