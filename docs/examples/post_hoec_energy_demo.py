#!/usr/bin/env python3
"""Demonstrate HOEC (higher-order elastic constants) fitting on synthetic data.

This tutorial is VASP-free. It demonstrates :mod:`mymetal.post.hoec_energy`,
which fits the strain-energy density ``u(xi) = [U(xi) - U(0)] / V0`` (GPa) of
every deformation mode to a polynomial and reads off the xi^2, xi^3, xi^4
coefficients ``P2, P3, P4``. Across all modes these are stacked into the
linear systems

    A2 @ SOEC = P2,   A3 @ TOEC = P3,   A4 @ FOEC = P4

solved by least squares for the 2nd/3rd/4th-order elastic constants.

The real workflow reads energies scraped from OUTCAR by
``pei_vasp_univ_post``; this demo synthesises the strain-energy curves
analytically from textbook cubic Cu constants and reproduces the
``np.polyfit`` + coefficient-extraction step (``fit_P``) directly, so no VASP
run is required. Three representative deformation modes are used: a uniaxial
mode, a shear mode, and a biaxial mode, exactly mirroring the cubic Wang-Li
Table I mode set used by the library function.
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
# Reference "true" constants (GPa): textbook ambient FCC Cu, cubic symmetry.
# These are *not* DFT data; they are the analytic inputs used to synthesise
# the per-mode strain-energy curves, so the fitted P2/P3/P4 can be checked
# against them.
# ---------------------------------------------------------------------------
C11 = 168.0
C12 = 121.0
C44 = 76.0

# 3rd-order (TOEC) and 4th-order (FOEC) reference constants (GPa) for the
# cubic modes used here; taken as representative Wang-Li-like values. The demo
# only checks self-consistency (fitted P_n == analytic P_n), not literature.
C111 = -1861.0
C112 = -845.0
C123 = -585.0
C144 =  -97.0
C166 = -376.0
C456 =   95.0
C1111 = 11460.0
C1112 =  3640.0
C1122 =  3150.0
C1123 =  -990.0
C1144 =   940.0
C1244 =  -420.0
C1456 =   -95.0
C4456 =   350.0

# Strain grid and synthetic-noise parameters.
STRAIN_GRID = np.array([-0.03, -0.02, -0.01, 0.0, 0.01, 0.02, 0.03])
NOISE_SEED = 42
NOISE_SIGMA_GPA = 5e-4   # tiny, well below HOEC precision


def _analytic_P2(mode: str) -> float:
    """Return the analytic xi^2 coefficient (GPa) of a cubic mode."""
    if mode == "uniaxial_x":
        return 0.5 * C11
    if mode == "shear_yz":
        return 0.5 * C44
    if mode == "biaxial_xy":
        return C11 + C12
    raise ValueError("unknown mode: " + mode)


def _analytic_P3(mode: str) -> float:
    """Return the analytic xi^3 coefficient (GPa) of a cubic mode."""
    if mode == "uniaxial_x":
        return (1.0 / 6.0) * C111
    if mode == "shear_yz":
        return (1.0 / 6.0) * C456
    if mode == "biaxial_xy":
        return (1.0 / 6.0) * (C111 + 2.0 * C112 + 2.0 * C123)
    raise ValueError("unknown mode: " + mode)


def _analytic_P4(mode: str) -> float:
    """Return the analytic xi^4 coefficient (GPa) of a cubic mode."""
    if mode == "uniaxial_x":
        return (1.0 / 24.0) * C1111
    if mode == "shear_yz":
        return (1.0 / 24.0) * C4456
    if mode == "biaxial_xy":
        return (1.0 / 24.0) * (C1111 + 4.0 * C1112 + 2.0 * C1122
                               + 4.0 * C1123)
    raise ValueError("unknown mode: " + mode)


def build_synthetic_data() -> dict[str, dict[str, np.ndarray]]:
    """Return synthetic strain-energy-density curves for the three modes.

    Each entry has keys ``strain`` (dimensionless xi) and ``u`` (energy
    density in GPa, ``U/V0``), built as

        u(xi) = P2*xi^2 + P3*xi^3 + P4*xi^4

    from the analytic cubic Wang-Li coefficients, plus a reproducible small
    noise so the polynomial fit is exercised on non-exact data, mirroring
    real DFT. Keys ``p2``, ``p3``, ``p4`` carry the analytic coefficients for
    the self-check.
    """
    rng = np.random.default_rng(NOISE_SEED)
    xi = STRAIN_GRID.copy()
    out: dict[str, dict[str, np.ndarray]] = {}
    for mode in ("uniaxial_x", "shear_yz", "biaxial_xy"):
        p2 = _analytic_P2(mode)
        p3 = _analytic_P3(mode)
        p4 = _analytic_P4(mode)
        u = p2 * xi ** 2 + p3 * xi ** 3 + p4 * xi ** 4
        u = u + rng.normal(0.0, NOISE_SIGMA_GPA, size=xi.shape)
        out[mode] = {
            "strain": xi.copy(),
            "u": u,
            "p2": float(p2),
            "p3": float(p3),
            "p4": float(p4),
        }
    return out


def fit_P_demo(xi: np.ndarray, u: np.ndarray, fitdeg: int = 4) -> tuple:
    """Fit u(xi) to a polynomial and return its low-order coefficients.

    This reproduces the polynomial-fit step of
    :func:`mymetal.post.hoec_energy.fit_P`: that function calls
    ``np.polyfit`` (via ``my_ployfit``) on each deformation mode and packs the
    leading coefficients into ``P0..P4``. Here we call ``np.polyfit`` directly
    so the demo stays self-contained.

    Args:
        xi (np.ndarray): Strain amplitudes.
        u (np.ndarray): Strain-energy densities (GPa).
        fitdeg (int): Polynomial degree (default 4, matching Wang-Li).

    Returns:
        tuple: ``(P0, P1, P2, P3, P4, rms, coeffs)`` with ``coeffs`` the full
        polynomial (highest power first) so the plot can draw the fitted model.
    """
    if len(xi) < fitdeg + 1:
        fitdeg = len(xi) - 1
    coeffs = np.polyfit(xi, u, fitdeg)        # highest power first
    lp = coeffs[::-1]                          # index = power
    P = [float(lp[k]) if k < len(lp) else 0.0 for k in range(5)]
    u_fit = np.polyval(coeffs, xi)
    rms = float(np.sqrt(np.mean((u_fit - u) ** 2)))
    return P[0], P[1], P[2], P[3], P[4], rms, coeffs


def render_figure(modes: dict[str, dict[str, np.ndarray]],
                  dict_P: dict[str, dict[str, float]],
                  path_image: Path) -> None:
    """Render a 2x2 figure: per-mode energy-strain curves + fitted polynomials."""
    path_image.parent.mkdir(parents=True, exist_ok=True)

    general_set_all_rcParams(
        backend="Agg",
        axes_grid=True,
        figure_subp=(2, 2),
        figure_one_figsize=(6.0, 5.0),
        figure_left=0.75,
        figure_top=0.55,
        figure_wspace=0.80,
        figure_hspace=0.80,
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

    fig, axes = plt.subplots(2, 2)
    axes = axes.flatten()

    colors = {"uniaxial_x": "tab:blue",
              "shear_yz": "tab:green",
              "biaxial_xy": "tab:red"}
    labels = {"uniaxial_x": "(a) uniaxial x",
              "shear_yz": "(b) yz shear",
              "biaxial_xy": "(c) biaxial xy"}
    xi_fine = np.linspace(STRAIN_GRID.min(), STRAIN_GRID.max(), 200)

    # ---- Panel 1-3: per-mode energy-strain curves with polynomial fits ----
    for i, mode in enumerate(("uniaxial_x", "shear_yz", "biaxial_xy")):
        ax = axes[i]
        m = modes[mode]
        ax.scatter(m["strain"], m["u"], color=colors[mode], marker="o",
                   zorder=3, label="synthetic data")
        coeffs = dict_P[mode]["coeffs"]
        ax.plot(xi_fine, np.polyval(coeffs, xi_fine),
                color=colors[mode], linestyle="--", linewidth=1.6,
                zorder=2, label="4th-order fit")
        ax.set_xlabel(r"strain $\xi$")
        ax.set_ylabel(r"$u = U/V_0$ (GPa)")
        ax.set_title(labels[mode])
        legend = ax.legend(loc="upper center", fontsize=9, framealpha=0.9)
        general_modify_legend(legend, linewidth=1.2)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))

    # ---- Panel 4: solved constants bar chart (2nd/3rd/4th order) ----
    ax = axes[3]
    p2_in = [_analytic_P2(m) for m in ("uniaxial_x", "shear_yz", "biaxial_xy")]
    p2_fit = [dict_P[m]["P2"] for m in ("uniaxial_x", "shear_yz", "biaxial_xy")]
    p3_in = [_analytic_P3(m) for m in ("uniaxial_x", "shear_yz", "biaxial_xy")]
    p3_fit = [dict_P[m]["P3"] for m in ("uniaxial_x", "shear_yz", "biaxial_xy")]
    p4_in = [_analytic_P4(m) for m in ("uniaxial_x", "shear_yz", "biaxial_xy")]
    p4_fit = [dict_P[m]["P4"] for m in ("uniaxial_x", "shear_yz", "biaxial_xy")]

    x = np.arange(3)
    width = 0.27
    ax.bar(x - width, p2_in, width, color="tab:gray", edgecolor="black",
           label="P2 analytic")
    ax.bar(x, p2_fit, width, color="tab:blue", edgecolor="black",
           label="P2 fit")
    # overlay P3, P4 as text so the bar panel stays readable while still
    # reporting every order.
    for xi_i, (p3f, p4f) in enumerate(zip(p3_fit, p4_fit)):
        ax.text(x[xi_i], max(p2_in[xi_i], p2_fit[xi_i]) * 1.02,
                "P3=%.2f\nP4=%.2f" % (p3f, p4f),
                ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(["uniaxial_x", "shear_yz", "biaxial_xy"], fontsize=9)
    ax.set_ylabel(r"$P_2$ coefficient (GPa)")
    ax.set_title("(d) P2: analytic vs fit (P3, P4 annotated)")
    legend = ax.legend(loc="upper right", fontsize=9, framealpha=0.9)
    general_modify_legend(legend, linewidth=1.2)
    ax.set_ylim(0, max(max(p2_in), max(p2_fit)) * 1.35)

    fig.savefig(path_image, bbox_inches="tight")
    plt.close(fig)

    # ---- non-blank guard ----
    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError("HOEC image was not created: " + str(path_image))
    image_rgb = plt.imread(path_image)[..., :3]
    if float(np.mean(np.abs(image_rgb - 1.0))) < 0.002:
        raise AssertionError("HOEC image is effectively blank: " + str(path_image))


def run_example(path_output: Path) -> dict[str, object]:
    """Generate synthetic data, fit per-mode polynomials, render the figure.

    Returns a dict carrying per-mode fitted P2/P3/P4, the analytic values,
    and relative errors for the self-check.
    """
    path_output = path_output.resolve()
    path_output.mkdir(parents=True, exist_ok=True)

    modes = build_synthetic_data()
    dict_P: dict[str, dict[str, float]] = {}
    for mode in ("uniaxial_x", "shear_yz", "biaxial_xy"):
        m = modes[mode]
        P0, P1, P2, P3, P4, rms, coeffs = fit_P_demo(m["strain"], m["u"], fitdeg=4)
        dict_P[mode] = {"P2": P2, "P3": P3, "P4": P4, "rms": rms,
                       "coeffs": coeffs}

    path_image = path_output / "post_hoec_energy_demo.png"
    render_figure(modes, dict_P, path_image)

    # ---- summary table ----
    print("=" * 72)
    print("HOEC energy-strain fitting (synthetic cubic Cu-like data)")
    print("  strain grid = " + str(STRAIN_GRID.tolist()))
    print("  noise sigma  = %.1e GPa  (seed=%d)" % (NOISE_SIGMA_GPA, NOISE_SEED))
    print("-" * 72)
    print("%14s %12s %12s %12s %12s" % ("mode", "P2 analyt", "P2 fit", "P3 analyt", "P3 fit"))
    rows: list[dict[str, object]] = []
    for mode in ("uniaxial_x", "shear_yz", "biaxial_xy"):
        p2a = _analytic_P2(mode)
        p2f = dict_P[mode]["P2"]
        p3a = _analytic_P3(mode)
        p3f = dict_P[mode]["P3"]
        p4a = _analytic_P4(mode)
        p4f = dict_P[mode]["P4"]
        print("%14s %12.3f %12.3f %12.3f %12.3f" % (mode, p2a, p2f, p3a, p3f))
        rows.append({
            "mode": mode,
            "P2_analyt": p2a, "P2_fit": p2f,
            "P3_analyt": p3a, "P3_fit": p3f,
            "P4_analyt": p4a, "P4_fit": p4f,
        })
    print("-" * 72)
    print("NOTE: data is synthetic (built from the input constants);")
    print("      fit reproduces fit_P polynomial math, not a real DFT run.")
    print("wrote: " + str(path_image))
    print("=" * 72)

    return {"rows": rows, "dict_P": dict_P, "path_image": path_image}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fit HOEC from synthetic energy-strain data (VASP-free).")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("hoec-energy-output"),
        help="Directory for the HOEC fitting PNG.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    result = run_example(parse_args().output)
    rows = result["rows"]

    # Deterministic assertions: fitted P2 within 5% of analytic inputs.
    for row in rows:
        p2a = float(row["P2_analyt"])
        p2f = float(row["P2_fit"])
        rel = abs(p2f - p2a) / abs(p2a) * 100.0 if p2a != 0.0 else 0.0
        assert rel < 5.0, (
            "%s P2 relative error %.3f%% >= 5%%" % (row["mode"], rel))
    # Sanity: three modes recovered.
    assert len(rows) == 3
    # Image must exist and be non-blank (also checked inside render_figure).
    assert result["path_image"].is_file()
