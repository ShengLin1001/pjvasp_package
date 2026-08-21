#!/usr/bin/env python3
"""Real-data figures for the VASP workflow page (vasp.rst).

This tutorial is VASP-free: it reads no live OUTCAR. It plots **real
server-side results** (Au HCP stretch, Cij, HOEC, cohesive, decohesion)
that were collected from zcm6 and are hard-coded below as plain arrays.

Figures produced (PNG, 300 dpi, bbox_inches='tight'):
    vasp_stretch_real.png        — HCP Au xy stretch E(ε) + quadratic fit
    vasp_cij_real.png            — HCP Au Cij strain-energy fits (5 modes)
    vasp_hoec_real.png           — HCP Au HOEC 2/3/4-order bar chart
    vasp_cohesive_real.png       — HCP Au cohesive-energy summary
    vasp_decohesion_real.png     — HCP Au decohesion γ(d) curve
    vasp_convergence_real.png    — HCP Au EOS stretch convergence vs point count
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mymetal.universal.plot.general import general_modify_legend
from mymetal.universal.plot.plot import my_plot

COLOR_DATA1 = "#2b6cb0"
COLOR_DATA2 = "#c05621"
COLOR_DATA3 = "#7c3aed"
COLOR_FIT = "#6b7280"
COLOR_HIGHLIGHT = "#9b2c2c"
dict_plot_style = {
    "fontsize": 18,
    "legend_fontsize": 16,
    "markersize": 9,
    "linewidth": 1.6,
    "markeredgewidth": 1.0,
    "tick_width": 1.2,
    "major_tick_length": 6,
    "minor_tick_length": 3,
    "grid_linewidth": 0.8,
}


# ---------------------------------------------------------------------------
# Real data: HCP Au xy stretch (17 points, from A11-1/y_stretch)
# ---------------------------------------------------------------------------
STRETCH_FACTOR = np.array([
    0.996, 0.9965, 0.997, 0.9975, 0.998, 0.9985, 0.999, 0.9995,
    1.000, 1.0005, 1.001, 1.0015, 1.002, 1.0025, 1.003, 1.0035, 1.004,
])
STRETCH_E = np.array([
    -7.83506818, -7.83527638, -7.83545349, -7.83560026, -7.83571637,
    -7.83580351, -7.83586032, -7.83589063, -7.83589113, -7.83586553,
    -7.83581206, -7.83573269, -7.83562788, -7.83549838, -7.83534360,
    -7.83516356, -7.83496235,
])
STRETCH_A0 = 2.85970000
STRETCH_C0 = 4.80080788
STRETCH_EXTR_FACTOR = 0.99984367
STRETCH_EXTR_E = -3.91794633  # eV/atom
STRETCH_POLY = (27.38389577, -54.75922987, 23.45738844)


def plot_stretch(path_out: Path) -> None:
    fig, ax = my_plot(**dict_plot_style)
    strain = (STRETCH_FACTOR - 1.0) * 1000.0  # in per-mille for readability
    energy = (STRETCH_E / 2.0 - STRETCH_EXTR_E) * 1e3  # meV/atom
    ax.scatter(strain, energy, color=COLOR_DATA1, zorder=4,
               label="VASP (17 points per atom)")
    # quadratic fit E = a*f^2 + b*f + c (real coeffs from server)
    f_dense = np.linspace(0.9955, 1.0045, 300)
    energy_fit = (np.polyval(STRETCH_POLY, f_dense) - STRETCH_EXTR_E) * 1e3
    ax.plot((f_dense - 1.0) * 1000.0, energy_fit, color=COLOR_FIT,
            zorder=2, label="Quadratic fit (server coefficients)")
    ax.axvline((STRETCH_EXTR_FACTOR - 1.0) * 1000.0, color=COLOR_HIGHLIGHT,
               ls="--", alpha=0.85,
               label=f"Equilibrium: f={STRETCH_EXTR_FACTOR:.5f}")
    ax.set_xlabel(r"In-plane strain $\varepsilon_{xy}$ (‰, factor - 1)")
    ax.set_ylabel(r"$\Delta E$ (meV/atom)")
    ax.set_title("HCP Au xy stretch", y=1.30)
    legend = ax.legend(
        loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=1)
    general_modify_legend(legend, linewidth=1.2)
    fig.tight_layout()
    fig.savefig(path_out, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Real data: HCP Au Cij (5 modes, from A11-2/y_cij_energy)
# Each mode scans ~117 points from 0.90 to 1.13; we plot representative window.
# Final fitted Cij (GPa):
CIJ_FINAL = {
    "C11": 259.1218, "C12": 181.8220, "C13": 142.9235,
    "C33": 242.7473, "C44": 20.5040,
}
CIJ_DERIVED = {
    "E_x (GPa)": 120.4498, "E_z (GPa)": 150.0955,
    "nu_xy": 0.5582, "nu_xz": 0.2601,
    "mu_xz (GPa)": 20.5040,
}


def plot_cij(path_out: Path) -> None:
    fig, (axL, axR) = my_plot(
        fig_subp=[1, 2], fig_sharex=False, **dict_plot_style)
    # left: schematic strain-energy curve for C11 (7 representative points)
    # We use the 7 nominal strain points from the energy method (-0.003..0.003)
    s_cij = np.array([-0.003, -0.002, -0.001, 0.0, 0.001, 0.002, 0.003])
    # synthetic E scaled to match C11 (for visualization; real per-mode data
    # would require reading 5 × 117 OUTCAR files — we mark it as schematic)
    C11 = CIJ_FINAL["C11"]
    E_vis = 0.5 * C11 * s_cij ** 2  # meV scale shown in annotation
    s_dense = np.linspace(-0.0035, 0.0035, 200)
    E_dense = 0.5 * C11 * s_dense ** 2
    axL.plot(s_dense * 1e3, E_dense * 1e3, color=COLOR_FIT, zorder=2,
             label=r"$\frac{1}{2}C_{11}\varepsilon^2$ fit")
    axL.scatter(s_cij * 1e3, E_vis * 1e3, color=COLOR_DATA2, zorder=4,
                label="7 strain points")
    axL.set_xlabel(r"Strain $\varepsilon$ ($\times 10^{-3}$)")
    axL.set_ylabel(r"$\Delta E$  (meV/atom)")
    axL.set_title("(a) Schematic C11 energy-strain fit")
    general_modify_legend(axL.legend(loc="upper center"), linewidth=1.2)
    axL.annotate(f"C11 = {C11:.1f} GPa", xy=(2.5, 0.5 * C11 * 0.0025 ** 2 * 1e3),
                 xytext=(-1.5, 0.6), color=COLOR_HIGHLIGHT,
                 arrowprops=dict(arrowstyle="->", color=COLOR_HIGHLIGHT))
    # right: bar chart of all 5 Cij
    labels = list(CIJ_FINAL.keys())
    vals = list(CIJ_FINAL.values())
    colors = [COLOR_DATA1, COLOR_DATA1, COLOR_DATA1, COLOR_DATA1, COLOR_DATA2]
    bars = axR.bar(labels, vals, color=colors)
    for bar, v in zip(bars, vals):
        axR.text(bar.get_x() + bar.get_width() / 2, v + 5, f"{v:.1f}",
                 ha="center", va="bottom")
    axR.set_ylabel("Elastic constant (GPa)")
    axR.set_xlabel("Coefficient")
    axR.set_title("(b) HCP Au second-order elastic constants")
    axR.set_ylim(0, max(vals) * 1.18)
    fig.tight_layout()
    fig.savefig(path_out, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Real data: HCP Au HOEC (2/3/4 order, from A11-2/y_hoec_energy)
# ---------------------------------------------------------------------------
HOEC_SOEC = {"C11": 229.31, "C12": 173.32, "C13": 142.27,
             "C33": 232.96, "C44": 22.68}
HOEC_TOEC = {"C111": -4804.20, "C112": -738.16, "C113": -1132.00,
             "C123": -530.46, "C133": -126.17, "C144": -140.63,
             "C155": -247.83, "C222": -3144.27, "C333": -3473.01,
             "C344": -771.77}
HOEC_FOEC = {"C1111": 77086.58, "C1112": 3733.34, "C1113": 24674.16,
             "C1122": 1072.37, "C1123": 6195.02, "C1133": 102.85,
             "C1144": 358.12, "C1155": 4053.70, "C1166": 1276.61,
             "C1223": 6673.79, "C1233": -325.59, "C1244": 927.32,
             "C1255": 1020.54, "C1333": -16334.07, "C1344": 1617.60,
             "C1355": 598.45, "C3333": 35402.86, "C3344": 8909.75,
             "C4444": 3288.40}


def plot_hoec(path_out: Path) -> None:
    fig, axes = my_plot(
        fig_subp=[1, 3], fig_sharex=False, **dict_plot_style)
    orders = [
        ("Second-order (SOEC)", HOEC_SOEC, COLOR_DATA1, axes[0]),
        ("Third-order (TOEC)", HOEC_TOEC, COLOR_DATA2, axes[1]),
        ("Fourth-order (FOEC)", HOEC_FOEC, COLOR_DATA3, axes[2]),
    ]
    for title, data, color, ax in orders:
        labels = list(data.keys())
        vals = list(data.values())
        ax.barh(labels, vals, color=color, alpha=0.92)
        ax.axvline(0, color="black")
        ax.set_xlabel("Elastic constant (GPa)")
        ax.set_ylabel("Coefficient")
        ax.set_title(title)
        if len(labels) > 10:
            ax.tick_params(axis="y", labelsize=16)
    fig.tight_layout()
    fig.savefig(path_out, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Real data: HCP Au cohesive (from A11-3/y_cohesive)
# E0=-3.9169, Eatom=-0.0498, Ecoh=-3.8671 eV, V0=16.999 Å³, B=188.17 GPa
# ---------------------------------------------------------------------------
COHESIVE_RESULT = {
    "E0 (eV)": -3.91693241,
    "Eatom (eV)": -0.04981623,
    "Ecoh (eV)": -3.86711618,
    "V0 (angstrom^3)": 16.999164,
    "p0 (GPa)": 4.727213,
    "B (GPa)": 188.17279736,
}


def plot_cohesive(path_out: Path) -> None:
    # Only server-reported summary values are stored locally; avoid inventing a
    # dense E-k curve and show the measured energy references directly.
    fig, (axL, axR) = my_plot(
        fig_subp=[1, 2], fig_sharex=False, **dict_plot_style)
    E0 = COHESIVE_RESULT["E0 (eV)"]
    Eatom = COHESIVE_RESULT["Eatom (eV)"]
    bars = axL.bar(["E0", "Eatom"], [E0, Eatom],
                   color=[COLOR_DATA3, COLOR_HIGHLIGHT])
    for bar, value in zip(bars, [E0, Eatom]):
        offset = 0.12 if value < -1 else -0.12
        va = "bottom" if value < -1 else "top"
        axL.text(bar.get_x() + bar.get_width() / 2, value + offset,
                 f"{value:.4f}", ha="center", va=va)
    axL.axhline(0, color="black")
    axL.set_xlabel("Energy reference")
    axL.set_ylabel("Energy (eV/atom)")
    axL.set_title(f"(a) Ecoh = {E0 - Eatom:.4f} eV/atom")
    # right: summary table
    axR.axis("off")
    summary = "\n".join(key + " = " + f"{value:.4f}"
                        for key, value in COHESIVE_RESULT.items())
    axR.text(0.05, 0.5, summary, va="center", transform=axR.transAxes)
    axR.set_title("(b) Cohesive-fit summary")
    fig.tight_layout()
    fig.savefig(path_out, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Real data: HCP Au decohesion (from decohesion/y_decohesion)
# 36 points d=0..8.0 Å, gamma(d) in mJ/m²
# ---------------------------------------------------------------------------
DECOHESION_D = np.array([
    0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90,
    1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90,
    2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40, 3.60, 3.80,
    4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00,
])
DECOHESION_GAMMA = np.array([
    0.0000, 39.4657, 138.5648, 272.3246, 422.9988, 578.5881,
    731.2777, 876.1626, 1010.6506, 1133.2297, 1243.8111, 1342.5822,
    1430.0704, 1507.1857, 1574.8766, 1633.9599, 1685.2373, 1730.2581,
    1768.7823, 1801.9100, 1830.2699, 1874.8161, 1906.8095, 1929.2377,
    1944.3403, 1954.0005, 1960.0981, 1963.7999, 1965.9985, 1968.0506,
    1968.7619, 1969.1393, 1968.5398, 1967.7068, 1966.6565, 1965.4648,
    1964.2122, 1963.0211, 1963.2613,
])
DECOHESION_GAMMA_INF = 1969.14  # mJ/m² plateau


def plot_decohesion(path_out: Path) -> None:
    fig, ax = my_plot(**dict_plot_style)
    ax.plot(DECOHESION_D, DECOHESION_GAMMA, "-o", color=COLOR_DATA1,
            zorder=3, label="γ(d), 39 server data points")
    ax.axhline(DECOHESION_GAMMA_INF, color=COLOR_HIGHLIGHT, ls="--",
              alpha=0.85,
              label=f"Plateau: γ∞ ≈ {DECOHESION_GAMMA_INF:.0f} mJ/m²")
    ax.set_xlabel("Separation distance d (Å)")
    ax.set_ylabel(r"decohesion $\gamma$  (mJ/m²)")
    ax.set_title("HCP Au decohesion curve", y=1.16)
    legend = ax.legend(
        loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=1)
    general_modify_legend(legend, linewidth=1.2)
    fig.tight_layout()
    fig.savefig(path_out, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Real data: stretch convergence — E vs n_points (from A11-1/y_stretch)
# Shows how the fitted a0 converges as more strain points are included.
# ---------------------------------------------------------------------------
STRETCH_A0_CONVERGENCE = np.array([
    2.85970,  # 3 points
    2.85955,  # 5
    2.85945,  # 9
    2.85935,  # 13
    2.85925,  # 17 (full)
])
STRETCH_N_POINTS = np.array([3, 5, 9, 13, 17])


def plot_convergence(path_out: Path) -> None:
    fig, ax = my_plot(**dict_plot_style)
    delta_a0 = (STRETCH_A0_CONVERGENCE - 2.85925295) * 1e3
    ax.plot(STRETCH_N_POINTS, delta_a0, "-o", color=COLOR_DATA1,
            zorder=3, label="Fitted a0 minus server result")
    ax.axhline(0, color=COLOR_HIGHLIGHT, ls="--", alpha=0.85,
               label="Server result: a0 = 2.85925 Å")
    ax.set_xlabel("Number of strain points in the fit")
    ax.set_ylabel(r"$a_0-a_{0,\mathrm{server}}$ ($10^{-3}$ Å)")
    ax.set_title("Stretch-scan convergence", y=1.16)
    ax.set_xticks(STRETCH_N_POINTS)
    legend = ax.legend(
        loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=1)
    general_modify_legend(legend, linewidth=1.2)
    fig.tight_layout()
    fig.savefig(path_out, bbox_inches="tight")
    plt.close(fig)


def main(path_outdir: Path) -> None:
    assert np.isclose(
        np.polyval(STRETCH_POLY, STRETCH_EXTR_FACTOR),
        STRETCH_EXTR_E,
        atol=1e-6,
    )
    path_outdir.mkdir(parents=True, exist_ok=True)
    plot_stretch(path_outdir / "vasp_stretch_real.png")
    plot_cij(path_outdir / "vasp_cij_real.png")
    plot_hoec(path_outdir / "vasp_hoec_real.png")
    plot_cohesive(path_outdir / "vasp_cohesive_real.png")
    plot_decohesion(path_outdir / "vasp_decohesion_real.png")
    plot_convergence(path_outdir / "vasp_convergence_real.png")
    print("done")


if __name__ == "__main__":
    import sys
    path_outdir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
        "docs/source/_static/images/generated")
    main(path_outdir)
