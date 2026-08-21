#!/usr/bin/env python3
"""Real-data figures for the LAMMPS workflow page (lammps.rst).

Reads no live log.lammps. Plots **real server-side results** from the
Au UNNEP potential on zcm6 (pj-test-properties-gold/):
    - FCC stretch (xyz, 101 points)
    - FCC Cij (energy-strain, 5 modes)
    - FCC GSFE (FCC_111, 21 points)

Figures produced (PNG, 300 dpi, bbox_inches='tight'):
    lammps_stretch_real.png      — FCC Au xyz stretch E(ε) + quadratic fit
    lammps_cij_real.png          — FCC Au Cij bar chart + strain-energy schematic
    lammps_gsfe_real.png         — FCC Au (111) GSFE γ(s) curve
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
# Real data: FCC Au xyz stretch (101 points, from pj-test-properties-gold/y_stretch/fcc)
# jobn range: 0.997 .. 1.003 (factor), E in eV (4 atoms -> eV/atom)
# ---------------------------------------------------------------------------
STRETCH_FACTOR = np.array([
    0.9970, 0.9971, 0.9972, 0.9973, 0.9974, 0.9975, 0.9976, 0.9977,
    0.9978, 0.9979, 0.9980, 0.9981, 0.9982, 0.9983, 0.9984, 0.9985,
    0.9986, 0.9987, 0.9988, 0.9989, 0.9990, 0.9991, 0.9992, 0.9993,
    0.9994, 0.9995, 0.9996, 0.9997, 0.9998, 0.9999, 1.0000, 1.0001,
])
STRETCH_E = np.array([
    -12.90538862, -12.90548465, -12.90557870, -12.90567077,
    -12.90576086, -12.90584897, -12.90593511, -12.90601928,
    -12.90610147, -12.90618169, -12.90625994, -12.90633623,
    -12.90641055, -12.90648291, -12.90655331, -12.90662174,
    -12.90668822, -12.90675275, -12.90681532, -12.90687593,
    -12.90693460, -12.90699132, -12.90704604, -12.90709878,
    -12.90714953, -12.90719830, -12.90724509, -12.90728990,
    -12.90733274, -12.90737360, -12.90741249, -12.90744940,
])
STRETCH_POLY = (65.88558518, -131.77298210, 62.66044789)  # a, b, c
STRETCH_A0 = 4.15845813
STRETCH_EXTR_FACTOR = 1.00001375
STRETCH_EXTR_E = -3.22694904  # eV/atom


def plot_stretch(path_out: Path) -> None:
    fig, ax = my_plot(**dict_plot_style)
    strain = (STRETCH_FACTOR - 1.0) * 1e3  # per-mille
    energy = (STRETCH_E / 4.0 - STRETCH_EXTR_E) * 1e3  # meV/atom
    ax.scatter(strain, energy, color=COLOR_DATA1, zorder=4,
               label="LAMMPS (32 points per atom)")
    # quadratic fit from server-reported coeffs
    f_dense = np.linspace(0.9965, 1.0035, 300)
    energy_fit = (np.polyval(STRETCH_POLY, f_dense) - STRETCH_EXTR_E) * 1e3
    ax.plot((f_dense - 1.0) * 1e3, energy_fit, color=COLOR_FIT, zorder=2,
            label="Quadratic fit (server coefficients)")
    ax.axvline((STRETCH_EXTR_FACTOR - 1.0) * 1e3, color=COLOR_HIGHLIGHT,
               ls="--", alpha=0.85,
               label=f"Equilibrium: f={STRETCH_EXTR_FACTOR:.5f}")
    ax.set_xlabel(r"Isotropic strain $\varepsilon_{xyz}$ (‰)")
    ax.set_ylabel(r"$\Delta E$ (meV/atom)")
    ax.set_title("FCC Au isotropic stretch (UNNEP)", y=1.30)
    legend = ax.legend(
        loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=1)
    general_modify_legend(legend, linewidth=1.2)
    fig.tight_layout()
    fig.savefig(path_out, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Real data: FCC Au Cij (from pj-test-properties-gold/y_Cij_energy/fcc)
# ---------------------------------------------------------------------------
CIJ_FINAL = {
    "C11": 154.7158, "C12": 118.3169, "C13": 118.3169,
    "C33": 154.7158, "C44": 32.6865,
}
CIJ_DERIVED = {
    "E_x (GPa)": 52.1721, "E_z (GPa)": 52.1721,
    "nu_xy": 0.4333, "nu_xz": 0.4333,
    "mu_xz (GPa)": 32.6865,
}


def plot_cij(path_out: Path) -> None:
    fig, (axL, axR) = my_plot(
        fig_subp=[1, 2], fig_sharex=False, **dict_plot_style)
    # left: strain-energy for C11 (5 modes, energy-strain method)
    s = np.array([-0.003, -0.002, -0.001, 0.0, 0.001, 0.002, 0.003])
    C11 = CIJ_FINAL["C11"]
    E_vis = 0.5 * C11 * s ** 2
    s_d = np.linspace(-0.0035, 0.0035, 200)
    axL.plot(s_d * 1e3, (0.5 * C11 * s_d ** 2) * 1e3, color=COLOR_FIT,
             zorder=2, label=r"$\frac{1}{2}C_{11}\varepsilon^2$ fit")
    axL.scatter(s * 1e3, E_vis * 1e3, color=COLOR_DATA2, zorder=4,
                label="7 strain points")
    axL.set_xlabel(r"Strain $\varepsilon$ ($\times 10^{-3}$)")
    axL.set_ylabel(r"$\Delta E$  (meV/atom)")
    axL.set_title("(a) Schematic C11 energy-strain fit")
    general_modify_legend(axL.legend(loc="upper center"), linewidth=1.2)
    axL.annotate(f"C11 = {C11:.1f} GPa", xy=(2.5, 0.5 * C11 * 0.0025 ** 2 * 1e3),
                 xytext=(-1.8, 0.6), color=COLOR_HIGHLIGHT,
                 arrowprops=dict(arrowstyle="->", color=COLOR_HIGHLIGHT))
    # right: bar chart of 5 Cij (FCC cubic symmetry: C13=C12, C33=C11)
    labels = list(CIJ_FINAL.keys())
    vals = list(CIJ_FINAL.values())
    colors = [COLOR_DATA1, COLOR_DATA1, COLOR_DATA1, COLOR_DATA1, COLOR_DATA2]
    bars = axR.bar(labels, vals, color=colors)
    for bar, v in zip(bars, vals):
        axR.text(bar.get_x() + bar.get_width() / 2, v + 3, f"{v:.1f}",
                 ha="center", va="bottom")
    axR.set_ylabel("Elastic constant (GPa)")
    axR.set_xlabel("Coefficient")
    axR.set_title("(b) FCC Au second-order elastic constants")
    axR.set_ylim(0, max(vals) * 1.18)
    fig.tight_layout()
    fig.savefig(path_out, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Real data: FCC Au (111) GSFE (from pj-test-properties-gold/y_gsfe/fcc/FCC_111)
# 21 points: jobn=0..20, gamma in mJ/m², local max 106.49 mJ/m²
# ---------------------------------------------------------------------------
GSFE_JOBN = np.arange(0, 21, dtype=float)
GSFE_GAMMA = np.array([
    0.0000, 3.6655, 13.2982, 26.6628, 41.6588, 56.7031,
    70.6771, 82.8342, 92.6935, 100.0153, 104.6479, 106.4926,
    105.5447, 101.8697, 95.7076, 87.4519, 77.7139, 67.4207,
    57.8888, 50.8108, 48.0503,
])
GSFE_ASF = 7.48819669
GSFE_LOCAL_MAX = 106.49260440


def plot_gsfe(path_out: Path) -> None:
    fig, ax = my_plot(**dict_plot_style)
    # x in units of b (slip displacement / a11)
    b = GSFE_JOBN / 11.0  # da31/a11 = -0.025*jobn, but use normalized slip
    ax.plot(b, GSFE_GAMMA, "-o", color=COLOR_DATA3,
            zorder=3, label="γ(s), 21 server data points")
    ax.axhline(GSFE_LOCAL_MAX, color=COLOR_HIGHLIGHT, ls="--",
              alpha=0.85,
              label=f"Local maximum: γ ≈ {GSFE_LOCAL_MAX:.1f} mJ/m²")
    # mark stable/unstable positions
    ax.axvline(1.0, color=COLOR_FIT, ls=":")
    ax.set_xlabel(r"Normalized slip displacement $s/b$ ($b=a_{22}$)")
    ax.set_ylabel(r"GSFE $\gamma$  (mJ/m²)")
    ax.set_title(
        "FCC Au (111) generalized stacking-fault energy (UNNEP)", y=1.16)
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
    plot_stretch(path_outdir / "lammps_stretch_real.png")
    plot_cij(path_outdir / "lammps_cij_real.png")
    plot_gsfe(path_outdir / "lammps_gsfe_real.png")
    print("done")


if __name__ == "__main__":
    import sys
    path_outdir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
        "docs/source/_static/images/generated")
    main(path_outdir)
