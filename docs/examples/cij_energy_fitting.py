#!/usr/bin/env python3
"""Fit second-order elastic constants (Cij) from synthetic energy-strain data.

This tutorial is VASP-free and fully deterministic. It demonstrates the
energy-strain method that :mod:`mymetal.post.Cij_energy` implements for real
LAMMPS/vasp deformation runs, but here the strained-cell energies are
synthetic (generated analytically from textbook cubic Cu constants), so no
external calculation is required.

Energy-strain method
--------------------
For small strains the elastic energy density obeys

    U / V0 = (1/2) * C_ij * eta_i * eta_j

(sum over Voigt indices i, j), where ``U`` is the energy relative to the
unstrained reference, ``V0`` the reference volume, and ``eta`` the strain.
For a cubic crystal only three constants are independent: ``C11``, ``C12``,
``C44``. Three deformation modes are sufficient to recover them:

* (a) uniaxial x strain ``eta_xx = eta`` -> ``U/V0 = (1/2) C11 eta^2``
* (b) y-z engineering shear ``eta_yz = eta`` -> ``U/V0 = (1/2) C44 eta^2``
      (the factor of 2 of engineering shear is absorbed in the Voigt sum)
* (c) x-y biaxial strain ``eta_xx = eta_yy = eta``
      -> ``U/V0 = (1/2) (2 C11 + 2 C12) eta^2 = (C11 + C12) eta^2``

Each curve is fit by ``np.polyfit(strain, energy_density, 2)``; the leading
quadratic coefficient is half the relevant elastic combination. This is the
same polynomial-fit math used by :func:`mymetal.post.Cij_energy.fit_cij_energy`,
which is tightly coupled to the LAMMPS directory layout (it reads
``y_post_data.txt`` and ``movie.lammpstrj`` from per-mode subdirectories). To
keep this example VASP-free and self-contained we reproduce that fitting step
directly on synthetic data and document the correspondence inline.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mymetal.universal.plot.general import general_set_all_rcParams


# ---------------------------------------------------------------------------
# Reference constants: textbook ambient Cu (cubic). These are *not* DFT data;
# they are the analytic inputs used to synthesise the strain-energy curves.
# ---------------------------------------------------------------------------
C11_REF_GPA = 168.0   # GPa
C12_REF_GPA = 121.0   # GPa
C44_REF_GPA = 76.0    # GPa
V0_A3 = 11.7          # reference volume per atom, angstrom^3 (FCC Cu ~ a^3/4)

# eV-to-GPa-A^3 conversion: 1 eV/A^3 = 160.21766 GPa. Used only to keep the
# printed units of U/V0 in GPa when we choose to express U in eV; here we
# build energy density directly in GPa, so this factor is documented for
# reference but not required in the fit.
EV_PER_A3_TO_GPA = 160.21766208

# Strain sampling grid (dimensionless) and synthetic-noise parameters.
STRAIN_GRID = np.array([-0.02, -0.01, 0.0, 0.01, 0.02])
NOISE_SEED = 42
NOISE_SIGMA_GPA = 1e-5   # tiny, well below physical Cij precision


def build_synthetic_data() -> dict[str, dict[str, np.ndarray]]:
    """Return synthetic strain-energy-density curves for the three cubic modes.

    Each entry has keys ``strain`` (dimensionless) and ``ed`` (energy density
    in GPa, i.e. U/V0). The analytic forms follow the cubic energy-strain
    relations documented in the module docstring. A reproducible small noise
    is added so the fit is exercised on non-exact data, mirroring real DFT.
    """
    rng = np.random.default_rng(NOISE_SEED)

    eta = STRAIN_GRID.copy()

    # Mode (a): uniaxial x.  U/V0 = 0.5 * C11 * eta^2
    ed_a = 0.5 * C11_REF_GPA * eta ** 2
    ed_a = ed_a + rng.normal(0.0, NOISE_SIGMA_GPA, size=eta.shape)

    # Mode (b): y-z engineering shear.  U/V0 = 0.5 * C44 * eta^2
    ed_b = 0.5 * C44_REF_GPA * eta ** 2
    ed_b = ed_b + rng.normal(0.0, NOISE_SIGMA_GPA, size=eta.shape)

    # Mode (c): x-y biaxial.  U/V0 = (C11 + C12) * eta^2
    ed_c = (C11_REF_GPA + C12_REF_GPA) * eta ** 2
    ed_c = ed_c + rng.normal(0.0, NOISE_SIGMA_GPA, size=eta.shape)

    return {
        "uniaxial_x":  {"strain": eta, "ed": ed_a,
                        "label": r"(a) uniaxial x:  $U/V_0 = \frac{1}{2} C_{11} \eta^2$",
                        "coeff": 0.5 * C11_REF_GPA},
        "shear_yz":    {"strain": eta, "ed": ed_b,
                        "label": r"(b) yz shear:  $U/V_0 = \frac{1}{2} C_{44} \eta^2$",
                        "coeff": 0.5 * C44_REF_GPA},
        "biaxial_xy":  {"strain": eta, "ed": ed_c,
                        "label": r"(c) biaxial xy:  $U/V_0 = (C_{11}+C_{12}) \eta^2$",
                        "coeff": (C11_REF_GPA + C12_REF_GPA)},
    }


def fit_cij_from_modes(modes: dict[str, dict[str, np.ndarray]]) -> dict[str, float]:
    """Fit quadratic energy-density vs strain and extract C11, C12, C44.

    This reproduces the polynomial-fit step of
    :func:`mymetal.post.Cij_energy.fit_cij_energy`: that function calls
    ``np.polyfit(e, ed, 2)`` on each deformation mode and reads the leading
    coefficient ``lp0 = param[0]`` to build the Cij vector. Here we apply the
    same ``np.polyfit`` to each synthetic mode and invert the mode->Cij
    relations:

    * uniaxial x :  lp0 = C11 / 2           ->  C11 = 2 * lp0
    * shear yz   :  lp0 = C44 / 2           ->  C44 = 2 * lp0
    * biaxial xy :  lp0 = (C11 + C12)       ->  C12 = lp0 - C11

    Returns a dict with keys ``C11``, ``C12``, ``C44`` (GPa) plus per-mode fit
    diagnostics.
    """
    out: dict[str, float] = {}

    # uniaxial x  -> C11
    p_a = np.polyfit(modes["uniaxial_x"]["strain"],
                     modes["uniaxial_x"]["ed"], 2)
    c11 = 2.0 * p_a[0]
    out["C11"] = float(c11)
    out["fit_uniaxial_x"] = float(p_a[0])

    # shear yz -> C44
    p_b = np.polyfit(modes["shear_yz"]["strain"],
                     modes["shear_yz"]["ed"], 2)
    c44 = 2.0 * p_b[0]
    out["C44"] = float(c44)
    out["fit_shear_yz"] = float(p_b[0])

    # biaxial xy -> (C11 + C12) ; subtract the already-fitted C11
    p_c = np.polyfit(modes["biaxial_xy"]["strain"],
                     modes["biaxial_xy"]["ed"], 2)
    c12 = p_c[0] - c11
    out["C12"] = float(c12)
    out["fit_biaxial_xy"] = float(p_c[0])

    return out


def render_figure(modes: dict[str, dict[str, np.ndarray]],
                  fitted: dict[str, float],
                  path_image: Path) -> None:
    """Render the 2-panel figure: strain-energy curves and Cij bar comparison."""
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
        axes_titlepad=8,
        lines_markersize=8,
        lines_linewidth=2.0,
    )

    fig, (ax_curves, ax_bar) = plt.subplots(1, 2)

    # ---- Panel A: three strain-energy curves with quadratic fits ----
    colors = {"uniaxial_x": "tab:blue",
              "shear_yz": "tab:green",
              "biaxial_xy": "tab:red"}
    eta_fine = np.linspace(STRAIN_GRID.min(), STRAIN_GRID.max(), 200)
    for key in ("uniaxial_x", "shear_yz", "biaxial_xy"):
        m = modes[key]
        ax_curves.scatter(m["strain"], m["ed"],
                          color=colors[key], marker="o",
                          label=m["label"], zorder=3)
        # refit here only for the drawn curve; matches fit_cij_from_modes math
        p = np.polyfit(m["strain"], m["ed"], 2)
        ax_curves.plot(eta_fine, np.polyval(p, eta_fine),
                       color=colors[key], linestyle="--", linewidth=1.6,
                       zorder=2)
    ax_curves.set_xlabel(r"strain $\eta$")
    ax_curves.set_ylabel(r"energy density $U/V_0$ (GPa)")
    ax_curves.set_title("(A) Synthetic strain-energy curves\n"
                        "(markers: data, dashed: quadratic fit)")
    ax_curves.legend(loc="upper center", fontsize=9, framealpha=0.9)
    ax_curves.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))

    # ---- Panel B: bar chart input vs fitted C11, C12, C44 ----
    labels = ["C11", "C12", "C44"]
    input_vals = [C11_REF_GPA, C12_REF_GPA, C44_REF_GPA]
    fitted_vals = [fitted["C11"], fitted["C12"], fitted["C44"]]
    x = np.arange(len(labels))
    width = 0.36
    ax_bar.bar(x - width / 2, input_vals, width,
               color="tab:gray", edgecolor="black", label="input")
    ax_bar.bar(x + width / 2, fitted_vals, width,
               color="tab:orange", edgecolor="black", label="fitted")
    # value annotations
    for xi, vi, vf in zip(x, input_vals, fitted_vals):
        ax_bar.text(xi - width / 2, vi + 2.0, f"{vi:.1f}",
                    ha="center", va="bottom", fontsize=9)
        ax_bar.text(xi + width / 2, vf + 2.0, f"{vf:.1f}",
                    ha="center", va="bottom", fontsize=9)
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels(labels)
    ax_bar.set_ylabel("elastic constant (GPa)")
    ax_bar.set_title("(B) Input vs fitted cubic Cij\n(synthetic data, seed=42)")
    ax_bar.legend(loc="upper right", fontsize=10, framealpha=0.9)
    ax_bar.set_ylim(0, max(max(input_vals), max(fitted_vals)) * 1.18)

    fig.suptitle("Cij energy-strain fitting (synthetic cubic Cu-like data)",
                 y=0.995, fontsize=13)
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

    path_image = path_output / "cij_energy_fitting.png"
    render_figure(modes, fitted, path_image)

    # ---- summary table ----
    input_map = {"C11": C11_REF_GPA, "C12": C12_REF_GPA, "C44": C44_REF_GPA}
    print("=" * 64)
    print("Cij energy-strain fitting (synthetic cubic Cu-like data)")
    print(f"  reference volume V0 = {V0_A3:.2f} A^3")
    print(f"  strain grid        = {STRAIN_GRID.tolist()}")
    print(f"  noise sigma        = {NOISE_SIGMA_GPA:.1e} GPa  (seed={NOISE_SEED})")
    print("-" * 64)
    print(f"{'constant':>10} {'input(GPa)':>12} {'fitted(GPa)':>13} "
          f"{'rel.err(%)':>12}")
    rows: list[dict[str, object]] = []
    for key in ("C11", "C12", "C44"):
        vi = input_map[key]
        vf = fitted[key]
        rel = abs(vf - vi) / vi * 100.0
        print(f"{key:>10} {vi:>12.3f} {vf:>13.3f} {rel:>12.4f}")
        rows.append({"constant": key, "input": vi, "fitted": vf, "rel_err_pct": rel})
    print("-" * 64)
    print("NOTE: data is synthetic (generated from the input constants),")
    print("      not real DFT. Fit reproduces fit_cij_energy polynomial math.")
    print("wrote: " + str(path_image))
    print("=" * 64)

    return {"rows": rows, "fitted": fitted, "path_image": path_image}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fit cubic Cij from synthetic energy-strain data (VASP-free).")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("cij-energy-output"),
        help="Directory for the Cij fitting PNG.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    result = run_example(parse_args().output)
    rows = result["rows"]

    # Deterministic assertions: fitted constants within 5% of inputs.
    for row in rows:
        assert row["rel_err_pct"] < 5.0, (
            f"{row['constant']} relative error {row['rel_err_pct']:.3f}% >= 5%")
    # Sanity: all three constants recovered.
    assert len(rows) == 3
    # The image must exist and be non-blank (also checked inside render_figure).
    assert result["path_image"].is_file()
