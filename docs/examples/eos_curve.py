#!/usr/bin/env python3
"""Fit and visualize Murnaghan / Birch-Murnaghan EOS from synthetic data.

This tutorial is VASP-free. It demonstrates the workflow that would follow a
real EOS calculation: take a series of (volume, energy) points, fit them with
an equation of state, and report the equilibrium volume, bulk modulus and
pressure derivative.

The energy data is synthesized from a Murnaghan EOS with Cu-like parameters
so the fit can be checked against the input. The script then fits the same
data with the ASE :class:`ase.eos.EquationOfState` engine using both the
``murnaghan`` and ``birchmurnaghan`` forms, and additionally runs a
:func:`scipy.optimize.curve_fit` on the Murnaghan form to expose the
pressure derivative ``B0'``. A two-panel figure shows the fits and their
residuals.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ase.eos import EquationOfState
from scipy.optimize import curve_fit

from mymetal.universal.plot.general import (
    general_modify_legend,
    general_set_all_rcParams,
)


# Reference "true" parameters used to synthesize the data. The values are
# chosen close to a real Cu Fm-3m EOS so the numbers look familiar; they are
# NOT a published benchmark and the tutorial does not claim so.
TRUE_E0_EV = -8.2
TRUE_V0_A3 = 11.7          # ~ Cu conventional-cell volume per atom
TRUE_B0_GPA = 140.0
TRUE_B0P = 5.2

# 1 eV/Å^3 = 160.21766 GPa
EV_PER_A3_TO_GPA = 160.21766


def murnaghan_energy(volume, e0, v0, b0_gpa, b0p):
    """Murnaghan EOS energy (PRB 28, 5480, 1983).

    E(V) = E0 + B0*V/B0' * [(V0/V)^B0'/(B0'-1) + 1] - B0*V0/(B0'-1)

    Args:
        volume (array_like): Volume per atom in Å^3.
        e0 (float): Equilibrium energy in eV.
        v0 (float): Equilibrium volume in Å^3.
        b0_gpa (float): Bulk modulus in GPa.
        b0p (float): Pressure derivative B0'.

    Returns:
        numpy.ndarray: Energy in eV.
    """
    volume = np.asarray(volume, dtype=float)
    b0 = b0_gpa / EV_PER_A3_TO_GPA
    eta = (v0 / volume) ** b0p
    return (
        e0
        + b0 * volume / b0p * (eta / (b0p - 1.0) + 1.0)
        - b0 * v0 / (b0p - 1.0)
    )


def birch_murnaghan_energy(volume, e0, v0, b0_ev_a3, b0p):
    """Third-order Birch-Murnaghan EOS energy (PRB 70, 224107)."""
    volume = np.asarray(volume, dtype=float)
    x = volume / v0
    xi = x ** (2.0 / 3.0)
    return e0 + 9.0 * b0_ev_a3 * v0 / 16.0 * (
        (xi - 1.0) ** 3 * b0p
        + (xi - 1.0) ** 2 * (6.0 - 4.0 * xi)
    )


def make_synthetic_ev(
    e0: float = TRUE_E0_EV,
    v0: float = TRUE_V0_A3,
    b0_gpa: float = TRUE_B0_GPA,
    b0p: float = TRUE_B0P,
    n_points: int = 13,
    strain_amp: float = 0.06,
    noise_std_ev: float = 0.0015,
    seed: int = 20260729,
) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    """Return a deterministic (volume, energy) table generated from a Murnaghan EOS.

    A small Gaussian noise is added to mimic the residual scatter of real DFT
    data. The seed is fixed so the fit assertions are reproducible.

    Args:
        e0 (float): Equilibrium energy in eV.
        v0 (float): Equilibrium volume in Å^3.
        b0_gpa (float): Bulk modulus in GPa.
        b0p (float): Pressure derivative B0'.
        n_points (int): Number of points in the scan.
        strain_amp (float): Half-range of the volume strain (V/V0 - 1).
        noise_std_ev (float): Standard deviation of the Gaussian noise in eV.
        seed (int): RNG seed.

    Returns:
        tuple: (volume_A3, energy_eV, meta) where meta stores the input
        parameters in the same units the fit reports.
    """
    rng = np.random.default_rng(seed)
    strains = np.linspace(-strain_amp, strain_amp, n_points)
    volumes = v0 * (1.0 + strains)
    energies = murnaghan_energy(volumes, e0, v0, b0_gpa, b0p)
    energies = energies + rng.normal(0.0, noise_std_ev, size=energies.shape)
    meta = {
        "e0_ev": e0,
        "v0_A3": v0,
        "b0_gpa": b0_gpa,
        "b0p": b0p,
    }
    return volumes, energies, meta


def fit_with_ase(
    volumes: np.ndarray,
    energies: np.ndarray,
    eos: str = "murnaghan",
) -> dict[str, float]:
    """Fit an EOS with :class:`ase.eos.EquationOfState` and return the params.

    Args:
        volumes (numpy.ndarray): Volume per atom in Å^3.
        energies (numpy.ndarray): Energy in eV.
        eos (str): EOS name passed to ASE (``murnaghan``, ``birchmurnaghan``,
            ``sjeos``, ``vinet`` ...).

    Returns:
        dict: ``e0_ev``, ``v0_A3``, ``b0_gpa`` and the ASE ``eos`` name.
    """
    fitter = EquationOfState(volumes.tolist(), energies.tolist(), eos=eos)
    v0, e0, b0_ev_a3 = fitter.fit()
    return {
        "eos": eos,
        "e0_ev": float(e0),
        "v0_A3": float(v0),
        "b0_gpa": float(b0_ev_a3 * EV_PER_A3_TO_GPA),
    }


def fit_murnaghan_with_b0p(
    volumes: np.ndarray,
    energies: np.ndarray,
) -> dict[str, float]:
    """Fit the Murnaghan EOS with scipy to also recover ``B0'``.

    ASE's :class:`EquationOfState` does not always expose the pressure
    derivative, so this helper runs :func:`scipy.optimize.curve_fit` on the
    same Murnaghan form with a robust initial guess built from the data.

    Args:
        volumes (numpy.ndarray): Volume per atom in Å^3.
        energies (numpy.ndarray): Energy in eV.

    Returns:
        dict: ``e0_ev``, ``v0_A3``, ``b0_gpa``, ``b0p``.
    """
    idx_min = int(np.argmin(energies))
    p0 = [
        float(energies[idx_min]),
        float(volumes[idx_min]),
        140.0,
        4.0,
    ]
    popt, _ = curve_fit(murnaghan_energy, volumes, energies, p0=p0, maxfev=20000)
    return {
        "e0_ev": float(popt[0]),
        "v0_A3": float(popt[1]),
        "b0_gpa": float(popt[2]),
        "b0p": float(popt[3]),
    }


def render_eos_figure(
    volumes: np.ndarray,
    energies: np.ndarray,
    fit_murn: dict[str, float],
    fit_bm: dict[str, float],
    fit_murn_b0p: dict[str, float],
    meta: dict[str, float],
    path_image: Path,
) -> None:
    """Render a two-panel EOS figure: raw data + fits, and residuals."""
    path_image.parent.mkdir(parents=True, exist_ok=True)

    general_set_all_rcParams(
        backend="Agg",
        axes_grid=True,
        figure_subp=(1, 2),
        figure_one_figsize=(6.4, 4.6),
        figure_left=0.75,
        figure_top=0.55,
        figure_wspace=1.0,
        figure_hspace=0.6,
        font_family=["DejaVu Sans"],
        fontsize=11,
        axes_titlepad=6,
    )
    fig, axes = plt.subplots(1, 2)

    v_dense = np.linspace(volumes.min(), volumes.max(), 400)
    e_murn = murnaghan_energy(
        v_dense,
        fit_murn_b0p["e0_ev"],
        fit_murn_b0p["v0_A3"],
        fit_murn_b0p["b0_gpa"],
        fit_murn_b0p["b0p"],
    )
    e_bm = birch_murnaghan_energy(
        v_dense,
        fit_bm["e0_ev"],
        fit_bm["v0_A3"],
        fit_bm["b0_gpa"] / EV_PER_A3_TO_GPA,
        fit_murn_b0p["b0p"],
    )

    ax = axes[0]
    ax.plot(volumes, energies, "o", color="black", label="synthetic data")
    ax.plot(v_dense, e_murn, "-", color="#1f77b4", label="Murnaghan fit")
    ax.plot(v_dense, e_bm, "--", color="#d62728", label="Birch-Murnaghan fit")
    ax.axvline(meta["v0_A3"], color="gray", linestyle=":", label="true $V_0$")
    ax.set_xlabel(r"Volume per atom (Å$^3$)")
    ax.set_ylabel(r"Energy per atom (eV)")
    ax.set_title("Equation of state", y=1.28)
    legend = ax.legend(
        loc="lower center", bbox_to_anchor=(0.5, 1.01), fontsize=9, ncol=2)
    general_modify_legend(legend, linewidth=1.0)

    e_murn_res = murnaghan_energy(
        volumes,
        fit_murn_b0p["e0_ev"],
        fit_murn_b0p["v0_A3"],
        fit_murn_b0p["b0_gpa"],
        fit_murn_b0p["b0p"],
    )
    e_bm_res = birch_murnaghan_energy(
        volumes,
        fit_bm["e0_ev"],
        fit_bm["v0_A3"],
        fit_bm["b0_gpa"] / EV_PER_A3_TO_GPA,
        fit_murn_b0p["b0p"],
    )
    ax = axes[1]
    ax.plot(volumes, energies - e_murn_res, "o", color="#1f77b4", label="Murnaghan")
    ax.plot(volumes, energies - e_bm_res, "s", color="#d62728", label="Birch-Murnaghan")
    ax.axhline(0.0, color="gray", linestyle=":")
    ax.set_xlabel(r"Volume per atom (Å$^3$)")
    ax.set_ylabel(r"Residual (eV)")
    ax.set_title("Fit residuals", y=1.18)
    legend = ax.legend(
        loc="lower center", bbox_to_anchor=(0.5, 1.01), fontsize=9, ncol=2)
    general_modify_legend(legend, linewidth=1.0)

    fig.tight_layout()
    fig.savefig(path_image, bbox_inches="tight")
    plt.close(fig)

    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError("EOS image was not created: " + str(path_image))


def run_example(path_output: Path) -> tuple[dict[str, float], dict[str, float], dict[str, float], dict[str, float]]:
    """Generate the data, fit both EOS, render the figure, and return the params."""
    path_output = path_output.resolve()
    path_output.mkdir(parents=True, exist_ok=True)

    volumes, energies, meta = make_synthetic_ev()
    fit_murn = fit_with_ase(volumes, energies, eos="murnaghan")
    fit_bm = fit_with_ase(volumes, energies, eos="birchmurnaghan")
    fit_murn_b0p = fit_murnaghan_with_b0p(volumes, energies)

    path_image = path_output / "eos_curve.png"
    render_eos_figure(volumes, energies, fit_murn, fit_bm, fit_murn_b0p, meta, path_image)

    print("true:    V0=%.4f Å^3  B0=%.2f GPa  B0'=%.3f  E0=%.6f eV"
          % (meta["v0_A3"], meta["b0_gpa"], meta["b0p"], meta["e0_ev"]))
    print("murn(ASE):    V0=%.4f Å^3  B0=%.2f GPa            E0=%.6f eV"
          % (fit_murn["v0_A3"], fit_murn["b0_gpa"], fit_murn["e0_ev"]))
    print("bm  (ASE):    V0=%.4f Å^3  B0=%.2f GPa            E0=%.6f eV"
          % (fit_bm["v0_A3"], fit_bm["b0_gpa"], fit_bm["e0_ev"]))
    print("murn(scipy):  V0=%.4f Å^3  B0=%.2f GPa  B0'=%.3f  E0=%.6f eV"
          % (fit_murn_b0p["v0_A3"], fit_murn_b0p["b0_gpa"], fit_murn_b0p["b0p"], fit_murn_b0p["e0_ev"]))
    print("wrote: " + str(path_image))
    return fit_murn, fit_bm, fit_murn_b0p, meta


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fit Murnaghan and Birch-Murnaghan EOS to synthetic data.")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("eos-output"),
        help="Directory for the EOS PNG.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    fit_murn, fit_bm, fit_murn_b0p, meta = run_example(parse_args().output)
    # Both ASE fits must recover the input volume within a small tolerance.
    assert abs(fit_murn["v0_A3"] - meta["v0_A3"]) < 0.05, fit_murn
    assert abs(fit_bm["v0_A3"] - meta["v0_A3"]) < 0.05, fit_bm
    # Bulk modulus is more sensitive to noise; allow a wider band.
    assert abs(fit_murn["b0_gpa"] - meta["b0_gpa"]) < 25.0, fit_murn
    assert abs(fit_bm["b0_gpa"] - meta["b0_gpa"]) < 25.0, fit_bm
    # The scipy Murnaghan fit must additionally recover B0' coarsely.
    assert abs(fit_murn_b0p["b0p"] - meta["b0p"]) < 3.0, fit_murn_b0p
