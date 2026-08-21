#!/usr/bin/env python3
"""Demonstrate ionic-relaxation convergence (energy and force vs step).

This tutorial is VASP-free. It demonstrates :mod:`mymetal.post.relax_convergence`,
which reads the per-ionic-step convergence data files produced by the bash
helper ``pei_vasp_univ_extract_convergence`` and renders one plot per
``y_dir`` sub-directory showing the ionic-relaxation trajectory:

    * ``|deltaE|`` relative to the last frame, in meV/atom
    * max force norm, in eV/Ang, together with the |EDIFFG| force criterion

The real workflow parses ``y_post_convergence_<name>.txt`` via
:func:`mymetal.post.relax_convergence.read_univ_post_convergence` and plots
through :func:`mymetal.universal.plot.workflow.my_plot_relax_convergence`; this
demo synthesises the ionic-step trajectory analytically (energy and force both
decay geometrically toward convergence, plus a reproducible small noise) and
draws the same two-panel log-y figure with the EDIFFG criterion line. No
external calculation is required.
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
# Reference "true" relaxation parameters. These are *not* DFT data; they are
# the analytic inputs used to synthesise the ionic-step trajectory.
# ---------------------------------------------------------------------------
N_STEPS = 20
NATOMS = 4
EDIFFG = -0.02          # eV/Ang force criterion (negative = force mode)
ISIF = 4
E_FINAL_EV = -14.96     # converged total energy, eV (~4 atoms * -3.74)
DE_PER_STEP = 0.55       # meV/atom decay factor per step (geometric)
F_INIT = 1.8             # initial max force, eV/Ang
F_DECAY = 0.72           # force decay factor per step (geometric)

NOISE_SEED = 42
NOISE_SIGMA_E = 1e-3      # eV
NOISE_SIGMA_F = 2e-3      # eV/Ang


def build_synthetic_data() -> dict[str, np.ndarray]:
    """Return a synthetic ionic-step relaxation trajectory.

    Energy approaches the converged value as a geometric decay (each step
    closes a fixed fraction of the remaining gap), and the max force norm
    decays geometrically toward zero. A reproducible small noise is added so
    the log-y plot is exercised on non-exact data, mirroring real ionic
    relaxation. The returned dict carries ``lframe`` (1-indexed step indices),
    ``lenergy`` (total energy in eV), and ``lforce`` (max force in eV/Ang).
    """
    rng = np.random.default_rng(NOISE_SEED)
    frames = np.arange(1, N_STEPS + 1, dtype=int)
    # energy: start above E_FINAL, decay geometrically
    e_gap = 0.040   # initial gap 40 meV/atom
    lenergy = np.zeros(N_STEPS)
    for i in range(N_STEPS):
        lenergy[i] = E_FINAL_EV + NATOMS * e_gap * (DE_PER_STEP ** i)
    lenergy = lenergy + rng.normal(0.0, NOISE_SIGMA_E, size=N_STEPS)
    # force: decay geometrically from F_INIT
    lforce = F_INIT * (F_DECAY ** (frames - 1))
    lforce = lforce + rng.normal(0.0, NOISE_SIGMA_F, size=N_STEPS)
    # force must stay positive for a clean log axis
    lforce = np.abs(lforce)
    return {"lframe": frames, "lenergy": lenergy, "lforce": lforce,
            "natoms": NATOMS, "ediffg": EDIFFG, "isif": ISIF,
            "e_final": E_FINAL_EV}


def render_figure(data: dict[str, np.ndarray],
                  path_image: Path) -> None:
    """Render a 2-panel log-y figure: |deltaE|/atom and max force vs step."""
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

    fig, (ax_e, ax_f) = plt.subplots(1, 2)

    lframe = data["lframe"]
    natoms = data["natoms"]
    e_final = data["e_final"]
    ediffg = data["ediffg"]
    isif = data["isif"]
    lenergy_rel = np.abs(data["lenergy"] - e_final) * 1000.0 / natoms
    lforce_abs = np.abs(data["lforce"])

    # ---- Panel A: |deltaE|/atom vs step, log y ----
    ax_e.plot(lframe, lenergy_rel, marker="o", linestyle="-",
              color="tab:blue", label=r"$|\Delta E|$ per atom")
    ax_e.set_yscale("log")
    ax_e.set_xlabel("ionic step")
    ax_e.set_ylabel(r"$|\Delta E|$ (meV/atom)")
    ax_e.set_title("(A) Energy convergence")
    legend = ax_e.legend(loc="upper right", fontsize=10, framealpha=0.9)
    general_modify_legend(legend, linewidth=1.2)

    # ---- Panel B: max force vs step, log y, with |EDIFFG| line ----
    ax_f.plot(lframe, lforce_abs, marker="o", linestyle="-",
              color="tab:blue", label="max force")
    ax_f.set_yscale("log")
    ax_f.axhline(y=-ediffg, color="tab:orange", linestyle="--",
                 linewidth=1.6, label="|EDIFFG| = %g eV/Ang" % float(-ediffg))
    ax_f.set_xlabel("ionic step")
    ax_f.set_ylabel("Max Force (eV/Ang)")
    ax_f.set_title("(B) Force convergence")
    legend = ax_f.legend(loc="upper right", fontsize=10, framealpha=0.9)
    general_modify_legend(legend, linewidth=1.2)
    # ISIF annotation, mirroring my_plot_relax_convergence
    ax_f.text(0.02, 0.98, "ISIF=%d" % int(isif),
              transform=ax_f.transAxes, ha="left", va="top",
              fontsize=10, color="black")

    fig.savefig(path_image, bbox_inches="tight")
    plt.close(fig)

    # ---- non-blank guard ----
    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError("relax-convergence image was not created: "
                             + str(path_image))
    image_rgb = plt.imread(path_image)[..., :3]
    if float(np.mean(np.abs(image_rgb - 1.0))) < 0.002:
        raise AssertionError("relax-convergence image is effectively blank: "
                              + str(path_image))


def run_example(path_output: Path) -> dict[str, object]:
    """Generate synthetic relaxation data, render the figure, return summary."""
    path_output = path_output.resolve()
    path_output.mkdir(parents=True, exist_ok=True)

    data = build_synthetic_data()

    path_image = path_output / "post_relax_convergence_demo.png"
    render_figure(data, path_image)

    # ---- summary ----
    lenergy_rel = np.abs(data["lenergy"] - data["e_final"]) * 1000.0 / NATOMS
    lforce_abs = np.abs(data["lforce"])
    print("=" * 64)
    print("Ionic-relaxation convergence demo (synthetic data)")
    print("  n_steps = %d   natoms = %d" % (N_STEPS, NATOMS))
    print("  EDIFFG  = %.4f eV/Ang   ISIF = %d" % (EDIFFG, ISIF))
    print("  noise sigma = %.1e eV / %.1e eV/Ang  (seed=%d)"
          % (NOISE_SIGMA_E, NOISE_SIGMA_F, NOISE_SEED))
    print("-" * 64)
    print("  first |deltaE|/atom = %.4f meV/atom" % float(lenergy_rel[0]))
    print("  last  |deltaE|/atom = %.4f meV/atom" % float(lenergy_rel[-1]))
    print("  first max force     = %.4f eV/Ang" % float(lforce_abs[0]))
    print("  last  max force     = %.4f eV/Ang" % float(lforce_abs[-1]))
    print("-" * 64)
    print("NOTE: data is synthetic (geometric decay + noise);")
    print("      plot mirrors my_plot_relax_convergence, not real DFT.")
    print("wrote: " + str(path_image))
    print("=" * 64)

    return {"data": data, "path_image": path_image}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot synthetic ionic-relaxation convergence (VASP-free).")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("post-relax-convergence-output"),
        help="Directory for the relaxation-convergence PNG.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    result = run_example(parse_args().output)
    data = result["data"]

    # Deterministic assertions: force decays below EDIFFG by the last step.
    lforce_abs = np.abs(data["lforce"])
    assert float(lforce_abs[-1]) < float(-EDIFFG), (
        "last force %.4f eV/Ang >= |EDIFFG| %.4f"
        % (float(lforce_abs[-1]), float(-EDIFFG)))
    # Energy delta must shrink from first to last (convergence).
    lenergy_rel = np.abs(data["lenergy"] - data["e_final"]) * 1000.0 / NATOMS
    assert float(lenergy_rel[-1]) < float(lenergy_rel[0]), (
        "energy did not converge: last %.4f >= first %.4f"
        % (float(lenergy_rel[-1]), float(lenergy_rel[0])))
    # Image must exist and be non-blank (also checked inside render_figure).
    assert result["path_image"].is_file()
