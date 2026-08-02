#!/usr/bin/env python3
"""Extended-XYZ trajectory + delimited-file I/O round trips.

This tutorial is VASP-free and deterministic. It demonstrates two unrelated
``mymetal`` I/O workflows so that the docs page can contrast a binary-free
trajectory format against a plain-text tabular format:

* **Part A — Extended XYZ trajectory I/O** (``mymetal.io.extxyz``).
  A small 5-frame MD-like trajectory is built from an FCC Cu bulk cell by
  applying a deterministic per-frame perturbation (seed = 42). The trajectory
  is written with ``ase.io.write(..., format='extxyz')`` and read back with
  ``extxyz_to_atomlist``. The number of frames, the per-frame chemical
  formula, and the trajectory of the first atom are printed and asserted.

* **Part B — Delimited file I/O** (``mymetal.io.general``).
  A synthetic convergence table (``encut`` from 300 to 550 eV with a
  monotonically converging ``energy`` and a small ``pressure`` column) is
  written with ``general_write`` and read back with ``general_read``. The
  round-tripped numeric values are checked against the original with
  ``numpy.allclose``.

A 2-panel figure is rendered:

* Left panel — first-atom trajectory across the 5 frames, overlaid with the
  frame-0 and frame-4 positions of all atoms (scatter, two colors).
* Right panel — ``encut`` vs ``energy`` line plot with markers, drawn from the
  round-tripped DataFrame, so the figure visibly exercises both I/O paths.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # MUST precede import pyplot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ase import Atoms
from ase.build import bulk
from ase.io import write as ase_write

from mymetal.io.extxyz import extxyz_to_atomlist
from mymetal.io.general import general_read, general_write
from mymetal.universal.plot.general import general_set_all_rcParams


# Textbook ambient lattice constant for FCC Cu (Å). Hardcoded so the script is
# fully deterministic and free of network access.
A_FCC_CU = 3.61
# Number of frames in the synthetic trajectory.
N_FRAMES = 5
# Per-frame position perturbation amplitude (Å). Small enough to stay in the
# cell, large enough to be visible in the figure.
PERTURB_AMP = 0.05
# Reproducibility seed for the perturbation generator.
SEED = 42


# ---------------------------------------------------------------------------
# Part A: Extended XYZ trajectory I/O
# ---------------------------------------------------------------------------


def build_trajectory() -> list[Atoms]:
    """Build a 5-frame MD-like trajectory of FCC Cu (deterministic, seed=42).

    The base cell is a 2x2x2 repeat of the cubic FCC Cu conventional cell so
    that the trajectory has enough atoms to look interesting in the figure.
    Each frame is a copy of the base with positions perturbed by a fixed
    random vector drawn from a seeded generator. Frames share the same cell
    and chemical identity; only positions move.
    """
    base = bulk("Cu", "fcc", a=A_FCC_CU, cubic=True).repeat((2, 2, 2))
    rng = np.random.default_rng(SEED)
    frames: list[Atoms] = []
    for _ in range(N_FRAMES):
        atoms = base.copy()
        delta = rng.normal(loc=0.0, scale=PERTURB_AMP, size=atoms.positions.shape)
        atoms.positions = atoms.positions + delta
        atoms.wrap()
        frames.append(atoms)
    return frames


def write_trajectory(frames: list[Atoms], path_xyz: Path) -> None:
    """Write the trajectory to an Extended XYZ file using ASE."""
    path_xyz.parent.mkdir(parents=True, exist_ok=True)
    # ASE writes all frames in a single extxyz file when passed a list.
    ase_write(str(path_xyz), frames, format="extxyz")


def read_trajectory(path_xyz: Path) -> list[Atoms]:
    """Read the trajectory back with mymetal.io.extxyz.extxyz_to_atomlist."""
    return extxyz_to_atomlist(str(path_xyz))


def summarize_trajectory(frames: list[Atoms]) -> dict[str, object]:
    """Deterministic summary used by the printed table and the assertions."""
    formulas = [atoms.get_chemical_formula() for atoms in frames]
    first_atom_pos = np.array([atoms.positions[0] for atoms in frames])
    return {
        "n_frames": len(frames),
        "formula_0": formulas[0],
        "formulas_all_same": len(set(formulas)) == 1,
        "first_atom_pos": first_atom_pos,
        "formulas": formulas,
    }


# ---------------------------------------------------------------------------
# Part B: general_read / general_write
# ---------------------------------------------------------------------------


def build_convergence_table() -> pd.DataFrame:
    """Build a synthetic ENCUT convergence table (deterministic).

    Columns: ``encut`` (eV, 300–550 in 50 eV steps), ``energy`` (eV, converging
    monotonically toward a reference), ``pressure`` (kB, a small oscillating
    residual). All values are deterministic; no randomness is used.
    """
    encut = np.array([300, 350, 400, 450, 500, 550], dtype=float)
    # Energy converges toward -15.6420 eV from above (typical ENCUT curve).
    energy = np.array([-15.4180, -15.5821, -15.6245, -15.6378, -15.6411, -15.6420])
    # Pressure oscillates around 0 with shrinking amplitude.
    pressure = np.array([2.13, -1.07, 0.48, -0.22, 0.09, -0.04])
    return pd.DataFrame(
        {
            "encut": encut,
            "energy": energy,
            "pressure": pressure,
        }
    )


def write_table(df: pd.DataFrame, path_txt: Path) -> None:
    """Write the DataFrame to a formatted text file with general_write."""
    path_txt.parent.mkdir(parents=True, exist_ok=True)
    general_write(filename=str(path_txt), dfc=df, if_write_col_num=True)


def read_table(path_txt: Path) -> pd.DataFrame:
    """Read the formatted text file back with general_read."""
    return general_read(filepath=str(path_txt), has_header=True)


def summarize_table(df_original: pd.DataFrame, df_read: pd.DataFrame) -> dict[str, object]:
    """Deterministic summary of the round-trip."""
    numeric_cols = [c for c in df_original.columns if df_original[c].dtype.kind in "iuf"]
    matches = {
        col: bool(np.allclose(df_original[col].to_numpy(), df_read[col].to_numpy(), atol=1e-6))
        for col in numeric_cols
    }
    return {
        "n_rows_original": len(df_original),
        "n_rows_read": len(df_read),
        "columns": list(df_original.columns),
        "numeric_matches": matches,
        "all_match": all(matches.values()),
    }


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------


def render_figure(
    frames: list[Atoms],
    df_read: pd.DataFrame,
    path_image: Path,
) -> None:
    """Render a 2-panel figure: trajectory (left) and convergence (right).

    Left panel overlays frame 0 (blue) and frame N-1 (orange) atom positions
    and draws the first atom's trajectory as a line with markers.
    Right panel plots ``encut`` vs ``energy`` from the round-tripped DataFrame
    with markers, so the figure visibly exercises both I/O paths.
    """
    path_image.parent.mkdir(parents=True, exist_ok=True)

    general_set_all_rcParams(
        backend="Agg",
        axes_grid=False,
        figure_subp=(1, 2),
        figure_one_figsize=(5.0, 4.4),
        figure_left=0.50,
        figure_top=0.50,
        figure_wspace=0.55,
        figure_hspace=0.45,
        font_family=["DejaVu Sans"],
        fontsize=10,
        axes_titlepad=6,
    )
    fig, (ax_traj, ax_conv) = plt.subplots(1, 2)

    # --- Left panel: trajectory ---
    pos0 = frames[0].positions
    posN = frames[-1].positions
    ax_traj.scatter(pos0[:, 0], pos0[:, 1], c="tab:blue", s=40, label="frame 0", zorder=3)
    ax_traj.scatter(posN[:, 0], posN[:, 1], c="tab:orange", s=40, marker="s",
                    label=f"frame {N_FRAMES - 1}", zorder=3)
    # First-atom trajectory line across all frames.
    first_pos = np.array([atoms.positions[0] for atoms in frames])
    ax_traj.plot(first_pos[:, 0], first_pos[:, 1], "-k", lw=1.2, zorder=2,
                 label="atom 0 trajectory")
    ax_traj.plot(first_pos[:, 0], first_pos[:, 1], "ko", ms=4, zorder=4)
    ax_traj.set_xlabel("x (Å)")
    ax_traj.set_ylabel("y (Å)")
    ax_traj.set_title(f"Cu trajectory ({N_FRAMES} frames, extxyz)")
    ax_traj.legend(loc="best", fontsize=8)
    ax_traj.set_aspect("equal", adjustable="datalim")

    # --- Right panel: convergence table ---
    encut = df_read["encut"].to_numpy(dtype=float)
    energy = df_read["energy"].to_numpy(dtype=float)
    ax_conv.plot(encut, energy, "-o", color="tab:green", lw=1.5, ms=6)
    ax_conv.set_xlabel("ENCUT (eV)")
    ax_conv.set_ylabel("energy (eV)")
    ax_conv.set_title("Convergence (general_write/read round-trip)")
    # Annotate the converged value.
    ax_conv.annotate(
        f"converged: {energy[-1]:.4f} eV",
        xy=(encut[-1], energy[-1]),
        xytext=(0.45, 0.30),
        textcoords="axes fraction",
        arrowprops=dict(arrowstyle="->", color="gray", lw=0.8),
        fontsize=8,
    )

    fig.savefig(path_image, bbox_inches="tight")
    plt.close(fig)

    # Non-blank self-check (mandatory per docs/examples convention).
    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError("image was not created: " + str(path_image))
    image_rgb = plt.imread(path_image)[..., :3]
    if float(np.mean(np.abs(image_rgb - 1.0))) < 0.002:
        raise AssertionError("image is effectively blank: " + str(path_image))


# ---------------------------------------------------------------------------
# Assertions
# ---------------------------------------------------------------------------


def run_assertions(traj_summary: dict[str, object], table_summary: dict[str, object]) -> None:
    """Deterministic invariant checks for both I/O workflows."""
    # Part A assertions.
    assert traj_summary["n_frames"] == N_FRAMES, (
        f"expected {N_FRAMES} frames, got {traj_summary['n_frames']}"
    )
    assert traj_summary["formulas_all_same"], (
        "frames have inconsistent formulas: " + str(traj_summary["formulas"])
    )
    # cubic FCC Cu conventional cell has 4 atoms; 2x2x2 repeat => 32 atoms.
    assert traj_summary["formula_0"] == "Cu32", (
        f"frame 0 formula should be Cu32 for 2x2x2 FCC Cu, got {traj_summary['formula_0']}"
    )
    first_pos = np.asarray(traj_summary["first_atom_pos"])
    assert first_pos.shape == (N_FRAMES, 3), (
        f"first-atom trajectory shape should be ({N_FRAMES}, 3), got {first_pos.shape}"
    )
    # The first atom must actually move across frames (perturbation is nonzero).
    drift = np.linalg.norm(first_pos[-1] - first_pos[0])
    assert drift > 0.0, "first atom did not move across the trajectory"

    # Part B assertions.
    assert table_summary["n_rows_original"] == 6, (
        f"expected 6 rows, got {table_summary['n_rows_original']}"
    )
    assert table_summary["n_rows_read"] == 6, (
        f"expected 6 rows after round-trip, got {table_summary['n_rows_read']}"
    )
    assert list(table_summary["columns"]) == ["encut", "energy", "pressure"], (
        "unexpected columns: " + str(table_summary["columns"])
    )
    assert table_summary["all_match"], (
        "round-trip numeric mismatch: " + str(table_summary["numeric_matches"])
    )


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def run_example(path_output: Path) -> dict[str, object]:
    """Run both I/O workflows, render the figure, and return a combined summary."""
    path_output = path_output.resolve()
    path_output.mkdir(parents=True, exist_ok=True)

    # --- Part A: extxyz trajectory ---
    print("================ Part A: extxyz trajectory I/O ================")
    frames = build_trajectory()
    path_xyz = path_output / "trajectory.xyz"
    write_trajectory(frames, path_xyz)
    frames_read = read_trajectory(path_xyz)
    traj_summary = summarize_trajectory(frames_read)
    print(f"extxyz file written : {path_xyz}")
    print(f"frames read back    : {traj_summary['n_frames']}")
    print(f"frame 0 formula     : {traj_summary['formula_0']}")
    print(f"all formulas same   : {traj_summary['formulas_all_same']}")
    first_pos = np.asarray(traj_summary["first_atom_pos"])
    print("first-atom positions across frames (x, y, z) [Å]:")
    for i, p in enumerate(first_pos):
        print(f"  frame {i}: ({p[0]:.6f}, {p[1]:.6f}, {p[2]:.6f})")

    # --- Part B: general_read / general_write ---
    print("================ Part B: general_read / general_write ================")
    df = build_convergence_table()
    path_txt = path_output / "convergence.txt"
    write_table(df, path_txt)
    df_read = read_table(path_txt)
    table_summary = summarize_table(df, df_read)
    print(f"text file written   : {path_txt}")
    print(f"DataFrame before round-trip ({len(df)} rows):")
    print(df.to_string(index=False))
    print(f"DataFrame after round-trip ({len(df_read)} rows):")
    print(df_read.to_string(index=False))
    print("numeric round-trip match per column:")
    for col, ok in table_summary["numeric_matches"].items():
        print(f"  {col:<10}: {'OK' if ok else 'MISMATCH'}")

    # --- Figure ---
    path_image = path_output / "io_extxyz_and_general.png"
    render_figure(frames_read, df_read, path_image)
    print("wrote: " + str(path_image))

    # --- Assertions ---
    run_assertions(traj_summary, table_summary)

    return {"trajectory": traj_summary, "table": table_summary, "image": path_image}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Demonstrate extxyz_to_atomlist and general_read/general_write I/O "
            "round trips (VASP-free, deterministic)."
        )
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("io-extxyz-and-general-output"),
        help="Directory for the .xyz / .txt / .png outputs.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    summary = run_example(parse_args().output)
    assert summary["trajectory"]["n_frames"] == N_FRAMES
    assert summary["table"]["all_match"]
    print("OK: all assertions passed.")
