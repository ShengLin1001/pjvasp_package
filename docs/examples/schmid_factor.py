#!/usr/bin/env python3
"""Compute and visualize FCC Schmid factors for arbitrary loading directions.

This tutorial is VASP-free. It demonstrates
:func:`mymetal.calculate.material_science.schmid.cal_fcc_schmid_factors`:

1. Pick a loading direction (e.g. ``[1, 1, 6]``) and compute the Schmid
   factor of every FCC slip system ``{111}<110>``.
2. Scan a set of standard tensile directions and report the maximum
   resolved Schmid factor for each, which is the quantity plasticity
   simulations care about.
3. Render a two-panel figure: a polar bar plot of the maximum Schmid
   factor along stereographic directions, and a horizontal bar chart of
   all activated slip systems for ``[1, 1, 6]``.

The script is fully deterministic. No DFT, no atomistic relaxation, no
external data file is required.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mymetal.calculate.material_science.schmid import cal_fcc_schmid_factors
from mymetal.universal.plot.general import general_set_all_rcParams


# Reference loading direction reused throughout the tutorial. ``[1, 1, 6]``
# is a classic example from FCC single-crystal deformation studies because
# it activates multiple slip systems with comparable Schmid factors.
REFERENCE_DIRECTION = [1, 1, 6]

# Standard tensile-axis directions used in the stereographic scan. Each
# entry is a crystallographic direction ``[u, v, w]`` in the cubic basis.
SCAN_DIRECTIONS = [
    [1, 0, 0],
    [1, 1, 0],
    [1, 1, 1],
    [1, 1, 2],
    [1, 1, 3],
    [1, 1, 6],
    [1, 2, 3],
    [1, 3, 4],
    [2, 1, 3],
    [3, 2, 1],
]

# Number of polar grid points used for the continuous stereographic scan.
POLAR_RESOLUTION = 36


def compute_reference_table() -> pd.DataFrame:
    """Return the per-slip-system Schmid factors for the reference direction."""
    return cal_fcc_schmid_factors(normal_orientation=REFERENCE_DIRECTION)


def compute_scan_table() -> pd.DataFrame:
    """Return the maximum Schmid factor per direction in ``SCAN_DIRECTIONS``."""
    rows: list[dict[str, object]] = []
    for direction in SCAN_DIRECTIONS:
        df = cal_fcc_schmid_factors(normal_orientation=direction)
        rows.append({
            "direction": direction,
            "max_schmid": float(df["schmid_factor"].max()),
            "n_active": int((df["schmid_factor"] > 1e-6).sum()),
        })
    return pd.DataFrame(rows)


def compute_polar_grid() -> pd.DataFrame:
    """Return a dense grid of (theta, r, schmid) values for the polar plot.

    ``theta`` is the in-plane angle (radians) of the loading direction
    projected onto the (001) stereographic triangle, and ``r`` is the
    radial coordinate along the [001]-[011]-[111] boundary.
    """
    rows: list[dict[str, float]] = []
    for i_theta in range(POLAR_RESOLUTION):
        theta = 2.0 * np.pi * i_theta / POLAR_RESOLUTION
        for i_r in range(1, POLAR_RESOLUTION + 1):
            r = float(i_r) / POLAR_RESOLUTION
            direction = [
                int(round(10.0 * r * np.cos(theta))),
                int(round(10.0 * r * np.sin(theta))),
                10,
            ]
            if direction == [0, 0, 0]:
                continue
            df = cal_fcc_schmid_factors(normal_orientation=direction)
            rows.append({
                "theta": theta,
                "r": r,
                "schmid": float(df["schmid_factor"].max()),
            })
    return pd.DataFrame(rows)


def render_figure(
    reference_table: pd.DataFrame,
    scan_table: pd.DataFrame,
    polar_grid: pd.DataFrame,
    path_image: Path,
) -> None:
    """Render the two-panel Schmid-factor figure."""
    path_image.parent.mkdir(parents=True, exist_ok=True)

    general_set_all_rcParams(
        backend="Agg",
        axes_grid=True,
        figure_subp=(1, 2),
        figure_one_figsize=(6.0, 5.0),
        figure_left=0.65,
        figure_top=0.55,
        figure_wspace=1.0,
        figure_hspace=0.6,
        font_family=["DejaVu Sans"],
        fontsize=11,
        axes_titlepad=6,
    )
    fig, axes = plt.subplots(1, 2, subplot_kw={"projection": None})

    # Left: polar bar of maximum Schmid factor on a (001) stereographic
    # triangle. Each bar is colored by its Schmid value.
    ax = axes[0]
    ax.remove()
    ax = fig.add_subplot(1, 2, 1, projection="polar")
    theta = polar_grid["theta"].to_numpy()
    radii = polar_grid["r"].to_numpy()
    schmid = polar_grid["schmid"].to_numpy()
    sc = ax.scatter(theta, radii, c=schmid, cmap="viridis",
                    s=18, vmin=0.0, vmax=0.5)
    ax.set_theta_zero_location("E")
    ax.set_theta_direction(1)
    ax.set_rlim(0.0, 1.0)
    ax.set_rticks([0.25, 0.5, 0.75, 1.0])
    ax.set_title("Max Schmid factor\n(001) stereographic scan", pad=12)
    cb = fig.colorbar(sc, ax=ax, orientation="vertical", pad=0.10, shrink=0.8)
    cb.set_label("Schmid factor")

    # Right: per-slip-system Schmid factors for the reference direction.
    ax = axes[1]
    df = reference_table.copy()
    df = df.sort_values("schmid_factor", ascending=True).reset_index(drop=True)
    y_pos = np.arange(len(df))
    labels = [
        f"plane {tuple(row['slip_plane'])}\nb {tuple(np.round(row['dislocation'], 2))}"
        for _, row in df.iterrows()
    ]
    colors = ["#d62728" if v > 1e-6 else "#cccccc" for v in df["schmid_factor"]]
    ax.barh(y_pos, df["schmid_factor"], color=colors)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Schmid factor")
    ax.set_title(f"All FCC slip systems\nloading direction = {REFERENCE_DIRECTION}")
    ax.axvline(0.0, color="gray", linestyle=":", linewidth=0.8)
    ax.set_xlim(0.0, 0.5)

    fig.savefig(path_image, bbox_inches="tight")
    plt.close(fig)

    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError("schmid image was not created: " + str(path_image))


def run_example(path_output: Path):
    """Compute Schmid factors and render the figure."""
    path_output = path_output.resolve()
    path_output.mkdir(parents=True, exist_ok=True)

    reference_table = compute_reference_table()
    scan_table = compute_scan_table()
    polar_grid = compute_polar_grid()

    path_image = path_output / "schmid_factor.png"
    render_figure(reference_table, scan_table, polar_grid, path_image)

    print(f"reference direction: {REFERENCE_DIRECTION}")
    print(f"  total slip systems: {len(reference_table)}")
    print(f"  active (m > 1e-6):   {(reference_table['schmid_factor'] > 1e-6).sum()}")
    print(f"  max Schmid factor:  {float(reference_table['schmid_factor'].max()):.6f}")
    print("")
    print("scan summary:")
    print("direction       max_schmid  n_active")
    for _, row in scan_table.iterrows():
        d = row["direction"]
        print(f"  [{d[0]} {d[1]} {d[2]}]       "
              f"{row['max_schmid']:.6f}    {int(row['n_active'])}")
    print("wrote: " + str(path_image))
    return {
        "reference_table": reference_table,
        "scan_table": scan_table,
        "polar_grid": polar_grid,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute and visualize FCC Schmid factors.")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("schmid-output"),
        help="Directory for the Schmid PNG.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    result = run_example(parse_args().output)
    # The reference direction must activate at least one slip system.
    assert (result["reference_table"]["schmid_factor"] > 1e-6).any()
    # Schmid factors must lie in [0, 0.5].
    assert result["reference_table"]["schmid_factor"].between(0.0, 0.5).all()
    # [1, 1, 1] is a "hard" axis: it must NOT activate slip on the (1, 1, 1)
    # plane itself, because the loading direction is parallel to the plane
    # normal. Other {111} planes are still activated.
    df_111 = cal_fcc_schmid_factors(normal_orientation=[1, 1, 1])
    mask_111_plane = df_111["slip_plane"].apply(
        lambda v: list(v) == [1, 1, 1])
    assert df_111.loc[mask_111_plane, "schmid_factor"].max() < 1e-6
    # The maximum Schmid factor along [1, 1, 6] must be > 0.4.
    df_116 = cal_fcc_schmid_factors(normal_orientation=[1, 1, 6])
    assert df_116["schmid_factor"].max() > 0.4
    # [1, 1, 6] must have a higher max Schmid factor than [1, 1, 1].
    assert df_116["schmid_factor"].max() > df_111["schmid_factor"].max()
