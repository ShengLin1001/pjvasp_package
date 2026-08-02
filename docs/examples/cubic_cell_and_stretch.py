#!/usr/bin/env python3
"""Find the orthorhombic (cubic) xy cell of an HCP film and stretch it.

This tutorial is VASP-free and deterministic. It demonstrates two film-building
helpers in :mod:`mymetal.build.film`:

* :func:`mymetal.build.film.findcubic.find_cubic` — convert the hexagonal
  primitive cell of an HCP (or FCC) film into an orthorhombic xy cell that
  preserves the xy-projected area (``a = sqrt(area/sqrt(3))``,
  ``b = sqrt(3) * a``, gamma = 90 deg).
* :func:`mymetal.build.film.stretch.stretch_along_direction_to_cell` — apply a
  uniaxial stretch factor to a chosen cell direction (here ``'x'``).

The script produces one figure with a 2x2 grid of panels:

* top-left  — primitive HCP Mg film (hexagonal xy cell, gamma = 120 deg)
* top-right — orthorhombic (cubic) xy cell of the same film, gamma = 90 deg
* bottom-left  — orthorhombic film stretched along x by 0.99
* bottom-right — orthorhombic film stretched along x by 1.01

Run with::

    python docs/examples/cubic_cell_and_stretch.py --output docs/_build/example-cubic-stretch
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ase import Atoms
from ase.visualize.plot import plot_atoms

from mymetal.build.film.extrfilm import cal_area
from mymetal.build.film.findcubic import find_cubic
from mymetal.build.film.stretch import generate_film, stretch_along_direction_to_cell
from mymetal.universal.plot.general import general_set_all_rcParams


# Textbook ambient Mg hcp lattice constant (Å). Hard-coded so the script is
# fully deterministic and free of network access.
A_HCP_MG = 3.21
MY_VACUUM = 10.0
NUM_LAYERS = 4

# Stretch factors for Panel B. 1.00 is included in the printed summary but is
# not drawn (it is identical to the orthorhombic panel already shown).
STRETCH_FACTORS = [0.99, 1.00, 1.01]
STRETCH_DIRECTION = "x"


def build_primitive_film() -> Atoms:
    """Build the primitive HCP Mg (0001) film used as the starting point."""
    film = generate_film(
        symbols="Mg",
        structure="hcp",
        num_layers=NUM_LAYERS,
        my_vacuum=MY_VACUUM,
        slice_plane=(0, 0, 1),
        a_hcp=A_HCP_MG,
    )
    return film


def build_orthorhombic_film(prim: Atoms) -> Atoms:
    """Apply find_cubic to obtain the orthorhombic (cubic) xy cell."""
    return find_cubic(prim, type="hcp")


def build_stretched_films(ortho: Atoms) -> list[Atoms]:
    """Stretch the orthorhombic film along x by each factor in STRETCH_FACTORS."""
    stretched: list[Atoms] = []
    for factor in STRETCH_FACTORS:
        stretched.append(
            stretch_along_direction_to_cell(
                ortho,
                stretch_factor=factor,
                stretch_direction=STRETCH_DIRECTION,
            )
        )
    return stretched


def cell_summary(atoms: Atoms) -> dict[str, object]:
    """Return a deterministic per-cell summary used by the table and assertions."""
    cell = atoms.cell
    lengths = np.linalg.norm(cell.array, axis=1)
    angles = cell.angles()
    return {
        "formula": atoms.get_chemical_formula(),
        "atoms": len(atoms),
        "a_A": round(float(lengths[0]), 6),
        "b_A": round(float(lengths[1]), 6),
        "c_A": round(float(lengths[2]), 6),
        "alpha_deg": round(float(angles[0]), 4),
        "beta_deg": round(float(angles[1]), 4),
        "gamma_deg": round(float(angles[2]), 4),
        "area_A2": round(float(cal_area(atoms)), 6),
    }


def print_summary_table(
    prim: Atoms,
    ortho: Atoms,
    stretched: list[Atoms],
) -> list[dict[str, object]]:
    """Print a human-readable summary table and return the row dicts."""
    rows: list[dict[str, object]] = []

    prim_row = cell_summary(prim)
    prim_row["label"] = "primitive (hex)"
    rows.append(prim_row)

    ortho_row = cell_summary(ortho)
    ortho_row["label"] = "orthorhombic (cubic)"
    rows.append(ortho_row)

    for factor, atoms in zip(STRETCH_FACTORS, stretched):
        row = cell_summary(atoms)
        row["label"] = f"stretch x * {factor:.2f}"
        rows.append(row)

    header = (
        "label                       formula  atoms  a_A      b_A      c_A      "
        "alpha   beta    gamma   area_A2"
    )
    print(header)
    for row in rows:
        print(
            f"{row['label']:<26}  {row['formula']:<7}  "
            f"{row['atoms']:<5}  "
            f"{row['a_A']:<7.4f}  "
            f"{row['b_A']:<7.4f}  "
            f"{row['c_A']:<7.4f}  "
            f"{row['alpha_deg']:<6.2f}  "
            f"{row['beta_deg']:<6.2f}  "
            f"{row['gamma_deg']:<6.2f}  "
            f"{row['area_A2']:<8.4f}"
        )
    return rows


def render_figure(
    prim: Atoms,
    ortho: Atoms,
    stretched: list[Atoms],
    path_image: Path,
) -> None:
    """Render the 2x2 comparison figure and verify it is non-blank."""
    path_image.parent.mkdir(parents=True, exist_ok=True)

    general_set_all_rcParams(
        backend="Agg",
        axes_grid=False,
        figure_subp=(2, 2),
        figure_one_figsize=(4.2, 4.6),
        figure_left=0.50,
        figure_top=0.50,
        figure_wspace=0.55,
        figure_hspace=0.45,
        font_family=["DejaVu Sans"],
        fontsize=10,
        axes_titlepad=6,
    )
    fig, axes = plt.subplots(2, 2)

    # Panel A (top row): primitive hexagonal vs orthorhombic cubic xy cell.
    # Use a slight top-down-ish rotation so the xy hexagon / rectangle is
    # clearly visible while the slab thickness is still hinted at.
    panel_specs = [
        (axes[0, 0], prim, "primitive HCP Mg film\n(hexagonal xy, gamma=120)"),
        (axes[0, 1], ortho, "orthorhombic xy cell\n(cubic, gamma=90)"),
        # Panel B (bottom row): two non-trivial stretches along x.
        (axes[1, 0], stretched[0], "stretch x * 0.99\n(compressed)"),
        (axes[1, 1], stretched[2], "stretch x * 1.01\n(expanded)"),
    ]
    for ax, atoms, title in panel_specs:
        plot_atoms(
            atoms,
            ax,
            rotation="5x,0y,0z",
            show_unit_cell=1,
            radii=0.6,
        )
        ax.set_title(title)
        ax.set_axis_off()

    fig.savefig(path_image, bbox_inches="tight")
    plt.close(fig)

    # Verify the image exists, is non-trivially sized, and is not near-blank.
    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError("image was not created: " + str(path_image))
    image_rgb = plt.imread(path_image)[..., :3]
    if float(np.mean(np.abs(image_rgb - 1.0))) < 0.002:
        raise AssertionError("image is effectively blank: " + str(path_image))


def run_assertions(
    prim: Atoms,
    ortho: Atoms,
    stretched: list[Atoms],
    rows: list[dict[str, object]],
) -> None:
    """Deterministic checks that the documented invariants hold."""

    # 1. Area is preserved by find_cubic (orthorhombic xy cell has the same
    #    xy-projected area as the primitive hexagonal cell).
    prim_area = cal_area(prim)
    ortho_area = cal_area(ortho)
    assert abs(prim_area - ortho_area) < 1e-6, (
        f"find_cubic changed the xy area: {prim_area} -> {ortho_area}"
    )

    # 2. The orthorhombic cell is genuinely orthorhombic in xy: gamma == 90.
    assert abs(rows[1]["gamma_deg"] - 90.0) < 1e-3, (
        f"orthorhombic gamma != 90: {rows[1]['gamma_deg']}"
    )
    # And the primitive hexagonal cell has gamma == 120.
    assert abs(rows[0]["gamma_deg"] - 120.0) < 1e-3, (
        f"primitive gamma != 120: {rows[0]['gamma_deg']}"
    )

    # 3. find_cubic satisfies b = sqrt(3) * a (orthorhombic HCP convention).
    a_ortho = rows[1]["a_A"]
    b_ortho = rows[1]["b_A"]
    assert abs(b_ortho - np.sqrt(3.0) * a_ortho) < 1e-3, (
        f"b != sqrt(3)*a: b={b_ortho}, sqrt(3)*a={np.sqrt(3.0)*a_ortho}"
    )

    # 4. Stretch along x scales only the x-component of the cell vectors
    #    (temp[:,0] *= factor) and leaves y, z untouched. Because the
    #    orthorhombic cell is axis-aligned, the cell parameter a scales by
    #    exactly the factor while b and c are unchanged.
    ortho_a = rows[1]["a_A"]
    ortho_b = rows[1]["b_A"]
    ortho_c = rows[1]["c_A"]
    for factor, atoms, row in zip(STRETCH_FACTORS, stretched, rows[2:]):
        expected_a = round(float(ortho_a * factor), 6)
        assert abs(row["a_A"] - expected_a) < 1e-4, (
            f"stretch x * {factor}: a={row['a_A']} != expected {expected_a}"
        )
        assert abs(row["b_A"] - ortho_b) < 1e-4, (
            f"stretch x * {factor}: b changed ({row['b_A']} != {ortho_b})"
        )
        assert abs(row["c_A"] - ortho_c) < 1e-4, (
            f"stretch x * {factor}: c changed ({row['c_A']} != {ortho_c})"
        )
        # Area scales by the stretch factor since only x moves.
        expected_area = round(float(ortho_area * factor), 6)
        assert abs(row["area_A2"] - expected_area) < 1e-4, (
            f"stretch x * {factor}: area={row['area_A2']} != expected {expected_area}"
        )


def run_example(path_output: Path) -> list[dict[str, object]]:
    """Build the films, render the figure, run assertions, return summary rows."""
    path_output = path_output.resolve()
    path_output.mkdir(parents=True, exist_ok=True)

    prim = build_primitive_film()
    ortho = build_orthorhombic_film(prim)
    stretched = build_stretched_films(ortho)

    rows = print_summary_table(prim, ortho, stretched)

    path_image = path_output / "cubic_cell_and_stretch.png"
    render_figure(prim, ortho, stretched, path_image)

    run_assertions(prim, ortho, stretched, rows)

    print("wrote: " + str(path_image))
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Find the orthorhombic (cubic) xy cell of an HCP Mg film and "
            "stretch it along x. VASP-free, deterministic."
        )
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("cubic-stretch-output"),
        help="Directory for the comparison PNG.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    rows = run_example(parse_args().output)
    # Sanity: every row has a positive atom count and a positive area.
    assert all(row["atoms"] > 0 for row in rows)
    assert all(row["area_A2"] > 0 for row in rows)
    print("OK: all assertions passed.")
