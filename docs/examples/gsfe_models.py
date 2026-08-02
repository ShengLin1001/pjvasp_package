#!/usr/bin/env python3
"""Build and compare GSFE (generalized stacking fault energy) supercells.

This tutorial is VASP-free and deterministic. It demonstrates
:func:`mymetal.build.bulk.gsfe.create_gsfe_model`, the one-stop builder for
generalized-stacking-fault supercells, on the three pure-ASE slip geometries
that do not require a running VASP environment:

* ``FCC_111``  -> Au, FCC (111) slip plane, Shockley partial system
* ``HCP_basal``  -> Mg, HCP (0001) basal plane
* ``HCP_prism1w`` -> Mg, HCP (10-10) prism I, wide spacing

The other ``gsfe_type`` options (``FCC_100``, ``HCP_pyr1w``, ``HCP_pyr2``) route
through ``vasp_create_*`` helpers that depend on :mod:`myvasp`, so they are
intentionally out of scope for this VASP-free example.

Each model is built, its chemical symbols relabelled to the physical element
(``create_gsfe_model`` defaults to ``Au``), a summary table is printed, and a
1x3 side-view comparison figure is rendered for the documentation page.
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

from mymetal.build.bulk.gsfe import create_gsfe_model
from mymetal.universal.plot.general import general_set_all_rcParams


# Lattice constants (Å). Textbook ambient values; the script is fully
# deterministic and needs no network access.
A_FCC_AU = 4.08
A_HCP_MG = 3.21
C_HCP_MG = 5.21

# Supercell sizes. ``create_gsfe_model`` accepts ``np.array``; the z-repeat sets
# how many fault layers the slab contains. These are the builder's documented
# defaults, restated here so the example is self-describing.
SIZE_FCC_111 = np.array([1, 1, 6])
SIZE_HCP_BASAL = np.array([1, 1, 6])
SIZE_HCP_PRISM1W = np.array([1, 1, 6])


def build_gsfe_models() -> list[Atoms]:
    """Return the three GSFE supercells in display order, relabelled.

    ``create_gsfe_model`` does not expose a ``symbol`` argument, so it always
    builds ``Au`` atoms. We relabel post-hoc to the physical element of each
    system (Au for the FCC case, Mg for the two HCP cases). Relabelling is a
    metadata-only operation: it changes neither positions nor the cell.
    """
    specs = [
        ("FCC_111", A_FCC_AU, None, SIZE_FCC_111, "Au"),
        ("HCP_basal", A_HCP_MG, C_HCP_MG, SIZE_HCP_BASAL, "Mg"),
        ("HCP_prism1w", A_HCP_MG, C_HCP_MG, SIZE_HCP_PRISM1W, "Mg"),
    ]
    cells: list[Atoms] = []
    for gsfe_type, a, c, size, symbol in specs:
        atoms = create_gsfe_model(gsfe_type=gsfe_type, a=a, c=c, size=size)
        atoms.set_chemical_symbols([symbol] * len(atoms))
        cells.append(atoms)
    return cells


def cell_label(idx: int) -> str:
    labels = [
        "Au FCC(111)",
        "Mg HCP(0001) basal",
        "Mg HCP(10-10) prism I (wide)",
    ]
    return labels[idx]


def summarize_cell(idx: int, cell: Atoms) -> dict[str, object]:
    """Return a deterministic summary used by the table and the assertions."""
    lengths = np.linalg.norm(cell.cell.array, axis=1)
    angles = cell.cell.angles()
    return {
        "label": cell_label(idx),
        "gsfe_type": ["FCC_111", "HCP_basal", "HCP_prism1w"][idx],
        "formula": cell.get_chemical_formula(),
        "atoms": len(cell),
        "a_A": round(float(lengths[0]), 4),
        "b_A": round(float(lengths[1]), 4),
        "c_A": round(float(lengths[2]), 4),
        "alpha_deg": round(float(angles[0]), 4),
        "beta_deg": round(float(angles[1]), 4),
        "gamma_deg": round(float(angles[2]), 4),
        "volume_A3": round(float(cell.get_volume()), 4),
        "pbc": cell.pbc.tolist(),
    }


def render_comparison(cells: list[Atoms], path_image: Path) -> None:
    """Render a 1x3 side-view comparison figure with the unit cell shown."""
    path_image.parent.mkdir(parents=True, exist_ok=True)

    general_set_all_rcParams(
        backend="Agg",
        axes_grid=False,
        figure_subp=(1, 3),
        figure_one_figsize=(3.6, 4.4),
        figure_left=0.50,
        figure_top=0.50,
        figure_wspace=0.55,
        figure_hspace=0.45,
        font_family=["DejaVu Sans"],
        fontsize=10,
        axes_titlepad=6,
    )
    fig, axes = plt.subplots(1, 3)
    for col, cell in enumerate(cells):
        # Repeat along x so the side view shows a couple of in-plane periods.
        display = cell.repeat((2, 1, 1))
        plot_atoms(
            display,
            axes[col],
            rotation="90x,0y,0z",
            show_unit_cell=1,
            radii=0.6,
        )
        axes[col].set_title(cell_label(col) + "\nside view")
        axes[col].set_axis_off()
    fig.savefig(path_image, bbox_inches="tight")
    plt.close(fig)

    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError("comparison image was not created: " + str(path_image))
    image_rgb = plt.imread(path_image)[..., :3]
    if float(np.mean(np.abs(image_rgb - 1.0))) < 0.002:
        raise AssertionError("comparison image is effectively blank: " + str(path_image))


def run_example(path_output: Path) -> list[dict[str, object]]:
    """Build, render, and return the per-model summary rows."""
    path_output = path_output.resolve()
    path_output.mkdir(parents=True, exist_ok=True)

    cells = build_gsfe_models()
    rows: list[dict[str, object]] = []
    for idx, cell in enumerate(cells):
        rows.append(summarize_cell(idx, cell))

    path_image = path_output / "gsfe_models.png"
    render_comparison(cells, path_image)

    print("label                          gsfe_type      formula  atoms  "
          "a_A     b_A     c_A      alpha   beta    gamma    volume_A3  pbc")
    for row in rows:
        print(
            f"{row['label']:<30}  {row['gsfe_type']:<13}  {row['formula']:<7}  "
            f"{row['atoms']:<5}  "
            f"{row['a_A']:<6.4f}  {row['b_A']:<6.4f}  {row['c_A']:<7.4f}  "
            f"{row['alpha_deg']:<6.2f}  {row['beta_deg']:<6.2f}  "
            f"{row['gamma_deg']:<7.2f}  {row['volume_A3']:<9.4f}  {row['pbc']}"
        )
    print("wrote: " + str(path_image))
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build and compare VASP-free GSFE supercells (FCC_111, "
                    "HCP_basal, HCP_prism1w).")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("gsfe-models-output"),
        help="Directory for the comparison PNG.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    rows = run_example(parse_args().output)
    # Deterministic assertions.
    assert all(row["atoms"] > 0 for row in rows), rows
    assert all(row["pbc"] == [True, True, True] for row in rows), rows
    # Formulas must match the relabelled elements.
    assert rows[0]["formula"].startswith("Au"), rows[0]
    assert rows[1]["formula"].startswith("Mg"), rows[1]
    assert rows[2]["formula"].startswith("Mg"), rows[2]
