#!/usr/bin/env python3
"""Periodic-table heatmap and van Arkel-Ketelaar triangle.

This tutorial is VASP-free. It demonstrates two visualization helpers in
:mod:`mymetal.universal.plot.plotting` that classify and render elemental /
binary-material data without any electronic-structure calculation:

* :func:`periodic_table_heatmap` -- overlays a numeric per-element property
  onto the periodic table. Here we feed it a small synthetic bulk-modulus
  data set (textbook ambient values, GPa) for ~15 pure elements.
* :func:`van_arkel_triangle` -- places binary compounds on the
  electronegativity-based van Arkel-Ketelaar triangle, separating ionic,
  metallic and covalent character. We pass eight binaries (NaCl, MgO,
  Al2O3, SiO2, GaAs, ZnS, CuBr, InP) as element pairs.

The script renders each figure to its own PNG and verifies both are
non-blank, then prints a summary table. Run with::

    python docs/examples/periodic_table_and_arkel.py \
        --output docs/_build/example-periodic-arkel
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # MUST precede import pyplot
import matplotlib.pyplot as plt
import numpy as np

from mymetal.universal.plot.general import general_set_all_rcParams
from mymetal.universal.plot.plotting import (
    periodic_table_heatmap,
    van_arkel_triangle,
)


# ---------------------------------------------------------------------------
# Part A: periodic-table heatmap input
# ---------------------------------------------------------------------------
# Bulk moduli (GPa) of pure elements -- textbook ambient experimental values.
# Synthetic-looking but physically plausible; kept inline so the script is
# fully deterministic and free of any network / DB access.
BULK_MODULUS_GPA: dict[str, int] = {
    "Fe": 211,   # BCC iron
    "Cu": 140,   # FCC copper
    "Al": 76,    # FCC aluminium
    "W": 310,    # BCC tungsten (stiffest pure metal)
    "Cr": 160,   # BCC chromium
    "Ni": 180,   # FCC nickel
    "Mo": 230,   # BCC molybdenum
    "Ti": 110,   # HCP titanium
    "Mg": 45,    # HCP magnesium
    "Au": 137,   # FCC gold
    "Ag": 100,   # FCC silver
    "Pb": 46,    # FCC lead (soft)
    "Zn": 70,    # HCP zinc
    "Si": 100,   # diamond cubic
    "Pt": 230,   # FCC platinum
}


# ---------------------------------------------------------------------------
# Part B: van Arkel triangle input
# ---------------------------------------------------------------------------
# Eight binary/ternary materials, represented as element-symbol pairs. The
# van_arkel_triangle helper only needs the two constituent elements per
# compound (it derives chi from pymatgen Element.X), so ternaries are reduced
# to their two distinct element types. Order: (less electronegative, more).
ARKEL_MATERIALS: list[list[str]] = [
    ["Na", "Cl"],  # NaCl  - rock-salt, strongly ionic
    ["Mg", "O"],   # MgO   - ionic oxide
    ["Al", "O"],   # Al2O3 - corundum, ionic-covalent oxide
    ["Si", "O"],   # SiO2  - covalent-ionic oxide
    ["Ga", "As"],  # GaAs  - III-V semiconductor, covalent
    ["Zn", "S"],   # ZnS   - II-VI, mixed ionic/covalent
    ["Cu", "Br"],  # CuBr  - mixed ionic/covalent
    ["In", "P"],   # InP   - III-V semiconductor, covalent
]

# Display labels (kept parallel to ARKEL_MATERIALS) for the summary table.
ARKEL_LABELS: list[str] = [
    "NaCl", "MgO", "Al2O3", "SiO2", "GaAs", "ZnS", "CuBr", "InP",
]


def build_formula_dict() -> dict[str, int]:
    """Return the synthetic bulk-modulus elemental data (deterministic)."""
    return dict(BULK_MODULUS_GPA)


def build_arkel_materials() -> list[list[str]]:
    """Return the list of binary element pairs for the van Arkel triangle."""
    return [list(pair) for pair in ARKEL_MATERIALS]


def summarize_arkel_row(idx: int, pair: list[str]) -> dict[str, object]:
    """Deterministic summary row: reduced formula + chi mean / difference."""
    from pymatgen.core import Element

    chi_a = float(Element(pair[0]).X)
    chi_b = float(Element(pair[1]).X)
    return {
        "label": ARKEL_LABELS[idx],
        "pair": "-".join(pair),
        "chi_mean": round((chi_a + chi_b) / 2.0, 4),
        "chi_diff": round(abs(chi_a - chi_b), 4),
    }


def render_periodic_heatmap(
    formula_dict: dict[str, int], path_image: Path
) -> None:
    """Render the periodic-table heatmap PNG and verify it is non-blank.

    periodic_table_heatmap (matplotlib backend, no pymatviz) creates its own
    figure via plt.subplots(); we capture plt.gcf() and save it.
    """
    path_image.parent.mkdir(parents=True, exist_ok=True)
    general_set_all_rcParams(backend="Agg", font_family=["DejaVu Sans"], fontsize=10)

    periodic_table_heatmap(
        elemental_data=formula_dict,
        cbar_label="Bulk modulus (GPa)",
        cbar_label_size=14,
        cmap="YlOrRd",
        value_format="%.0f",
        value_fontsize=9,
        symbol_fontsize=12,
        pymatviz=False,  # force the matplotlib static backend -> PNG-friendly
    )
    fig = plt.gcf()
    fig.savefig(path_image, bbox_inches="tight")
    plt.close(fig)
    _assert_image_nonblank(path_image, "periodic heatmap")


def render_arkel_triangle(
    materials: list[list[str]], path_image: Path
) -> None:
    """Render the van Arkel-Ketelaar triangle PNG and verify non-blank.

    van_arkel_triangle draws onto the current axes (plt.gca()) and does NOT
    create its own figure, so we create one first.
    """
    path_image.parent.mkdir(parents=True, exist_ok=True)
    general_set_all_rcParams(backend="Agg", font_family=["DejaVu Sans"], fontsize=10)

    fig, ax = plt.subplots(figsize=(8, 7))
    van_arkel_triangle(materials, annotate=True)
    fig.savefig(path_image, bbox_inches="tight")
    plt.close(fig)
    _assert_image_nonblank(path_image, "van Arkel triangle")


def _assert_image_nonblank(path_image: Path, label: str) -> None:
    """Read back a saved PNG and assert it exists, is non-trivial, non-blank."""
    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError(f"{label} image was not created: {path_image}")
    image_rgb = plt.imread(path_image)[..., :3]
    deviation = float(np.mean(np.abs(image_rgb - 1.0)))
    if deviation < 0.002:
        raise AssertionError(f"{label} image is effectively blank: {path_image}")


def run_example(path_output: Path) -> list[dict[str, object]]:
    """Build inputs, render both figures, print a summary, return rows."""
    path_output = path_output.resolve()
    path_output.mkdir(parents=True, exist_ok=True)

    # --- Part A: periodic table heatmap ----------------------------------
    formula_dict = build_formula_dict()
    print("Part A: periodic_table_heatmap input (bulk modulus, GPa):")
    for el, val in formula_dict.items():
        print(f"  {el:>3}: {val}")
    assert len(formula_dict) >= 15, "formula_dict must have >= 15 entries"
    path_heatmap = path_output / "periodic_table_heatmap.png"
    render_periodic_heatmap(formula_dict, path_heatmap)
    print("wrote: " + str(path_heatmap))

    # --- Part B: van Arkel triangle --------------------------------------
    materials = build_arkel_materials()
    print("\nPart B: van_arkel_triangle input (element pairs):")
    rows: list[dict[str, object]] = []
    for idx, pair in enumerate(materials):
        row = summarize_arkel_row(idx, pair)
        rows.append(row)
        print(f"  {row['label']:<6} ({row['pair']:<5})  "
              f"chi_mean={row['chi_mean']:.4f}  chi_diff={row['chi_diff']:.4f}")
    assert len(materials) == 8, "expected 8 van Arkel materials"
    path_arkel = path_output / "van_arkel_triangle.png"
    render_arkel_triangle(materials, path_arkel)
    print("wrote: " + str(path_arkel))

    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render a periodic-table heatmap and a van Arkel triangle."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("example-periodic-arkel-output"),
        help="Directory for the two output PNGs.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    arkel_rows = run_example(parse_args().output)
    # Deterministic assertions on the returned rows.
    assert len(arkel_rows) == 8
    # NaCl must have the largest chi difference (most ionic) of this set.
    max_diff = max(float(r["chi_diff"]) for r in arkel_rows)
    assert abs(max_diff - float(arkel_rows[0]["chi_diff"])) < 1e-9, (
        "NaCl should be the most ionic (largest chi_diff) in this set"
    )
    # Every material must have a positive chi difference (two distinct elements).
    assert all(float(r["chi_diff"]) > 0.0 for r in arkel_rows)
    print("OK: all assertions passed.")
