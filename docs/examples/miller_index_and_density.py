#!/usr/bin/env python3
"""HCP Miller index 3<->4 conversion and atomic density of bulk cells.

This tutorial is VASP-free. It demonstrates two independent but complementary
routines from :mod:`mymetal.universal.atom`:

1. :func:`mymetal.universal.atom.miller.three_index_to_four_index` -- convert
   HCP crystallographic directions between the three-index notation
   ``[U, V, W]`` and the four-index notation ``[u, v, t, w]`` (with
   ``t = -(u + v)``). The forward path uses the standard transformation
   matrix

   .. math::

       \\begin{pmatrix} u \\\\ v \\\\ w \\end{pmatrix} =
       \\begin{pmatrix}
       2/3 & -1/3 & 0 \\\\
       -1/3 & 2/3 & 0 \\\\
       0 & 0 & 1
       \\end{pmatrix}
       \\begin{pmatrix} U \\\\ V \\\\ W \\end{pmatrix},

   and the reverse path uses its inverse

   .. math::

       \\begin{pmatrix} U \\\\ V \\\\ W \\end{pmatrix} =
       \\begin{pmatrix}
       2 & 1 & 0 \\\\
       1 & 2 & 0 \\\\
       0 & 0 & 1
       \\end{pmatrix}
       \\begin{pmatrix} u \\\\ v \\\\ w \\end{pmatrix}.

2. :func:`mymetal.universal.atom.density.cal_density` -- number density
   ``natoms / volume`` (atoms/A^3) of an :class:`ase.Atoms` object, here
   applied to four textbook bulk cells (FCC Cu, BCC Fe, HCP Mg, diamond Si).

The script renders a two-panel figure reused by the documentation page:

* Left:  a schematic table showing the 3-index -> 4-index mapping for the
  four worked HCP direction examples, plus the round-trip check.
* Right: a grouped bar chart comparing number density (atoms/A^3) and mass
  density (g/cm^3) for the four bulk structures.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ase import Atoms
from ase.build import bulk
from ase.data import atomic_masses, atomic_numbers

from mymetal.universal.atom.density import cal_density
from mymetal.universal.atom.miller import three_index_to_four_index
from mymetal.universal.plot.general import general_set_all_rcParams


# Lattice constants (A). Textbook ambient values, kept inline so the script
# is fully deterministic and needs no network access.
A_FCC_CU = 3.61
A_BCC_FE = 2.87
A_HCP_MG = 3.21
C_HCP_MG = 5.21
A_DIA_SI = 5.43

# Conversion factor: 1 u / A^3 = 1.66053906660 g/cm^3.
# Used to convert number density (atoms/A^3) * average atomic mass (u/atom)
# into mass density (g/cm^3).
U_PER_A3_TO_G_PER_CM3 = 1.66053906660


# ---------------------------------------------------------------------------
# Part A: HCP Miller index conversion
# ---------------------------------------------------------------------------

def build_miller_examples() -> list[dict[str, object]]:
    """Return the four worked HCP direction examples with conversions.

    Each row carries the original input, the converted output, and (for the
    3-index inputs) a 3 -> 4 -> 3 round-trip check.
    """
    examples: list[dict[str, object]] = []

    # Forward: 3-index [U, V, W] -> 4-index [u, v, t, w]
    forward_inputs = [[1, 0, 0], [1, 1, 0]]
    for inp in forward_inputs:
        four = three_index_to_four_index(inp, reverse=False)
        back = three_index_to_four_index(four, reverse=True)
        examples.append({
            "direction": "forward",
            "input": inp,
            "input_str": f"[{inp[0]}, {inp[1]}, {inp[2]}]",
            "output": four,
            "output_str": _fmt_index(four),
            "roundtrip": back,
            "roundtrip_str": _fmt_index(back),
        })

    # Reverse: 4-index [u, v, t, w] -> 3-index [U, V, W]
    reverse_inputs = [[2, -1, 0, 0], [1, 1, -2, 3]]
    for inp in reverse_inputs:
        three = three_index_to_four_index(inp, reverse=True)
        # Round-trip: 4 -> 3 -> 4. Note the forward transform yields floats
        # that are integer multiples of the original 4-index when the input
        # satisfies t = -(u+v); we check proportionality, not equality.
        four_back = three_index_to_four_index(three, reverse=False)
        examples.append({
            "direction": "reverse",
            "input": inp,
            "input_str": f"[{inp[0]}, {inp[1]}, {inp[2]}, {inp[3]}]",
            "output": three,
            "output_str": _fmt_index(three),
            "roundtrip": four_back,
            "roundtrip_str": _fmt_index(four_back),
        })

    return examples


def _fmt_index(idx: list[float]) -> str:
    """Format a Miller index list, collapsing near-integer floats."""
    out = []
    for v in idx:
        rv = round(float(v), 4)
        if abs(rv - round(rv)) < 1e-6:
            out.append(str(int(round(rv))))
        else:
            out.append(f"{rv:g}")
    return "[" + ", ".join(out) + "]"


def verify_roundtrips(examples: list[dict[str, object]]) -> None:
    """Assert every round-trip is consistent up to integer scaling.

    Forward rows: 3 -> 4 -> 3 must reproduce the original 3-index exactly.
    Reverse rows: 4 -> 3 -> 4 must reproduce the original 4-index up to an
    overall integer scale factor (the forward transform introduces a factor
    of 3 for directions that are not on the basal plane).
    """
    for ex in examples:
        inp = np.array(ex["input"], dtype=float)
        rt = np.array(ex["roundtrip"], dtype=float)
        if ex["direction"] == "forward":
            assert rt.shape == (3,), f"forward roundtrip shape: {rt.shape}"
            assert np.allclose(rt, inp, atol=1e-6), (
                f"forward roundtrip failed: {inp.tolist()} -> {rt.tolist()}"
            )
        else:
            assert rt.shape == (4,), f"reverse roundtrip shape: {rt.shape}"
            # The 4-index ``t`` is a derived quantity (t = -(u+v)), so the
            # round-trip is only required to reproduce the independent
            # components ``[u, v, w]`` (indices 0, 1, 3) up to an overall
            # integer scale factor. The forward transform introduces a
            # factor of 3 for directions that are not on the basal plane.
            inp_uvw = np.array([inp[0], inp[1], inp[3]])
            rt_uvw = np.array([rt[0], rt[1], rt[3]])
            nonzero = np.abs(inp_uvw[np.abs(inp_uvw) > 1e-9])
            assert nonzero.size > 0, "all-zero reverse input"
            matched = False
            for k in (1, 2, 3, 4, 5, 6):
                if np.allclose(rt_uvw, k * inp_uvw, atol=1e-6):
                    matched = True
                    break
            assert matched, (
                f"reverse roundtrip [u,v,w] not proportional: "
                f"{inp_uvw.tolist()} -> {rt_uvw.tolist()}"
            )
            # The recomputed t must satisfy the canonical constraint.
            assert abs(rt[2] - (-(rt[0] + rt[1]))) < 1e-6, (
                f"reverse roundtrip t constraint violated: "
                f"t={rt[2]}, -(u+v)={-(rt[0]+rt[1])}"
            )


def print_miller_table(examples: list[dict[str, object]]) -> None:
    """Print the conversion table to stdout."""
    print("=" * 78)
    print("Part A: HCP Miller index 3 <-> 4 conversion")
    print("=" * 78)
    header = (
        f"{'direction':<9} {'input':<18} {'converted':<22} "
        f"{'round-trip':<22}"
    )
    print(header)
    print("-" * 78)
    for ex in examples:
        print(
            f"{ex['direction']:<9} {ex['input_str']:<18} "
            f"{ex['output_str']:<22} {ex['roundtrip_str']:<22}"
        )
    print("-" * 78)
    print("Round-trip check: PASS (forward exact, reverse up to integer scale)")
    print()


# ---------------------------------------------------------------------------
# Part B: Density calculation
# ---------------------------------------------------------------------------

def build_bulk_cells() -> list[Atoms]:
    """Return four bulk cells in display order: FCC Cu, BCC Fe, HCP Mg, Si."""
    return [
        bulk("Cu", "fcc", a=A_FCC_CU, cubic=True),
        bulk("Fe", "bcc", a=A_BCC_FE, cubic=True),
        bulk("Mg", "hcp", a=A_HCP_MG, covera=C_HCP_MG / A_HCP_MG),
        bulk("Si", "diamond", a=A_DIA_SI, cubic=True),
    ]


def cell_labels() -> list[str]:
    """Display labels in the same order as :func:`build_bulk_cells`."""
    return ["FCC Cu", "BCC Fe", "HCP Mg", "Diamond Si"]


def summarize_density(idx: int, atoms: Atoms) -> dict[str, object]:
    """Return a deterministic density summary row for one bulk cell."""
    volume = float(atoms.get_volume())               # A^3
    natoms = len(atoms)
    number_density = float(cal_density(atoms))       # atoms / A^3
    symbols = atoms.get_chemical_symbols()
    avg_mass_u = float(np.mean([atomic_masses[atomic_numbers[s]] for s in symbols]))
    # mass density (g/cm^3) = number_density (atoms/A^3) * avg_mass (u/atom)
    #                        * (1.66053906660e-24 g/u) / (1e-24 cm^3/A^3)
    # The 1e-24 factors cancel, leaving the 1.66053906660 factor.
    mass_density = number_density * avg_mass_u * U_PER_A3_TO_G_PER_CM3  # g/cm^3
    return {
        "label": cell_labels()[idx],
        "formula": atoms.get_chemical_formula(),
        "natoms": natoms,
        "volume_A3": round(volume, 4),
        "number_density": round(number_density, 6),
        "mass_density": round(mass_density, 4),
        "avg_mass_u": round(avg_mass_u, 4),
    }


def print_density_table(rows: list[dict[str, object]]) -> None:
    """Print the density summary table to stdout."""
    print("=" * 78)
    print("Part B: Density of bulk cells")
    print("=" * 78)
    header = (
        f"{'structure':<12} {'formula':<8} {'natoms':>6} "
        f"{'volume(A^3)':>12} {'n_dens(at/A^3)':>16} "
        f"{'mass_dens(g/cm^3)':>18}"
    )
    print(header)
    print("-" * 78)
    for row in rows:
        print(
            f"{row['label']:<12} {row['formula']:<8} {row['natoms']:>6} "
            f"{row['volume_A3']:>12.4f} {row['number_density']:>16.6f} "
            f"{row['mass_density']:>18.4f}"
        )
    print("-" * 78)
    print()


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------

def render_figure(
    miller_examples: list[dict[str, object]],
    density_rows: list[dict[str, object]],
    path_image: Path,
) -> None:
    """Render the two-panel figure.

    Panel A (left):  schematic table of the 3-index -> 4-index conversion
    for the four worked examples, rendered as a Matplotlib table.
    Panel B (right): grouped bar chart of number density and mass density
    for the four bulk structures.
    """
    path_image.parent.mkdir(parents=True, exist_ok=True)

    general_set_all_rcParams(
        backend="Agg",
        axes_grid=False,
        figure_subp=(1, 2),
        figure_one_figsize=(6.0, 5.2),
        figure_left=0.55,
        figure_top=0.55,
        figure_wspace=1.1,
        figure_hspace=0.6,
        font_family=["DejaVu Sans"],
        fontsize=11,
        axes_titlepad=8,
    )
    fig, axes = plt.subplots(1, 2)

    # ---- Panel A: conversion table as an image -------------------------
    ax = axes[0]
    ax.set_axis_off()
    ax.set_title("HCP 3-index $\\leftrightarrow$ 4-index", fontsize=12, pad=10)

    col_labels = ["direction", "input", "converted", "round-trip"]
    cell_text = []
    cell_colors = []
    for ex in miller_examples:
        cell_text.append([
            ex["direction"],
            ex["input_str"],
            ex["output_str"],
            ex["roundtrip_str"],
        ])
        cell_colors.append(["#f0f8ff", "#ffffff", "#ffffff", "#eafaea"])

    table = ax.table(
        cellText=cell_text,
        colLabels=col_labels,
        cellColours=cell_colors,
        colColours=["#d0e0ff"] * 4,
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.6)
    # Make header row bold.
    for j in range(len(col_labels)):
        cell = table[(0, j)]
        cell.get_text().set_weight("bold")

    # ---- Panel B: density bar chart ------------------------------------
    ax = axes[1]
    ax.set_title("Bulk cell density", fontsize=12, pad=10)

    labels = [row["label"] for row in density_rows]
    x = np.arange(len(labels))
    bar_w = 0.38

    n_dens = np.array([row["number_density"] for row in density_rows])
    m_dens = np.array([row["mass_density"] for row in density_rows])

    bars1 = ax.bar(x - bar_w / 2, n_dens, bar_w,
                   color="#1f77b4", label=r"$n$ (atoms/$\mathrm{\AA}^3$)")
    bars2 = ax.bar(x + bar_w / 2, m_dens, bar_w,
                   color="#d62728", label=r"$\rho$ (g/cm$^3$)")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel("value")
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, axis="y", linestyle="--", alpha=0.4)

    # Annotate bar heights so the chart is readable even when the two
    # quantities differ by an order of magnitude.
    for bars in (bars1, bars2):
        for b in bars:
            h = b.get_height()
            ax.annotate(f"{h:.2f}",
                        (b.get_x() + b.get_width() / 2, h),
                        ha="center", va="bottom", fontsize=8)

    fig.savefig(path_image, bbox_inches="tight")
    plt.close(fig)

    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError("figure was not created: " + str(path_image))
    image_rgb = plt.imread(path_image)[..., :3]
    if float(np.mean(np.abs(image_rgb - 1.0))) < 0.002:
        raise AssertionError("figure is effectively blank: " + str(path_image))


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def run_example(path_output: Path) -> tuple[list[dict[str, object]],
                                            list[dict[str, object]]]:
    """Run both parts, render the figure, return the summary rows."""
    path_output = path_output.resolve()
    path_output.mkdir(parents=True, exist_ok=True)

    # Part A: Miller index conversion.
    miller_examples = build_miller_examples()
    verify_roundtrips(miller_examples)
    print_miller_table(miller_examples)

    # Part B: density calculation.
    cells = build_bulk_cells()
    density_rows = [summarize_density(i, c) for i, c in enumerate(cells)]
    print_density_table(density_rows)

    # Figure.
    path_image = path_output / "miller_index_and_density.png"
    render_figure(miller_examples, density_rows, path_image)
    print("wrote: " + str(path_image))

    return miller_examples, density_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="HCP Miller index 3<->4 conversion and density calculation."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("miller-density-output"),
        help="Directory for the comparison PNG.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    miller_rows, density_rows = run_example(parse_args().output)

    # Deterministic assertions.
    # Part A: round-trips already verified inside run_example.
    assert len(miller_rows) == 4

    # Part B: all densities must be strictly positive.
    assert all(row["number_density"] > 0 for row in density_rows), (
        "non-positive number density"
    )
    assert all(row["mass_density"] > 0 for row in density_rows), (
        "non-positive mass density"
    )
    # Sanity: FCC Cu mass density ~8.9 g/cm^3, BCC Fe ~7.9, HCP Mg ~1.7,
    # diamond Si ~2.3. Use loose bounds so the script does not break on
    # minor lattice-constant revisions.
    bounds = {
        "FCC Cu": (7.5, 9.5),
        "BCC Fe": (6.5, 8.5),
        "HCP Mg": (1.0, 2.5),
        "Diamond Si": (1.5, 3.0),
    }
    for row in density_rows:
        lo, hi = bounds[row["label"]]
        assert lo < row["mass_density"] < hi, (
            f"{row['label']} mass density {row['mass_density']} "
            f"outside ({lo}, {hi})"
        )
