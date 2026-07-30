#!/usr/bin/env python3
"""Move, fix, and selectively constrain atoms in an ASE structure.

This tutorial is VASP-free. It demonstrates three universal helpers that
operate on :class:`ase.Atoms` objects:

* :func:`mymetal.universal.atom.moveatom.move_atoms` -- translate a
  selected subset of atoms (by symbol, by Cartesian position range, or
  by all atoms) in either fractional or Cartesian coordinates.
* :func:`mymetal.universal.atom.fixatom.fixatoms` -- attach an
  :class:`ase.constraints.FixAtoms` constraint from a boolean mask or an
  index list, so VASP/SLURM runs keep the chosen atoms frozen.
* :func:`mymetal.io.vasp.my_write_vasp` with selective dynamics -- write
  the constraint back into a ``POSCAR`` so the next job reads it.

The script builds an Au(111) slab, freezes its bottom half, translates
the top layer by a fractional shift, and writes a ``POSCAR`` whose
``Selective dynamics`` block matches the constraint. It then renders a
two-panel before/after figure with frozen atoms marked as faded disks.
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

from mymetal.build.film.stretch import generate_film
from mymetal.io.vasp import my_read_vasp, my_write_vasp
from mymetal.universal.atom.fixatom import fixatoms
from mymetal.universal.atom.moveatom import move_atoms
from mymetal.universal.plot.general import general_set_all_rcParams


A_FCC_AU = 4.08
NUM_LAYERS = 6
VACUUM = 15.0


def build_slab() -> Atoms:
    """Return a 6-layer Au(111) slab used as the manipulation reference."""
    return generate_film(
        symbols="Au",
        structure="fcc",
        num_layers=NUM_LAYERS,
        my_vacuum=VACUUM,
        slice_plane=(1, 1, 1),
        a_fcc=A_FCC_AU,
    )


def fix_bottom_half(atoms: Atoms) -> Atoms:
    """Return a copy of ``atoms`` with the bottom half frozen."""
    z_min = float(atoms.positions[:, 2].min())
    z_max = float(atoms.positions[:, 2].max())
    z_mid = 0.5 * (z_min + z_max)
    mask = [atom.position[2] < z_mid for atom in atoms]
    constrained = fixatoms(atoms, mask=mask)
    return constrained


def shift_top_layer(atoms: Atoms, shift: list[float] = None) -> Atoms:
    """Return a copy of ``atoms`` whose top layer is shifted by ``shift``.

    ``shift`` is given in fractional coordinates. The default
    ``[0.25, 0.0, 0.0]`` moves the top-layer Au atoms by one quarter of
    the first cell vector, which corresponds to a bridge-site shift on
    Au(111).
    """
    if shift is None:
        shift = [0.25, 0.0, 0.0]
    z_max = float(atoms.positions[:, 2].max())
    # Use a 0.1 A window around the topmost atom to define "top layer".
    return move_atoms(
        atoms,
        translate_matrix=np.asarray(shift),
        if_scale_position=True,
        atom_type=None,
        position_range=((-np.inf, np.inf), (-np.inf, np.inf),
                        (z_max - 0.1, z_max + 0.1)),
    )


def render_figure(
    slab: Atoms,
    slab_constrained: Atoms,
    slab_shifted: Atoms,
    mask: list[bool],
    path_image: Path,
) -> None:
    """Render a 1x3 before/during/after figure."""
    path_image.parent.mkdir(parents=True, exist_ok=True)

    general_set_all_rcParams(
        backend="Agg",
        axes_grid=False,
        figure_subp=(1, 3),
        figure_one_figsize=(4.4, 4.6),
        figure_left=0.55,
        figure_top=0.50,
        figure_wspace=1.0,
        figure_hspace=0.6,
        font_family=["DejaVu Sans"],
        fontsize=11,
        axes_titlepad=6,
    )
    fig, axes = plt.subplots(1, 3)

    titles = [
        "Reference slab\n(no constraints)",
        "Bottom half frozen\n(FixAtoms mask)",
        "Top layer shifted\n(bridge-site move)",
    ]
    panels = [slab, slab_constrained, slab_shifted]
    for ax, atoms, title in zip(axes, panels, titles):
        plot_atoms(
            atoms.repeat((3, 1, 1)),
            ax,
            rotation="90x,0y,0z",
            show_unit_cell=0,
            radii=0.7,
        )
        ax.set_title(title)
        ax.set_axis_off()

    # Overlay: mark the frozen atoms with a faint translucent rectangle.
    # We compute the bounding box of the frozen atoms in the side view
    # (after the 3x repeat), then draw it on the second panel.
    ax = axes[1]
    repeated = slab_constrained.repeat((3, 1, 1))
    mask_repeated = list(mask) * 3
    frozen_pos = repeated.positions[np.asarray(mask_repeated)]
    if len(frozen_pos) > 0:
        x_min, x_max = (frozen_pos[:, 0].min(),
                        frozen_pos[:, 0].max())
        z_min, z_max = (frozen_pos[:, 2].min(),
                        frozen_pos[:, 2].max())
        ax.add_patch(plt.Rectangle(
            (x_min - 0.3, z_min - 0.3),
            x_max - x_min + 0.6,
            z_max - z_min + 0.6,
            fill=True,
            facecolor="#d62728",
            alpha=0.15,
            edgecolor="#d62728",
            linewidth=1.2,
            linestyle="--",
        ))

    # Overlay: draw an arrow on the third panel showing the shift
    # direction in the side view.
    ax = axes[2]
    repeated = slab_shifted.repeat((3, 1, 1))
    z_max = float(repeated.positions[:, 2].max())
    top_atom = repeated.positions[:, 2].argmax()
    x0 = float(repeated.positions[top_atom, 0])
    # Translate the fractional shift [0.25, 0, 0] into Cartesian.
    shift_cart = 0.25 * np.asarray(slab_shifted.cell.array[0])
    dx = float(shift_cart[0])
    ax.annotate(
        "",
        xy=(x0 + dx, z_max),
        xytext=(x0, z_max),
        arrowprops=dict(arrowstyle="->", color="#d62728", linewidth=2.0),
    )

    fig.savefig(path_image, bbox_inches="tight")
    plt.close(fig)

    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError("manipulation image was not created: " + str(path_image))
    image_rgb = plt.imread(path_image)[..., :3]
    if float(np.mean(np.abs(image_rgb - 1.0))) < 0.002:
        raise AssertionError("manipulation image is effectively blank: " + str(path_image))


def run_example(path_output: Path):
    """Build the slab, fix and shift atoms, write POSCAR, render figure."""
    path_output = path_output.resolve()
    path_output.mkdir(parents=True, exist_ok=True)

    slab = build_slab()
    z_min = float(slab.positions[:, 2].min())
    z_max = float(slab.positions[:, 2].max())
    z_mid = 0.5 * (z_min + z_max)
    mask = [atom.position[2] < z_mid for atom in slab]

    slab_constrained = fix_bottom_half(slab)
    slab_shifted = shift_top_layer(slab_constrained)

    path_poscar = path_output / "POSCAR"
    my_write_vasp(path_poscar, slab_shifted, label="Au(111) selective dynamics")
    slab_read, _ = my_read_vasp(path_poscar)

    path_image = path_output / "atom_manipulation.png"
    render_figure(slab, slab_constrained, slab_shifted, mask, path_image)

    print(f"slab formula: {slab.get_chemical_formula()}")
    print(f"slab atoms: {len(slab)}")
    print(f"frozen atoms: {sum(mask)}")
    print(f"moving atoms: {len(slab) - sum(mask)}")
    print(f"top layer shift (fractional): [0.25, 0.0, 0.0]")
    print(f"POSCAR round-trip: ok")
    print(f"  formula: {slab_read.get_chemical_formula()}")
    print(f"  atoms:   {len(slab_read)}")
    print(f"  constraint: {type(slab_read.constraints[0]).__name__}")
    print("wrote: " + str(path_poscar))
    print("wrote: " + str(path_image))
    return {
        "slab": slab,
        "slab_constrained": slab_constrained,
        "slab_shifted": slab_shifted,
        "mask": mask,
        "slab_read": slab_read,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Move, fix, and write selective-dynamics POSCAR for Au(111).")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("atom-manipulation-output"),
        help="Directory for POSCAR and the manipulation PNG.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    result = run_example(parse_args().output)
    slab = result["slab"]
    slab_shifted = result["slab_shifted"]
    slab_read = result["slab_read"]
    mask = result["mask"]

    # The slab must have an even number of atoms so half/half works.
    assert len(slab) % 2 == 0
    # The frozen mask must freeze exactly half the slab.
    assert sum(mask) == len(slab) // 2
    # The constraint must survive the POSCAR round trip.
    assert len(slab_read.constraints) >= 1
    # The shift must move only the top-layer atoms; the bottom layer
    # must be unchanged in Cartesian coordinates.
    z_max = float(slab.positions[:, 2].max())
    top_idx = [i for i, atom in enumerate(slab) if atom.position[2] > z_max - 0.1]
    bot_idx = [i for i, atom in enumerate(slab) if atom.position[2] < z_max - 0.1]
    for i in bot_idx:
        assert np.allclose(slab.positions[i], slab_shifted.positions[i], atol=1e-10)
    # The top layer must move by the chosen shift in fractional coords.
    scaled_before = slab.get_scaled_positions()
    scaled_after = slab_shifted.get_scaled_positions()
    for i in top_idx:
        diff = scaled_after[i] - scaled_before[i]
        assert abs(diff[0] - 0.25) < 1e-6
        assert abs(diff[1]) < 1e-6
        assert abs(diff[2]) < 1e-6
