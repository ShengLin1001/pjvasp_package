#!/usr/bin/env python3
"""Build, stretch, and visualize a Cu(111) slab along multiple directions.

This tutorial is VASP-free. It demonstrates how
:func:`mymetal.build.film.stretch.stretch_list_along_direction_to_cell`
applies a list of stretch factors along chosen cell directions, and how the
strain energy of each mode would look for an elastic solid with Cu-like
elastic constants.

The script:

1. Builds a 4-layer Cu(111) slab.
2. Stretches it along ``x`` (in-plane uniaxial), ``z`` (out-of-plane
   uniaxial) and ``xy`` (in-plane biaxial) using a symmetric strain list.
3. Synthesizes a strain-energy curve for each direction from cubic elastic
   constants (no VASP involved); the constants are clearly labeled as
   literature-like values for teaching purposes.
4. Renders a two-panel figure: the deformed cells at the largest compressive
   strain, and the three strain-energy curves.
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

from mymetal.build.film.stretch import (
    generate_film,
    stretch_list_along_direction_to_cell,
)
from mymetal.universal.plot.general import general_set_all_rcParams


# Cu conventional lattice constant (Å). Reused for consistency with the
# fcc_surfaces tutorial.
A_FCC_CU = 3.61
NUM_LAYERS = 4
VACUUM = 12.0

# Symmetric stretch factors. 0.03 means a 3% tensile strain along the chosen
# direction; the same magnitude on the negative side gives 3% compression.
STRETCH_FACTORS = [0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03]
DIRECTIONS = ["x", "z", "xy"]
DIRECTION_LABELS = {
    "x": "in-plane uniaxial (x)",
    "z": "out-of-plane uniaxial (z)",
    "xy": "in-plane biaxial (xy)",
}

# Cu-like cubic elastic constants (GPa). Used only to synthesize a physically
# reasonable strain-energy curve for teaching; the script does not claim these
# are a published benchmark.
C11_GPA = 168.0
C12_GPA = 121.0
EV_PER_A3_TO_GPA = 160.21766


def build_slab() -> Atoms:
    """Return the 4-layer Cu(111) slab used as the reference configuration."""
    return generate_film(
        symbols="Cu",
        structure="fcc",
        num_layers=NUM_LAYERS,
        my_vacuum=VACUUM,
        slice_plane=(1, 1, 1),
        a_fcc=A_FCC_CU,
    )


def stretch_along(slab: Atoms, direction: str) -> list[Atoms]:
    """Return one stretched Atoms per entry in :data:`STRETCH_FACTORS`."""
    return stretch_list_along_direction_to_cell(
        slab,
        stretch_factor_list=STRETCH_FACTORS,
        stretch_direction_list=[direction],
    )


def synthesize_energy(
    atoms_list: list[Atoms],
    reference: Atoms,
    direction: str,
    seed: int = 20260729,
) -> np.ndarray:
    """Return a synthetic strain-energy curve in eV per atom.

    The model treats the slab as a cubic elastic solid whose energy density
    follows the Voigt-reduced form ``u = 1/2 * C_IJ * eta_I * eta_J``. Only
    the entries that the chosen direction activates are kept, and a small
    Gaussian noise is added so the curves look like real DFT data. The model
    is not a substitute for relaxation: it intentionally ignores Poisson
    coupling so the three directions produce visibly different slopes.

    Args:
        atoms_list (list): Stretched Atoms produced by :func:`stretch_along`.
        reference (Atoms): The un-stretched slab, used for the reference
            volume and lattice constants.
        direction (str): One of ``x``, ``z`` or ``xy``.
        seed (int): RNG seed for the noise.

    Returns:
        numpy.ndarray: Energy per atom in eV for each entry in ``atoms_list``.
    """
    rng = np.random.default_rng(seed)
    v0_per_atom = reference.get_volume() / len(reference)
    c11 = C11_GPA / EV_PER_A3_TO_GPA
    c12 = C12_GPA / EV_PER_A3_TO_GPA

    strains = np.array(STRETCH_FACTORS) - 1.0
    if direction == "x":
        # eta1 = strain, eta2 = eta3 = 0 -> u = 1/2 * C11 * eta1^2
        energy_density = 0.5 * c11 * strains ** 2
    elif direction == "z":
        # eta3 = strain -> u = 1/2 * C11 * eta3^2 (here c33 == c11 for cubic)
        energy_density = 0.5 * c11 * strains ** 2
    elif direction == "xy":
        # eta1 = eta2 = strain -> u = 1/2 * (2*C11 + 2*C12) * strain^2
        #                          = (C11 + C12) * strain^2
        energy_density = (c11 + c12) * strains ** 2
    else:
        raise ValueError("direction must be one of 'x', 'z', 'xy'")

    energy_per_atom = energy_density * v0_per_atom
    noise = rng.normal(0.0, 0.0005, size=energy_per_atom.shape)
    return energy_per_atom + noise


def render_figure(
    slab: Atoms,
    stretched_per_direction: dict[str, list[Atoms]],
    energies_per_direction: dict[str, np.ndarray],
    path_image: Path,
) -> None:
    """Render the deformed cells and the strain-energy curves."""
    path_image.parent.mkdir(parents=True, exist_ok=True)

    general_set_all_rcParams(
        backend="Agg",
        axes_grid=True,
        figure_subp=(2, 3),
        figure_one_figsize=(4.4, 4.4),
        figure_left=0.65,
        figure_top=0.55,
        figure_wspace=1.0,
        figure_hspace=0.6,
        font_family=["DejaVu Sans"],
        fontsize=11,
        axes_titlepad=6,
    )
    fig, axes = plt.subplots(2, 3)

    # Top row: side view of the most-compressed cell (strain = -3%) for each
    # direction, plus the reference cell outline so the deformation is visible.
    for col, direction in enumerate(DIRECTIONS):
        ax = axes[0, col]
        atoms_list = stretched_per_direction[direction]
        idx_compress = int(np.argmin(STRETCH_FACTORS))
        deformed = atoms_list[idx_compress]
        plot_atoms(
            deformed.repeat((3, 1, 1)),
            ax,
            rotation="90x,0y,0z",
            show_unit_cell=0,
            radii=0.7,
        )
        ax.set_title(DIRECTION_LABELS[direction] + f"\nstrain = {STRETCH_FACTORS[idx_compress]-1:+.2f}")
        ax.set_axis_off()

    # Bottom row: strain-energy curves.
    strains = np.array(STRETCH_FACTORS) - 1.0
    for col, direction in enumerate(DIRECTIONS):
        ax = axes[1, col]
        energies = energies_per_direction[direction]
        ax.plot(strains * 100, energies, "o-", color="#1f77b4")
        ax.axhline(0.0, color="gray", linestyle=":")
        ax.axvline(0.0, color="gray", linestyle=":")
        ax.set_xlabel("strain (%)")
        ax.set_ylabel("E - E0 (eV/atom)")
        ax.set_title(DIRECTION_LABELS[direction])
        ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))

    fig.savefig(path_image, bbox_inches="tight")
    plt.close(fig)

    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError("stretch image was not created: " + str(path_image))


def run_example(path_output: Path) -> dict[str, np.ndarray]:
    """Build, stretch, synthesize energies, render, and return the energies."""
    path_output = path_output.resolve()
    path_output.mkdir(parents=True, exist_ok=True)

    slab = build_slab()
    stretched_per_direction: dict[str, list[Atoms]] = {}
    energies_per_direction: dict[str, np.ndarray] = {}
    for direction in DIRECTIONS:
        atoms_list = stretch_along(slab, direction)
        stretched_per_direction[direction] = atoms_list
        energies_per_direction[direction] = synthesize_energy(
            atoms_list, slab, direction, seed=20260729 + DIRECTIONS.index(direction),
        )

    path_image = path_output / "biaxial_stretch.png"
    render_figure(slab, stretched_per_direction, energies_per_direction, path_image)

    print("direction  strains(%)  E-E0 (meV/atom)")
    strains_pct = [(f - 1.0) * 100 for f in STRETCH_FACTORS]
    for direction in DIRECTIONS:
        energies = energies_per_direction[direction] * 1000.0
        energy_str = "  ".join(f"{e:7.3f}" for e in energies)
        print(f"{direction:<9}  {strains_pct}  {energy_str}")
    print("wrote: " + str(path_image))
    return energies_per_direction


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build, stretch and visualize a Cu(111) slab along x/z/xy.")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("biaxial-stretch-output"),
        help="Directory for the stretch PNG.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    energies_per_direction = run_example(parse_args().output)
    # Every direction must produce one energy per stretch factor.
    for direction in DIRECTIONS:
        assert len(energies_per_direction[direction]) == len(STRETCH_FACTORS)
    # The reference (strain = 0) must be the lowest energy within each curve
    # up to the small synthetic noise.
    idx_zero = STRETCH_FACTORS.index(1.0)
    for direction in DIRECTIONS:
        e = energies_per_direction[direction]
        assert e[idx_zero] < e[0] + 1e-3, (direction, e)
        assert e[idx_zero] < e[-1] + 1e-3, (direction, e)
    # The biaxial direction must be stiffer than the uniaxial one at the
    # largest strain (C11 + C12 > C11).
    assert energies_per_direction["xy"][-1] > energies_per_direction["x"][-1]
