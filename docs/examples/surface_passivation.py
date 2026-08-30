"""Surface passivation tutorial example.

Builds an Au(111) slab, then adds hydrogen atoms to passivate the top and
bottom surfaces using ``add_hydro_atoms``.  Also demonstrates the dangling-bond
detection utilities (``get_neighbors``, ``find_matching_atom_in_bulk``,
``find_unsaturated_neighbors``) on a small Si bulk/slab pair.

This example does NOT run VASP, does NOT need POTCAR, and does NOT call sbatch.
"""

import argparse
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.build import bulk, surface
from ase.io import write

from mymetal.build.film.stretch import generate_film
from mymetal.build.film.hydroxyl import (
    add_hydro_atoms,
    get_neighbors,
    find_matching_atom_in_bulk,
    find_unsaturated_neighbors,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Surface passivation example")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("docs/_build/example-passivation"),
        help="Output directory",
    )
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)

    # ---- Part 1: Au(111) slab + H passivation ----
    slab = generate_film(
        symbols="Au",
        structure="fcc",
        num_layers=6,
        my_vacuum=15.0,
        slice_plane=(1, 1, 1),
        a_fcc=4.08,
    )
    z_max = float(slab.get_positions()[:, 2].max())
    z_min = float(slab.get_positions()[:, 2].min())
    assert slab.get_chemical_formula() == "Au6"
    print(f"Au(111) slab: {slab.get_chemical_formula()}, "
          f"z range {z_min:.3f}–{z_max:.3f} Å")

    # Passivate top surface: add H 1.8 Å above each Au surface atom
    passivated = add_hydro_atoms(
        slab.copy(),
        add_symbol="H",
        added_symbol="Au",
        surf_range=[z_max - 2.0, z_max + 5.0],
        shift_distance=[0.0, 0.0, 1.8],
        surf_direction=2,
    )
    assert "H" in passivated.get_chemical_symbols()
    assert len(passivated) == 7  # 6 Au + 1 H
    print(f"Top passivated: {passivated.get_chemical_formula()}")

    # Passivate bottom surface: add H 1.8 Å below each Au surface atom
    passivated = add_hydro_atoms(
        passivated.copy(),
        add_symbol="H",
        added_symbol="Au",
        surf_range=[z_min - 5.0, z_min + 2.0],
        shift_distance=[0.0, 0.0, -1.8],
        surf_direction=2,
    )
    assert len(passivated) == 8  # 6 Au + 2 H
    print(f"Both surfaces passivated: {passivated.get_chemical_formula()}")

    write(out / "POSCAR_passivated", passivated, format="vasp", direct=True)
    print(f"Wrote {out / 'POSCAR_passivated'}")

    # ---- Part 2: Dangling-bond detection on Si ----
    si_bulk = bulk("Si", "diamond", a=5.43, cubic=True)
    # Use one repeat so slab atoms stay within bulk coordinate range
    si_slab = si_bulk.repeat((1, 1, 1))
    print(f"\nSi bulk: {si_bulk.get_chemical_formula()}, "
          f"{len(si_bulk)} atoms")
    print(f"Si slab (1 repeat): {si_slab.get_chemical_formula()}, "
          f"{len(si_slab)} atoms")

    # Demonstrate get_neighbors on first Si atom
    positions = si_bulk.get_positions()
    cell = si_bulk.get_cell()
    pbc = [True, True, True]
    neighbors, offsets, distances = get_neighbors(
        0, positions, np.array(cell), pbc, cutoff=2.5,
    )
    print(f"Si atom 0 has {len(neighbors)} neighbors within 2.5 Å")
    assert len(neighbors) == 4  # diamond Si: 4-fold coordination
    print(f"  neighbor distances (Å): {[f'{d:.3f}' for d in distances]}")

    # Demonstrate find_matching_atom_in_bulk
    idx_match = find_matching_atom_in_bulk(
        si_bulk.positions[0], si_bulk.positions, tolerance=1e-5,
    )
    assert idx_match == 0
    print(f"Si bulk atom 0 matches bulk index {idx_match}")

    # Demonstrate find_unsaturated_neighbors
    # Simulate a missing neighbor by removing one offset
    offset_bulk = offsets  # all 4 neighbors in bulk
    offset_slab = offsets[:3]  # only 3 in slab (one dangling bond)
    unsaturated = find_unsaturated_neighbors(
        offset_bulk, offset_slab, tolerance=1e-3,
    )
    assert len(unsaturated) == 1
    print(f"Found {len(unsaturated)} unsaturated (dangling) bond(s)")
    print(f"  missing offset: {unsaturated[0]}")

    print("\nAll assertions passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
