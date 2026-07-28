#!/usr/bin/env python3
"""Recalculate one tracked Au(111) surface-energy fixture."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from mymetal.build.film.extrfilm import cal_area
from mymetal.calculate.calenergy.surfenergy import cal_surface_energy
from mymetal.io.vasp import my_read_vasp


PATH_REPOSITORY = Path(__file__).resolve().parents[2]
PATH_DEFAULT_CASE = (
    PATH_REPOSITORY
    / "mymetal"
    / "example"
    / "test-surface-energy"
    / "fcc"
    / "1.000-2.8485"
)
BULK_ENERGY_EV = -47.04775262
SLAB_ENERGY_EV = -46.24580704
SURFACE_FACTOR = 2


def calculate_surface_energy(path_case: Path) -> tuple[float, float, float]:
    """Return area, eV/Å², and J/m² for the selected repository fixture."""
    path_bulk = path_case / "bulk" / "CONTCAR"
    path_slab = path_case / "full_relaxed_surface_111" / "CONTCAR"
    if not path_bulk.is_file() or not path_slab.is_file():
        raise FileNotFoundError("fixture CONTCAR is missing under " + str(path_case))

    bulk, _ = my_read_vasp(path_bulk)
    slab, _ = my_read_vasp(path_slab)
    area = float(cal_area(slab))
    gamma_manual = (
        SLAB_ENERGY_EV
        - BULK_ENERGY_EV / len(bulk) * len(slab)
    ) / (SURFACE_FACTOR * area)
    gamma_ev_a2 = float(cal_surface_energy(
        bulk_energy=BULK_ENERGY_EV,
        bulk_atoms_number=len(bulk),
        relaxed_surface_energy=SLAB_ENERGY_EV,
        surface_atoms_number=len(slab),
        area=area,
        energy_unit="eV",
        factor=SURFACE_FACTOR,
    ))
    gamma_j_m2 = float(cal_surface_energy(
        bulk_energy=BULK_ENERGY_EV,
        bulk_atoms_number=len(bulk),
        relaxed_surface_energy=SLAB_ENERGY_EV,
        surface_atoms_number=len(slab),
        area=area,
        energy_unit="J",
        factor=SURFACE_FACTOR,
    ))

    assert len(bulk) == 12
    assert len(slab) == 12
    assert np.isclose(gamma_ev_a2, gamma_manual, atol=1e-12)
    assert np.isclose(gamma_ev_a2, 0.05706263510343279, atol=1e-12)
    return area, gamma_ev_a2, gamma_j_m2


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Recalculate the tracked Au(111) surface-energy example.")
    parser.add_argument(
        "-path_case",
        type=Path,
        default=PATH_DEFAULT_CASE,
        help="Case directory containing bulk/ and full_relaxed_surface_111/.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    path_case = parse_args().path_case.resolve()
    area_a2, surface_energy_ev_a2, surface_energy_j_m2 = (
        calculate_surface_energy(path_case)
    )
    print("bulk atoms: 12")
    print("slab atoms: 12")
    print(f"area_A2: {area_a2:.10f}")
    print("factor: " + str(SURFACE_FACTOR))
    print(f"manual_eV_A2: {surface_energy_ev_a2:.10f}")
    print(f"api_eV_A2: {surface_energy_ev_a2:.10f}")
    print(f"api_J_m2: {surface_energy_j_m2:.10f}")
    print("check: ok")
