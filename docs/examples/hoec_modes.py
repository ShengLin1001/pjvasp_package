"""Higher-order elastic constants (HOEC) example.

Demonstrates the ``mymetal.calculate.calmechanics.hoec`` module:
runs the built-in self-test, inspects the cubic and hexagonal mode tables,
and computes deformation gradients for representative strain modes.

This example does NOT run VASP, does NOT need POTCAR, and does NOT call sbatch.
It only depends on ``numpy`` and the Python standard library.
"""

import argparse
from pathlib import Path

import numpy as np
from ase.build import bulk

from mymetal.calculate.calmechanics.hoec import (
    get_model,
    get_hoec_modes,
    get_strain_list,
    get_deformation_gradient,
    check_symmetry,
    selftest_hoec,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="HOEC example")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("docs/_build/example-hoec"),
        help="Output directory",
    )
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)

    # ---- Part 1: Self-test (validates against Wang-Li Table I) ----
    print("Running selftest_hoec...")
    selftest_hoec()
    print("self-test PASSED\n")

    # ---- Part 2: Cubic mode table ----
    model_cubic = get_model("cubic")
    modes_cubic = get_hoec_modes("cubic")
    names2 = model_cubic.names(2)
    names3 = model_cubic.names(3)
    names4 = model_cubic.names(4)
    print(f"Cubic model: {len(names2)} SOEC, "
          f"{len(names3)} TOEC, {len(names4)} FOEC")
    print(f"Cubic independent constants:")
    print(f"  SOEC: {names2}")
    print(f"  TOEC: {names3}")
    print(f"  FOEC: {names4}")
    print(f"Cubic mode table ({len(modes_cubic)} modes):")
    for name, d in modes_cubic.items():
        print(f"  {name}: d = {d}")
    assert len(modes_cubic) == 11  # A..K

    # ---- Part 3: Hexagonal mode table ----
    model_hex = get_model("hex")
    modes_hex = get_hoec_modes("hex")
    print(f"\nHex model: {len(model_hex.names(2))} SOEC, "
          f"{len(model_hex.names(3))} TOEC, {len(model_hex.names(4))} FOEC")
    print(f"Hex mode table ({len(modes_hex)} modes):")
    for name, d in list(modes_hex.items())[:5]:
        print(f"  {name}: d = {d}")
    assert len(modes_hex) == 23  # M01..M23

    # ---- Part 4: Strain list and deformation gradient ----
    xi_list = get_strain_list(emax=0.06, de=0.02)
    print(f"\nStrain list (emax=0.06, de=0.02): {xi_list}")
    assert xi_list == [-0.06, -0.04, -0.02, 0.0, 0.02, 0.04, 0.06]

    # Mode A = uniaxial x strain: d = (1, 0, 0, 0, 0, 0)
    d_A = modes_cubic["A"]
    F = get_deformation_gradient(d_A, xi=0.02)
    print(f"\nMode A (uniaxial x), xi=0.02:")
    print(f"  d = {d_A}")
    print(f"  F =")
    for row in F:
        print(f"    [{row[0]:+.8f}, {row[1]:+.8f}, {row[2]:+.8f}]")

    # Verify: eta = xi * diag(1,0,0) => E = 1/2(F^T F - I) should equal eta
    E = 0.5 * (F.T @ F - np.eye(3))
    expected_eta = np.diag([0.02, 0, 0])
    print(f"  Lagrangian strain E =")
    for row in E:
        print(f"    [{row[0]:+.8f}, {row[1]:+.8f}, {row[2]:+.8f}]")
    assert np.allclose(E, expected_eta, atol=1e-10), \
        f"E mismatch: {E} vs {expected_eta}"
    print(f"  E matches xi*d within 1e-10 ✅")

    # ---- Part 5: Symmetry detection ----
    cu = bulk("Cu", "fcc", a=3.61, cubic=True)
    sym = check_symmetry(symmetry="auto", cell=cu.cell)
    print(f"\nDetected symmetry for FCC Cu: {sym}")
    assert sym == "cubic"

    # Mode F = pure shear: d = (0, 0, 0, 2, 2, 2)
    d_F = modes_cubic["F"]
    F_shear = get_deformation_gradient(d_F, xi=0.03)
    E_shear = 0.5 * (F_shear.T @ F_shear - np.eye(3))
    print(f"\nMode F (pure shear yz/xz/xy), xi=0.03:")
    print(f"  d = {d_F}")
    # eta4=2*eps_yz, eta5=2*eps_xz, eta6=2*eps_xy
    # d = (0,0,0,2,2,2) => eta has shear entries = xi * d / 2 = 0.03
    expected_shear = np.array([
        [0, 0.03, 0.03],
        [0.03, 0, 0.03],
        [0.03, 0.03, 0],
    ])
    assert np.allclose(E_shear, expected_shear, atol=1e-10)
    print(f"  Shear strain matches xi*d/2 within 1e-10 ✅")

    print("\nAll assertions passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
