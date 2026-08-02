#!/usr/bin/env python3
"""Regenerate deterministic images used by the Sphinx manual.

The script imports the ``render_*`` / ``run_example`` helpers from each
documentation example and writes the resulting PNG into
``docs/source/_static/images/generated/``. It is invoked by the CI workflow
after the examples themselves have run, and is also the entry point for
local documentation builds.
"""

from __future__ import annotations

import sys
from pathlib import Path


PATH_DOCS = Path(__file__).resolve().parents[1]
PATH_EXAMPLES = PATH_DOCS / "examples"
PATH_IMAGES = PATH_DOCS / "source" / "_static" / "images" / "generated"
sys.path.insert(0, str(PATH_EXAMPLES))


def main() -> None:
    PATH_IMAGES.mkdir(parents=True, exist_ok=True)

    # Au(111) slab: reused by the index page and the getting-started page.
    from getting_started_au111 import build_au111_slab, render_slab
    render_slab(build_au111_slab(), PATH_IMAGES / "au111_slab.png")
    print("generated: " + str(PATH_IMAGES / "au111_slab.png"))

    # FCC Cu(100)/(110)/(111) comparison.
    from fcc_surfaces import build_fcc_slabs, render_comparison
    render_comparison(build_fcc_slabs(), PATH_IMAGES / "fcc_surfaces.png")
    print("generated: " + str(PATH_IMAGES / "fcc_surfaces.png"))

    # HCP Mg low-index surfaces comparison.
    from hcp_surfaces import build_hcp_cells, render_comparison as render_hcp
    render_hcp(build_hcp_cells(), PATH_IMAGES / "hcp_surfaces.png")
    print("generated: " + str(PATH_IMAGES / "hcp_surfaces.png"))

    # EOS curve figure.
    from eos_curve import make_synthetic_ev, fit_with_ase, fit_murnaghan_with_b0p, render_eos_figure
    volumes, energies, meta = make_synthetic_ev()
    fit_murn = fit_with_ase(volumes, energies, eos="murnaghan")
    fit_bm = fit_with_ase(volumes, energies, eos="birchmurnaghan")
    fit_murn_b0p = fit_murnaghan_with_b0p(volumes, energies)
    render_eos_figure(
        volumes, energies, fit_murn, fit_bm, fit_murn_b0p, meta,
        PATH_IMAGES / "eos_curve.png",
    )
    print("generated: " + str(PATH_IMAGES / "eos_curve.png"))

    # Biaxial stretch figure.
    from biaxial_stretch import build_slab, stretch_along, synthesize_energy, render_figure, DIRECTIONS
    slab = build_slab()
    stretched_per_direction = {d: stretch_along(slab, d) for d in DIRECTIONS}
    energies_per_direction = {
        d: synthesize_energy(stretched_per_direction[d], slab, d, seed=20260729 + DIRECTIONS.index(d))
        for d in DIRECTIONS
    }
    render_figure(
        slab, stretched_per_direction, energies_per_direction,
        PATH_IMAGES / "biaxial_stretch.png",
    )
    print("generated: " + str(PATH_IMAGES / "biaxial_stretch.png"))

    # Bulk structures comparison (FCC / BCC / HCP / diamond).
    from bulk_structures import build_bulk_cells, render_comparison
    render_comparison(build_bulk_cells(), PATH_IMAGES / "bulk_structures.png")
    print("generated: " + str(PATH_IMAGES / "bulk_structures.png"))

    # k-point sampling figure (Monkhorst-Pack vs Gamma + RK scan).
    from kpoints_sampling import build_slab as build_kp_slab, get_mp_and_gamma, rk_scan, compare_reciprocal, render_figure as render_kp
    slab = build_kp_slab()
    mpk, gk = get_mp_and_gamma()
    rk_result = rk_scan(slab, [20, 40, 60, 80, 100, 120])
    render_kp(mpk, gk, rk_result, PATH_IMAGES / "kpoints_sampling.png")
    print("generated: " + str(PATH_IMAGES / "kpoints_sampling.png"))

    # FCC Schmid factor figure.
    from schmid_factor import compute_reference_table, compute_scan_table, compute_polar_grid, render_figure as render_schmid
    render_schmid(
        compute_reference_table(),
        compute_scan_table(),
        compute_polar_grid(),
        PATH_IMAGES / "schmid_factor.png",
    )
    print("generated: " + str(PATH_IMAGES / "schmid_factor.png"))

    # Neighbor distance / RDF figure.
    from neighbor_distances import build_structures, collect_distances, build_rdf, render_figure as render_neighbor
    structures = build_structures()
    rdfs = [build_rdf(collect_distances(atoms)) for atoms in structures]
    render_neighbor(rdfs, PATH_IMAGES / "neighbor_distances.png")
    print("generated: " + str(PATH_IMAGES / "neighbor_distances.png"))

    # Atom manipulation (move + fix + selective dynamics) figure.
    from atom_manipulation import build_slab as build_am_slab, fix_bottom_half, shift_top_layer, render_figure as render_am
    slab = build_am_slab()
    z_min = float(slab.positions[:, 2].min())
    z_max = float(slab.positions[:, 2].max())
    z_mid = 0.5 * (z_min + z_max)
    mask = [atom.position[2] < z_mid for atom in slab]
    slab_constrained = fix_bottom_half(slab)
    slab_shifted = shift_top_layer(slab_constrained)
    render_am(slab, slab_constrained, slab_shifted, mask, PATH_IMAGES / "atom_manipulation.png")
    print("generated: " + str(PATH_IMAGES / "atom_manipulation.png"))

    # Strain / deformation matrix figure.
    from strain_deformation import build_reference_cell, make_uniaxial_x_strain, make_biaxial_xy_strain, make_simple_shear_xy, compute_strain_pack, render_figure as render_strain
    cell_ref = build_reference_cell()
    packs = [
        compute_strain_pack(cell_ref, make_uniaxial_x_strain(cell_ref, 0.05)),
        compute_strain_pack(cell_ref, make_biaxial_xy_strain(cell_ref, 0.05)),
        compute_strain_pack(cell_ref, make_simple_shear_xy(cell_ref, 0.10)),
    ]
    titles = [
        "Uniaxial x\nstrain = +0.05",
        "Biaxial xy\nstrain = +0.05",
        "Simple shear xy\ngamma = +0.10",
    ]
    render_strain(packs, titles, PATH_IMAGES / "strain_deformation.png")
    print("generated: " + str(PATH_IMAGES / "strain_deformation.png"))

    # GSFE model comparison (FCC_111 / HCP_basal / HCP_prism1w).
    from gsfe_models import build_gsfe_models, render_comparison as render_gsfe
    render_gsfe(build_gsfe_models(), PATH_IMAGES / "gsfe_models.png")
    print("generated: " + str(PATH_IMAGES / "gsfe_models.png"))

    # Cubic cell finding + uniaxial stretch.
    from cubic_cell_and_stretch import (
        build_primitive_film, build_orthorhombic_film, build_stretched_films,
        render_figure as render_cubic,
    )
    prim = build_primitive_film()
    ortho = build_orthorhombic_film(prim)
    stretched = build_stretched_films(ortho)
    render_cubic(prim, ortho, stretched, PATH_IMAGES / "cubic_cell_and_stretch.png")
    print("generated: " + str(PATH_IMAGES / "cubic_cell_and_stretch.png"))

    # Deformation matrix + Hermite Normal Form.
    from deformation_and_hnf import build_deformation_pair, build_hnf_table, render_figure as render_deform
    render_deform(PATH_IMAGES / "deformation_and_hnf.png", build_deformation_pair(), build_hnf_table())
    print("generated: " + str(PATH_IMAGES / "deformation_and_hnf.png"))

    # Reciprocal lattice vectors + RK-based k-point mesh.
    from reciprocal_lattice import build_bulk_cells, reciprocal_pair, summarize_cell, render_figure as render_recip
    recip_rows = []
    for idx, atoms in enumerate(build_bulk_cells()):
        info = summarize_cell(idx, atoms)
        recip_rows.append(info)
    render_recip(recip_rows, PATH_IMAGES / "reciprocal_lattice.png")
    print("generated: " + str(PATH_IMAGES / "reciprocal_lattice.png"))


if __name__ == "__main__":
    main()
