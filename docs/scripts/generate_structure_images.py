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

    # Extended XYZ trajectory + general_read/general_write I/O.
    from io_extxyz_and_general import (
        build_trajectory, write_trajectory, read_trajectory,
        build_convergence_table, write_table, read_table,
        render_figure as render_io,
    )
    import tempfile
    with tempfile.TemporaryDirectory() as tmpd:
        tmpp = Path(tmpd)
        frames = build_trajectory()
        write_trajectory(frames, tmpp / "trajectory.xyz")
        frames_back = read_trajectory(tmpp / "trajectory.xyz")
        df = build_convergence_table()
        write_table(df, tmpp / "data.txt")
        df_read = read_table(tmpp / "data.txt")
        render_io(frames, df_read, PATH_IMAGES / "io_extxyz_and_general.png")
    print("generated: " + str(PATH_IMAGES / "io_extxyz_and_general.png"))

    # Cij elastic constant energy-strain fitting (synthetic data).
    from cij_energy_fitting import build_synthetic_data, fit_cij_from_modes, render_figure as render_cij
    modes = build_synthetic_data()
    fitted = fit_cij_from_modes(modes)
    render_cij(modes, fitted, PATH_IMAGES / "cij_energy_fitting.png")
    print("generated: " + str(PATH_IMAGES / "cij_energy_fitting.png"))

    # HCP Miller index 3<->4 conversion + density calculation.
    from miller_index_and_density import (
        build_miller_examples, build_bulk_cells, summarize_density,
        render_figure as render_miller,
    )
    miller_ex = build_miller_examples()
    density_rows = [summarize_density(idx, c) for idx, c in enumerate(build_bulk_cells())]
    render_miller(miller_ex, density_rows, PATH_IMAGES / "miller_index_and_density.png")
    print("generated: " + str(PATH_IMAGES / "miller_index_and_density.png"))

    # Periodic table heatmap + van Arkel triangle.
    from periodic_table_and_arkel import (
        build_formula_dict, build_arkel_materials,
        render_periodic_heatmap, render_arkel_triangle,
    )
    render_periodic_heatmap(build_formula_dict(), PATH_IMAGES / "periodic_table_heatmap.png")
    print("generated: " + str(PATH_IMAGES / "periodic_table_heatmap.png"))
    render_arkel_triangle(build_arkel_materials(), PATH_IMAGES / "van_arkel_triangle.png")
    print("generated: " + str(PATH_IMAGES / "van_arkel_triangle.png"))

    # Slurm script generation (dry-run, no sbatch).
    from slurm_script_generation import build_header, build_base_script, render_figure as render_slurm
    render_slurm(build_header(), build_base_script(PATH_IMAGES / "dummy.sh"),
                 PATH_IMAGES / "slurm_script_generation.png")
    print("generated: " + str(PATH_IMAGES / "slurm_script_generation.png"))

    # vasp_universal overview: runner flow, directory scan, clean-up comparison,
    # INCAR operations, exit-code convention.
    from vasp_universal_overview import render_overview as render_vasp_univ
    render_vasp_univ(PATH_IMAGES / "vasp_universal_overview.png")
    print("generated: " + str(PATH_IMAGES / "vasp_universal_overview.png"))

    # vasp_workflow_bulk + neb_utils overview: workflow classification,
    # lifecycle, plot_all dispatcher, NEB toolchain, strain types.
    from vasp_workflow_bulk_overview import render_overview as render_vasp_wf
    render_vasp_wf(PATH_IMAGES / "vasp_workflow_bulk_overview.png")
    print("generated: " + str(PATH_IMAGES / "vasp_workflow_bulk_overview.png"))

    # slurm_utils overview: three modes, directory discovery, retry strategy,
    # preset registry, monitor & useful commands.
    from slurm_utils_overview import render_overview as render_slurm
    render_slurm(PATH_IMAGES / "slurm_utils_overview.png")
    print("generated: " + str(PATH_IMAGES / "slurm_utils_overview.png"))

    # lmp_utils overview: runner flow, template comparison, .mod dependencies,
    # GSFE slip systems, sed template substitution.
    from lmp_utils_overview import render_overview as render_lmp
    render_lmp(PATH_IMAGES / "lmp_utils_overview.png")
    print("generated: " + str(PATH_IMAGES / "lmp_utils_overview.png"))

    # n2p2_utils overview: full workflow, active learning SF selection,
    # SF generation params, pipeline runner, clean_train strategy.
    from n2p2_utils_overview import render_overview as render_n2p2
    render_n2p2(PATH_IMAGES / "n2p2_utils_overview.png")
    print("generated: " + str(PATH_IMAGES / "n2p2_utils_overview.png"))

    # Plot gallery: one PNG per mymetal.universal.plot module, all VASP-free.
    from plot_gallery_demo import (
        render_general, render_plot, render_colorbar,
        render_convergence, render_relax_convergence, render_kpar_ncore,
        render_stretch, render_energy_components,
        render_interlayer_distance, render_zpositions, render_rdf,
        render_learning_curve, render_compare, render_rmse_by_tag,
        render_epoch_rmse, render_pretty_plot, render_periodic_heatmap,
        render_van_arkel, render_render_info, render_ppt_info, render_dos,
    )
    render_general(PATH_IMAGES / "plot_gallery_general.png")
    print("generated: " + str(PATH_IMAGES / "plot_gallery_general.png"))
    render_plot(PATH_IMAGES / "plot_gallery_plot.png")
    print("generated: " + str(PATH_IMAGES / "plot_gallery_plot.png"))
    render_colorbar(PATH_IMAGES / "plot_gallery_colorbar.png")
    print("generated: " + str(PATH_IMAGES / "plot_gallery_colorbar.png"))
    render_convergence(PATH_IMAGES / "plot_gallery_convergence.png")
    print("generated: " + str(PATH_IMAGES / "plot_gallery_convergence.png"))
    render_relax_convergence(PATH_IMAGES / "plot_gallery_relax.png")
    print("generated: " + str(PATH_IMAGES / "plot_gallery_relax.png"))
    render_kpar_ncore(PATH_IMAGES / "plot_gallery_kpar_ncore.png")
    print("generated: " + str(PATH_IMAGES / "plot_gallery_kpar_ncore.png"))
    render_stretch(PATH_IMAGES / "plot_gallery_stretch.png")
    print("generated: " + str(PATH_IMAGES / "plot_gallery_stretch.png"))
    render_energy_components(PATH_IMAGES / "plot_gallery_energy.png")
    print("generated: " + str(PATH_IMAGES / "plot_gallery_energy.png"))
    render_interlayer_distance(PATH_IMAGES / "plot_gallery_interlayer.png")
    print("generated: " + str(PATH_IMAGES / "plot_gallery_interlayer.png"))
    render_zpositions(PATH_IMAGES / "plot_gallery_zpositions.png")
    print("generated: " + str(PATH_IMAGES / "plot_gallery_zpositions.png"))
    render_rdf(PATH_IMAGES / "plot_gallery_rdf.png")
    print("generated: " + str(PATH_IMAGES / "plot_gallery_rdf.png"))
    render_learning_curve(PATH_IMAGES / "plot_gallery_learning.png")
    print("generated: " + str(PATH_IMAGES / "plot_gallery_learning.png"))
    render_compare(PATH_IMAGES / "plot_gallery_compare.png")
    print("generated: " + str(PATH_IMAGES / "plot_gallery_compare.png"))
    render_rmse_by_tag(PATH_IMAGES / "plot_gallery_rmse_tag.png")
    print("generated: " + str(PATH_IMAGES / "plot_gallery_rmse_tag.png"))
    render_epoch_rmse(PATH_IMAGES / "plot_gallery_epoch_rmse.png")
    print("generated: " + str(PATH_IMAGES / "plot_gallery_epoch_rmse.png"))
    render_pretty_plot(PATH_IMAGES / "plot_gallery_pretty.png")
    print("generated: " + str(PATH_IMAGES / "plot_gallery_pretty.png"))
    render_periodic_heatmap(PATH_IMAGES / "plot_gallery_periodic.png")
    print("generated: " + str(PATH_IMAGES / "plot_gallery_periodic.png"))
    render_van_arkel(PATH_IMAGES / "plot_gallery_arkel.png")
    print("generated: " + str(PATH_IMAGES / "plot_gallery_arkel.png"))
    render_render_info(PATH_IMAGES / "plot_gallery_render.png")
    print("generated: " + str(PATH_IMAGES / "plot_gallery_render.png"))
    render_ppt_info(PATH_IMAGES / "plot_gallery_ppt.png")
    print("generated: " + str(PATH_IMAGES / "plot_gallery_ppt.png"))
    render_dos(PATH_IMAGES / "plot_gallery_dos.png")
    print("generated: " + str(PATH_IMAGES / "plot_gallery_dos.png"))

    # Real server-data figures used by the VASP/LAMMPS/n2p2 manual pages.
    from vasp_workflow_real_data import main as render_vasp_workflow_real
    from lammps_workflow_real_data import main as render_lammps_workflow_real
    from n2p2_workflow_real_data import main as render_n2p2_workflow_real
    render_vasp_workflow_real(PATH_IMAGES)
    render_lammps_workflow_real(PATH_IMAGES)
    render_n2p2_workflow_real(PATH_IMAGES)


if __name__ == "__main__":
    main()
