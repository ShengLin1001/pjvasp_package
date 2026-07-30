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


if __name__ == "__main__":
    main()
