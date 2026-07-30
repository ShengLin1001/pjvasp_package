"""
mymetal.build.film

This subpackage provides functions for building and analyzing thin-film (slab)
structures, including generation, stretching, primitive/cubic-cell finding,
heterostructure construction, and surface passivation.

Modules:
    - stretch: Generates thin-film slabs from bulk and stretches the unit cell
      along specified directions.
    - extrfilm: Extracts lattice parameters, thickness, total energy and
      surface area from a film structure.
    - findprim: Finds the primitive cell of a thin-film structure via spglib.
    - findcubic: Finds the optimal cubic-like in-plane cell for a film.
    - findhetero: Builds heterostructures from two films via the hetbuilder
      package (optional dependency).
    - hydroxyl: Adds hydrogen/adsorbates to passivate dangling bonds on slab
      surfaces.
"""
