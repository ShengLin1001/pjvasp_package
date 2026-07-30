"""
mymetal.ml.n2p2.calculate

Calculation-side helpers for n2p2 symmetry-function (SF) analysis and training
post-processing.

Modules:
    - sf: Computes radial (g2) and angular (g3/g9) symmetry functions for
      atomic structures via ASE, and writes SF settings/parameter strings.
    - cur: CUR-decomposition-based selection of the most informative symmetry
      functions from per-SF ``function.data`` files produced by nnp-scaling.
    - post: Parses n2p2 ``nnp-train`` output (learning-curve, trainpoints,
      trainforces) and converts normalized units back to physical units.
"""
