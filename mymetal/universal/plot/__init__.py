"""
plot module

This module contains the plot class and its methods.

Submodules:
    - plot: The package's house-style plot helper (my_plot and friends).
    - general: General plot functions.
    - render: Renders images from OVITO pipelines.
    - atominfo: Plots interlayer distances, z-positions and radial
      distribution functions.
    - energy: Plots energy components with polynomial fits and comparisons.
    - workflow: Plots workflow results (convergence, Cij energy, KPAR/NCORE,
      relax convergence, NEB, stretch).
    - n2p2: Plots for n2p2 training post-processing (learning curve, DFT-vs-NNP
      scatter, per-tag RMSE).
    - plotting: Utilities for nicer matplotlib plots (vendored from pymatgen).
    - oldplotdos: Legacy DOS/COHP plotting helpers.
"""