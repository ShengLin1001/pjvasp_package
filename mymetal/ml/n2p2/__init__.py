"""
mymetal.n2p2

This module provides a wrapper for the N2P2 neural network potential.
It includes functions to generate datasets from VASP OUTCAR files and to read datasets in NNP format.

Submodules:
    - dataset: Contains functions to generate datasets from VASP OUTCAR files and to read datasets in NNP format.
    - workflow: High-level PeiN2p2 orchestrator chaining the full n2p2 training/testing pipeline.
"""

from mymetal.ml.n2p2.workflow import PeiN2p2
