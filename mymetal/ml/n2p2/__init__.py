"""
mymetal.n2p2

This module provides a wrapper for the N2P2 neural network potential.
It includes functions to generate datasets from VASP OUTCAR files and to read datasets in NNP format.

Submodules:
    - dataset: Contains functions to generate datasets from VASP OUTCAR files and to read datasets in NNP format.
    - workflow: High-level PeiN2p2 orchestrator chaining the full n2p2 training/testing pipeline.
    - workflow_md: PeiN2p2MD subclass driving MD active learning (LAMMPS heating runs ->
      committee divergence via nnp-dataset -> DFT labelling -> back into the training set).
"""

from mymetal.ml.n2p2.workflow import PeiN2p2
from mymetal.ml.n2p2.workflow_md import PeiN2p2MD
