"""
mymetal.ml.n2p2.sfparamgen

Symmetry-function parameter generation for n2p2, vendored from the n2p2
``symfunc_paramgen`` tool. It provides the ``SymFuncParamGenerator`` class for
generating, storing and writing radial/angular symmetry-function parameter sets
in the format required by n2p2.

Modules:
    - sfparamgen: The ``SymFuncParamGenerator`` class and its helpers.

The ``generate-acsf-*.py`` files in this directory are standalone executable
scripts (their file names contain hyphens and are not importable as modules);
they drive ``SymFuncParamGenerator`` for radial and narrow/wide angular ACSF
sets.
"""
