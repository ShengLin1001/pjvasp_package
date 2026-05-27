Dependencies
============

Core scientific dependencies
----------------------------

The main package uses ``numpy``, ``scipy``, ``ase``, ``spglib``, ``pandas`` and
``matplotlib`` across the common workflows.

Companion package
-----------------

``myalloy`` and ``myvasp`` are supplied by:

.. code-block:: text

   https://github.com/ShengLin1001/myalloy_package

Install with:

.. code-block:: bash

   python -m pip install "git+https://github.com/ShengLin1001/myalloy_package.git@master"

Optional dependencies
---------------------

Some functions require optional packages such as ``hetbuilder``, ``ovito``,
``pymatgen``, ``torch`` or ``torchvision``. The manual may mock these imports
when the package is only needed for runtime execution rather than for API
documentation.
