Overview
========

``pjvasp_package`` is a research-oriented materials-computation workspace.
The importable Python package is ``mymetal``; surrounding directories provide
VASP, LAMMPS and scheduler workflow helpers.

Major components
----------------

``mymetal.build``
   Build bulk, film, surface and heterostructure models.

``mymetal.io``
   Read and write VASP-style structures and intermediate datasets.

``mymetal.calculate``
   Compute surface energies, strain/deformation quantities, lattice mismatch
   metrics and other reusable analysis values.

``mymetal.post``
   Summarise completed VASP calculations and extract convergence, energy,
   pressure, force and warning information from calculation directories.

``mymetal.ml``
   Prepare data for machine-learning-potential workflows, especially n2p2.

External helper packages
------------------------

Some VASP workflow functions depend on ``myalloy`` and ``myvasp`` from
``myalloy_package``:

.. code-block:: text

   https://github.com/ShengLin1001/myalloy_package

The online manual installs that package directly from GitHub when building the
API documentation.
