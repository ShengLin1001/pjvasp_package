mymetal Manual
--------------

Version: 1.0.0

``mymetal`` is the Python package at the core of ``pjvasp_package``.  The
manual is organized like a scientific-code manual rather than a single landing
page: start with the User Guide, then move to workflow-specific chapters, and
use the reference sections when you need API or script details.

Manual contents
===============

The sections below mirror the layout of a simulation-code manual:

* **User Guide** explains the project, installation, quick starts and examples.
* **Workflow Guide** documents the VASP, LAMMPS, SLURM and n2p2 workflows.
* **Reference** collects API pages, script directories, dependencies and
  development notes.

.. toctree::
   :maxdepth: 2
   :numbered:
   :caption: User Guide

   user_guide/overview
   user_guide/install
   user_guide/quickstart
   user_guide/examples
   user_guide/troubleshooting

.. toctree::
   :maxdepth: 2
   :numbered:
   :caption: Workflow Guide

   manual/workflows
   manual/vasp
   manual/lammps
   manual/slurm
   manual/n2p2

.. toctree::
   :maxdepth: 2
   :caption: Reference

   api
   reference/scripts
   reference/dependencies
   reference/development



Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
