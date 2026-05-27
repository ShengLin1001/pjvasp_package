Installation
------------

Requirements
------------

``mymetal`` is developed as a research package for Python 3.10. Its main
scientific workflows use NumPy, SciPy, ASE, spglib, pymatgen and matplotlib.
Specialised functions may require additional packages such as ``hetbuilder``,
``ovito`` or PyTorch.

Install from source
-------------------

Clone the repository and install it in editable mode:

.. code-block:: bash

   git clone https://github.com/ShengLin1001/pjvasp_package.git
   cd pjvasp_package
   python -m pip install -r requirements.txt
   python -m pip install -e .

On HPC systems, choose a Python 3 environment compatible with the scientific
dependencies before installing. VASP executables, pseudopotentials and batch
scheduler commands are external to the Python package.

Build this documentation
------------------------

The online documentation is built and deployed by GitHub Actions. A local
documentation build only requires the small documentation dependency set:

.. code-block:: bash

   python -m pip install -r docs/requirements.txt
   python -m sphinx -b html -W --keep-going docs/source docs/_build/html

Open ``docs/_build/html/index.html`` after the build.
