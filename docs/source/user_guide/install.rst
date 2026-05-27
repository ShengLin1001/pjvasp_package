Installation
============

Python package
--------------

Use Python 3.10 when possible:

.. code-block:: bash

   git clone https://github.com/ShengLin1001/pjvasp_package.git
   cd pjvasp_package
   python -m pip install -r requirements.txt
   python -m pip install -e .

Runtime helper package
----------------------

Several historical VASP helpers import ``myalloy`` or ``myvasp``. Install them
from the companion repository:

.. code-block:: bash

   python -m pip install "git+https://github.com/ShengLin1001/myalloy_package.git@master"

Documentation build
-------------------

The GitHub Pages workflow installs a smaller documentation environment plus
``myalloy_package`` without its transitive dependencies:

.. code-block:: bash

   python -m pip install -r docs/requirements.txt
   python -m pip install --no-deps "git+https://github.com/ShengLin1001/myalloy_package.git@master"
   python -m sphinx -b html -W --keep-going docs/source docs/_build/html
