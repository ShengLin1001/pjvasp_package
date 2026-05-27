Development Notes
=================

Build documentation locally
---------------------------

.. code-block:: bash

   python -m pip install -r docs/requirements.txt
   python -m pip install --no-deps "git+https://github.com/ShengLin1001/myalloy_package.git@master"
   python -m sphinx -b html -W --keep-going docs/source docs/_build/html

Validation
----------

For Python code changes, run:

.. code-block:: bash

   python -m compileall mymetal

For documentation changes, keep the Sphinx build warning-free because the
GitHub Pages workflow uses ``-W --keep-going``.
