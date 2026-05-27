Troubleshooting
===============

Import errors for ``myvasp`` or ``myalloy``
-------------------------------------------

Install the companion package:

.. code-block:: bash

   python -m pip install "git+https://github.com/ShengLin1001/myalloy_package.git@master"

The local repository also contains a ``myvasp/`` workflow-script directory.
The importable helper package is supplied by ``myalloy_package``.

Missing HPC commands
--------------------

On the CentOS HPC platform, check available modules before installing tools
manually:

.. code-block:: bash

   module avail

Documentation warnings
----------------------

The online manual treats warnings as errors. If a new API module imports an
optional runtime package, add the package to ``docs/requirements.txt`` or list
it in ``autodoc_mock_imports`` when it is not needed for documentation.
