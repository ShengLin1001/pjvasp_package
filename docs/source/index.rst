mymetal Documentation
---------------------

.. toctree::
   :maxdepth: 2
   :caption: Documentation

   installation
   quickstart
   workflows
   api
   autoapi/index

``mymetal`` is the Python package at the core of ``pjvasp_package``. It
provides reusable tools for computational materials workflows: structure
construction, VASP-oriented I/O, property calculations, post-processing,
plotting utilities, and preparation of machine-learning-potential data.

The repository also contains HPC workflow scripts in ``vasp_utils/``,
``myvasp/`` and ``slurm_utils/``. Those scripts depend on the job scheduler
and directory conventions of the target cluster; the API reference here
focuses on the reusable Python package.

Key areas
---------

* :doc:`quickstart` introduces structure generation, VASP output, surface
  energy calculation and n2p2 dataset conversion.
* :doc:`workflows` describes how the package fits into VASP batch workflows.
* :doc:`api` covers the principal public modules from current Python
  docstrings.



Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
