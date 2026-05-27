Workflow Guide
--------------

The Python package and shell workflow directories serve different roles:

* ``mymetal/`` supplies importable structure, I/O, calculation and
  post-processing utilities.
* ``vasp_utils/`` and ``myvasp/`` supply scripts for preparing, submitting and
  analysing VASP calculation directories.
* ``slurm_utils/`` provides scheduler-oriented templates.
* ``mymetal/example/`` contains notebooks and small calculation datasets.

Typical VASP directory layout
-----------------------------

Many post-processing and batch scripts expect a reference calculation plus
derived structures:

.. code-block:: text

   project/
   |-- y_full_relax/
   |-- y_src/
   |   |-- INCAR
   |   |-- KPOINTS
   |   |-- POTCAR
   |   `-- poscars2/
   `-- y_dir/
       |-- 0.997/
       |-- 1.000/
       `-- 1.003/

Example calculation path
------------------------

1. Build a bulk, film or strained series with ``mymetal.build``.
2. Write ``POSCAR`` files through :mod:`mymetal.io.vasp`.
3. Submit jobs with the relevant cluster-specific workflow script.
4. Summarise completed VASP outputs with :mod:`mymetal.post`.
5. Derive properties through :mod:`mymetal.calculate`, or export data for
   machine-learning workflows through :mod:`mymetal.ml.n2p2.dataset`.

The submission helpers are environment-specific: verify scheduler commands,
VASP paths and pseudopotential inputs on the target HPC system before running
batch jobs.
