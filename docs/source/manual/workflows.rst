Workflow Guide
==============

The package is designed around a repeated pattern:

1. Build or transform structures with ``mymetal.build``.
2. Write input files with ``mymetal.io`` and workflow scripts.
3. Submit calculations with the cluster-specific VASP/LAMMPS helpers.
4. Collect outputs with ``mymetal.post``.
5. Feed summarised data into plotting, mechanics, energy or ML workflows.

Common directory layout
-----------------------

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

This convention appears in many VASP batch and post-processing utilities.
