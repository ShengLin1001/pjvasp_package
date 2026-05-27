VASP Workflows
==============

VASP-facing tools are spread across ``mymetal``, ``vasp_utils`` and
``myvasp``:

``mymetal.io.vasp``
   Python-level POSCAR/CONTCAR readers and writers.

``mymetal.post``
   Python post-processing functions for OUTCAR-heavy calculation directories.

``vasp_utils/vasp_universal``
   Reusable shell helpers for submission, cleanup, continuation and result
   collection.

``myvasp``
   Historical shell scripts and VTST-style helper scripts. Importable
   ``myvasp`` Python helpers come from ``myalloy_package``.

Before running VASP workflow scripts, verify scheduler commands, pseudopotential
paths and VASP executable paths on the target cluster.
