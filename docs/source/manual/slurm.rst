SLURM and Scheduler Helpers
===========================

Scheduler templates live under ``slurm_utils/`` and related workflow
directories. They are environment-specific by design.

Recommended checks before submission
------------------------------------

* Confirm the target partition, account and walltime.
* Confirm VASP or LAMMPS executable paths.
* Confirm that source directories contain ``INCAR``, ``KPOINTS``, ``POTCAR``
  and submit scripts when required.
* Run post-processing commands only after expected output files exist.
