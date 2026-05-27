LAMMPS Workflows
================

LAMMPS templates live in ``lmp_utils/``. The current repository uses them as
starting points for material-mechanics workflows such as:

* stretch calculations
* generalized stacking fault energy calculations
* elastic-constant energy calculations

The templates are intentionally separate from the ``mymetal`` Python package:
edit them for the target potential, unit style and cluster scheduler before
production runs.
