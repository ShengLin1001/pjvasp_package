#!/public3/home/scg6928/mysoft/env/pyenv/dft/bin/python
# pei_vasp_plot_convergence.py
#  >>> Plot the ENCUT / KPOINTS convergence curves of a finished y_convergence tree.

import argparse
import os

from mymetal.post.convergence import post_convergence

# Takes no options. It is run from INSIDE y_convergence (that is how pei_vasp_plot_all
# calls it) and chdir's one level up, which is where post_convergence expects to start.
# Parse first, so --help answers without moving the caller's cwd.
EPILOG = """\
most common:
  pei_vasp_plot_all -convergence                    # usual entry point: every y_convergence found
  cd y_convergence && pei_vasp_plot_convergence.py  # plot a single tree by hand

notes:
  takes no options. run it from INSIDE y_convergence: it chdir's one level up, to the
  level of y_full_relax, and reads y_convergence_encuts / y_convergence_kpoints from there.
  run it only once the convergence jobs have finished.
"""

parser = argparse.ArgumentParser(
    description="Plot ENCUT/KPOINTS convergence curves from a finished y_convergence tree.",
    epilog=EPILOG, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.parse_args()

os.chdir('..')
post_convergence()