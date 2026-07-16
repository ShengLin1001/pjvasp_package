#!/public3/home/scg6928/mysoft/env/pyenv/dft/bin/python
# pei_vasp_plot_stretch.py
#  >>> Plot the strain-energy curve of a finished y_stretch tree and fit the
#      equilibrium lattice constant.

import argparse
import os

from mymetal.post.stretch import post_stretch

# Takes no options. It is run from INSIDE y_stretch (that is how pei_vasp_plot_all calls
# it) and chdir's one level up, which is where post_stretch expects to start. Parse first,
# so --help answers without moving the caller's cwd.
EPILOG = """\
most common:
  pei_vasp_plot_all -stretch                    # usual entry point: every y_stretch found
  cd y_stretch && pei_vasp_plot_stretch.py      # plot a single tree by hand

notes:
  takes no options. run it from INSIDE y_stretch: it chdir's one level up, to the level
  of y_full_relax, and reads y_dir/* from there.
  run it only once the stretched jobs have finished.
"""

parser = argparse.ArgumentParser(
    description="Plot the strain-energy curve from a finished y_stretch tree.",
    epilog=EPILOG, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.parse_args()

os.chdir('..')
post_stretch()