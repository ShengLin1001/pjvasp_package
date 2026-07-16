#!/public3/home/scg6928/mysoft/env/pyenv/dft/bin/python
# pei_vasp_plot_neb.py
#  >>> Plot the migration barrier of a finished y_neb tree.

import argparse
import os

from mymetal.post.neb import post_neb

# Takes no options. It is run from INSIDE y_neb (that is how pei_vasp_plot_all calls it)
# and chdir's one level up, which is where post_neb expects to start. Parse first, so
# --help answers without moving the caller's cwd.
EPILOG = """\
most common:
  pei_vasp_plot_all -neb                 # usual entry point: every y_neb found
  cd y_neb && pei_vasp_plot_neb.py       # plot a single tree by hand

notes:
  takes no options. run it from INSIDE y_neb: it chdir's one level up, to the level of
  y_full_relax, and reads the image dirs from there.
  run it only once the NEB images have finished.
"""

parser = argparse.ArgumentParser(
    description="Plot the NEB migration barrier from a finished y_neb tree.",
    epilog=EPILOG, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.parse_args()

os.chdir('..')
post_neb()