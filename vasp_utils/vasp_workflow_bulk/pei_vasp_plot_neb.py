#!/public3/home/scg6928/mysoft/env/pyenv/dft/bin/python
# pei_vasp_plot_neb.py
#  >>> Plot the migration barrier of a finished y_neb tree.

import argparse
import os

from mymetal.post.neb import post_neb

# It is run from INSIDE y_neb (that is how pei_vasp_plot_all calls it) and chdir's one
# level up, which is where post_neb expects to start. Parse first, so --help answers
# without moving the caller's cwd.
EPILOG = """\
most common:
  pei_vasp_plot_all -neb                 # usual entry point: every y_neb found
  cd y_neb && pei_vasp_plot_neb.py       # plot a single tree by hand
  cd y_neb && pei_vasp_plot_neb.py -calc lammps
                                         # LAMMPS NEB (pei_lmp_univ_post -neb, no OUTCAR)

notes:
  run it from INSIDE y_neb: it chdir's one level up, to the level of y_full_relax, and
  reads the image dirs from there. -calc lammps skips the OUTCAR-trajectory plot (p_neb_full).
  run it only once the NEB images have finished.
"""

parser = argparse.ArgumentParser(
    description="Plot the NEB migration barrier from a finished y_neb tree.",
    epilog=EPILOG, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-calc", default="vasp", choices=["lammps", "vasp"],
                    help="Engine that ran the NEB (default vasp).")
args = parser.parse_args()

# the tree is named after the cwd (y_neb, or any -outdir name)
dirsurf = os.path.basename(os.getcwd())
os.chdir('..')
post_neb(dirsurf=dirsurf, calc=args.calc)