#!/public3/home/scg6928/mysoft/env/pyenv/dft/bin/python
# pei_vasp_plot_convergence.py
#  >>> Post-process every scan axis of a finished y_convergence tree.

import argparse
import os

from mymetal.post.convergence import post_convergence

# Takes no options. It is run from INSIDE y_convergence (that is how pei_vasp_plot_all
# calls it) and chdir's one level up, which is where post_convergence expects to start.
# Parse first, so --help answers without moving the caller's cwd.
EPILOG = """\
most common:
  pei_vasp_plot_all -convergence                    # usual entry point: every y_convergence found
  cd y_convergence && pei_vasp_plot_convergence.py  # post-process a single tree by hand

notes:
  takes no options. run it from INSIDE y_convergence: it chdir's one level up, to the
  level of y_full_relax, and post_convergence reads the tree from there.
  it auto-detects which scan axes the tree actually carries and processes each into its
  OWN subtree, so a run that only scanned one axis needs no extra flags:
    y_convergence_encuts     ENCUT curve            -> p_post_convergence.{txt,pdf}
    y_convergence_kpoints    explicit k-mesh curve  -> p_post_convergence.{txt,pdf}
    y_convergence_klength    automatic k-mesh curve -> p_post_convergence.{txt,pdf}
    y_convergence_ediff      EDIFF curve (log x)    -> p_post_convergence.{txt,pdf}
    y_convergence_kpar_ncore parallel timing        -> p_post_kpar_ncore.{txt,pdf}
  run it only once the jobs have finished.
"""

parser = argparse.ArgumentParser(
    description="Post-process every scan axis of a finished y_convergence tree.",
    epilog=EPILOG, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.parse_args()

os.chdir('..')
post_convergence()