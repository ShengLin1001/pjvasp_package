#!/public3/home/scg6928/mysoft/env/pyenv/dft/bin/python

# Plot the per-ionic-step relaxation convergence (energy/atom + max force) for
# every y_dir sub-directory. The data files are produced by the bash helper
# pei_vasp_univ_extract_convergence; here we only read and plot them.
#
# Run from the directory that contains y_dir (same place as pei_vasp_univ_post).
#
# J. Pei, 2026-07-05

import matplotlib
matplotlib.use('Agg')

from mymetal.post.relax_convergence import my_univ_post_convergence

my_univ_post_convergence()
