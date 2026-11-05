#!/usr/bin/env python3
# J. Pei, 2026-11-04
# pei_vasp_plot_thermal_expansion.py
#  >>> Post-process a y_thermal_expansion tree (LAMMPS): production averages per
#      temperature -> p_post_thermal_expansion.txt + .pdf. Plots this tree only; comparing
#      several potentials or experiment is left to the caller. Logic:
#      mymetal.post.thermal_expansion.

import argparse
import subprocess
import sys
from pathlib import Path

from mymetal.post.thermal_expansion import post_thermal_expansion

EPILOG = """\
most common:
  pei_vasp_plot_all -thermal_expansion           # every y_thermal_expansion below
  cd y_thermal_expansion && pei_vasp_plot_thermal_expansion.py

writes into the tree:
  p_post_thermal_expansion.txt   one row per finished temperature + the 0 K relaxation
  p_post_thermal_expansion.pdf   |a1| |a3|, |a3|/|a1|, V, E, MSD, CNA fractions vs T
  y_post_data.txt, y_dir/*/CONTCAR   pei_lmp_univ_post (last snapshot of each case)

notes:
  run it from INSIDE the tree; unfinished cases are listed and left out.
  -skip_univ_post skips pei_lmp_univ_post (CONTCARs are only for viewing).
"""

parser = argparse.ArgumentParser(
    description="Plot lattice parameters vs temperature from a y_thermal_expansion tree.",
    epilog=EPILOG, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-skip_univ_post", action="store_true",
                    help="Do not run pei_lmp_univ_post first.")
args = parser.parse_args()
sys.stdout.reconfigure(line_buffering=True)

path_work = Path.cwd().resolve()
print("🟢 ===== pei_vasp_plot_thermal_expansion in %s" % path_work)
if not args.skip_univ_post:
    # its exit code is 1 while cases are unfinished; the table still covers the finished ones
    subprocess.run(["pei_lmp_univ_post"], cwd=path_work)

df, lunfinished = post_thermal_expansion(path_work)

print("================ 📊 summary")
print("✅ finished   : %d" % (len(df) - 1))
print("❌ unfinished : %d" % len(lunfinished))
for name in lunfinished:
    print("   ❌ %s" % name)
df_md = df[df["T_set"] > 0]
print("📍 hdnnp EW in production: %s" % (" ".join("%d:%g" % (t, ew) for t, ew in
                                                 zip(df_md["T_set"], df_md["ew_prod"]))))
sys.exit(0 if not lunfinished else 1)
