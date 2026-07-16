#!/public3/home/scg6928/mysoft/env/pyenv/dft/bin/python
# J. Pei, 2025-12-07
import os
import argparse
import shutil
import numpy as np
from mymetal.build.workflow.general import cp_vaspfiles
from mymetal.build.film.stretch import stretch_list_along_direction_to_cell, adjust_lattice_for_volume_conservation
from mymetal.build.workflow.general import cp_vaspfiles, rm_i, compare_three_lattices
from mymetal.io.vasp import my_write_vasp, my_read_vasp
from mymetal.universal.print.print import confirm_prepare_outdir

# This script is to run eos calculations for stretched structures
# parameters
# -isif 4 or 2(default)
# -init 0.94 
# -final 1.06
# -interval 0.02
# -ratios None  , here is the volume change ratio
# --init/--final/--interval are VOLUME ratios here (not strains as in pei_vasp_run_stretch.py),
# and the cell is scaled by ratio^(1/3) in all three directions.
EPILOG = """\
most common:
  pei_vasp_run_eos.py                             # default V/V0 = 0.94 .. 1.06, step 0.02 (7 pts)
  pei_vasp_run_eos.py --init 0.90 --final 1.10 --interval 0.01
                                                  # wider + denser, for a stiffer fit
  pei_vasp_run_eos.py --isif 2                    # fixed cell shape (default 4 relaxes it)
  pei_vasp_run_eos.py --ratios 0.94,0.96,0.98,1.0,1.02,1.04,1.06
                                                  # explicit list instead of init/final/interval

then:
  submit y_eos/y_dir/* with your usual flow (pei_vasp_univ_sbatch_jobs), then fit E(V).

notes:
  run at the level of y_full_relax; values are VOLUME ratios (0.94 = 94% of V0), and each
  cell is scaled by ratio^(1/3) along xyz. --ratios overrides --init/--final/--interval.
  it only prepares dirs -- it never sbatches.
  TODO: pei_vasp_plot_all has no -eos option yet, so the fit is still manual.
"""

parser = argparse.ArgumentParser(
    description="Run EOS calculations for stretched structures.",
    epilog=EPILOG, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--isif", type = int, choices=[2,4], default = 4, help="ISIF tag in VASP calculations.")
parser.add_argument("--init", type = float, default = 0.94,      help = "Initial volume ratio (e.g., 0.94 for 94%% of the original volume)")
parser.add_argument("--final", type = float, default = 1.06,    help = "Final   volume ratio (e.g., 1.06 for 106%% of the original volume)")
parser.add_argument("--interval", type = float, default = 0.02, help = "Volume ratio interval (e.g., 0.02 for 2%% change in volume)")
# other methods to set up stretch list
parser.add_argument("--ratios", type = str, default = None, help = "Comma-separated volume ratio list (e.g., 0.94,0.96,0.98,1.0,1.02,1.04,1.06)")
parser.add_argument("--deleteold", action = "store_true", help="Whether to delete existing directory automatically.")
args = parser.parse_args()

print("ISIF                  :", args.isif)
print("Initial volume ratio  :", args.init)
print("Final   volume ratio  :", args.final)
print("Volume ratio interval :", args.interval)
print("Volume ratio list     :", args.ratios)
print("Delete old            :", args.deleteold)

myroot = os.getcwd()
os.makedirs("y_full_relax_temp", exist_ok=True)
cp_vaspfiles(cpfiles = ["CONTCAR", "KPOINTS", "INCAR", "POTCAR", "Y_CONSTR_LATT", "CHGCAR"], oldir = "y_full_relax", newdir = "y_full_relax_temp")

# ratio lists
if args.ratios is not None:
    ratio_list = [float(x) for x in args.ratios.split(",")]
else:
    if args.init > args.final:
        raise ValueError("Error: --init should be less than or equal to --final")
    ratio_list = []
    ratio = args.init
    while ratio <= args.final + 1e-8:
        ratio_list.append(round(ratio, 8))
        ratio += args.interval

# stretch list 
stretch_list = [x**(1/3) for x in ratio_list]
atoms, lc = my_read_vasp(os.path.join("y_full_relax_temp", "CONTCAR"))
atomsc = atoms.copy()
volume_ref = atomsc.get_volume()
films_stretch = stretch_list_along_direction_to_cell(atomsc, stretch_factor_list=stretch_list, stretch_direction_list=["xyz"])
print('=========================')
print('Ratio list   :', ratio_list)
print('Stretch list :', stretch_list)
print("Volume ref   :", volume_ref)
print(f"Total {len(films_stretch)} structures generated.\n")

workdir = os.path.join(myroot, "y_eos")
# existing output: ask before deleting (blank/No/no-tty aborts, --deleteold skips)
confirm_prepare_outdir(workdir, force=args.deleteold)
os.makedirs(os.path.join(workdir, "y_dir"), exist_ok=True)


for film, stretch, ratio in zip(films_stretch, stretch_list, ratio_list):
    dist_path = os.path.join(workdir, f"y_dir/{ratio:.8f}/")
    filmc = film.copy()
    volume = filmc.get_volume()

    print("==========================")
    print(f"Volume ratio : {ratio:.8f}")
    print(f"Stretch factor: {stretch:.8f}")
    print(f"Directory     : {dist_path}")

    # Check volume consistency
    if abs(volume - volume_ref * ratio) > 1e-10:
        raise ValueError(f"Error in volume calculation: calculated volume {volume}, expected volume {volume_ref * ratio}")

    os.makedirs(dist_path, exist_ok=True)
    cp_vaspfiles(cpfiles=[ "KPOINTS", "INCAR", "POTCAR", "Y_CONSTR_LATT", "CHGCAR"], oldir = "y_full_relax_temp", newdir = dist_path)

    os.chdir(dist_path)
    if args.isif == 2:
        os.system("pei_vasp_univ_find_and_change -isif 2 && echo -e '0\n0.0  0.0  0.0  0.0  0.0  0.0' > Y_CONSTR_LATT")
    elif args.isif == 4:
        os.system("pei_vasp_univ_find_and_change -isif 4 && echo -e '0\n0.0  0.0  0.0  0.0  0.0  0.0' > Y_CONSTR_LATT")
    
    # don't output CHGCAR
    os.system("pei_vasp_univ_find_and_change -lcharg F")

    os.chdir(myroot)

    my_write_vasp(os.path.join(dist_path, 'POSCAR'), filmc, lattice_scale_factor = lc * stretch, label = f"Stretched by {stretch:.8f}, ratio {ratio:.8f}. ")


shutil.rmtree("y_full_relax_temp")

print("\n==========================")
print("Please check the input files and POSCARs in each directory.")
print("Then you can run VASP (pei_vasp_univ_sbatch_jobs) in each directory to get the energies.")
print("After that, you can use 'pei_vasp_univ_plot_all -eos (TODO)' to get the plot.")
print('TODO: implement the -eos option in pei_vasp_univ_plot_all script.')
print("==========================")


