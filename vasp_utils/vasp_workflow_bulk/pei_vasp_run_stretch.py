#!/public3/home/scg6928/mysoft/env/pyenv/dft/bin/python
# J. Pei, 2025-09-21
import os
import numpy as np
import argparse
import shutil
from mymetal.build.film.stretch import stretch_list_along_direction_to_cell, adjust_lattice_for_volume_conservation
from mymetal.build.workflow.general import cp_vaspfiles, rm_i, compare_three_lattices
from mymetal.io.vasp import my_write_vasp, my_read_vasp

# parse arguments
parser = argparse.ArgumentParser(description="Stretch unit cells to find the equilibrium lattice constants.")
parser.add_argument("--type", type = str, choices=["xyz", "xy", "xz", "yz", "x", "y", "z"], required=True, help="Stretch direction")
parser.add_argument("--init", type = float, default = -4/1000,      help = "Initial strain  (e.g., -4/1000 for -0.4%%)")
parser.add_argument("--final", type = float, default = 4/1000,      help = "Final   strain  (e.g.,  4/1000 for  0.4%%)")
parser.add_argument("--interval", type = float, default = 0.5/1000, help = "Strain interval (e.g., 0.5/1000 for 0.05%%)")
# other methods to set up stretch list
parser.add_argument("--strains", type = str, default = None, help = "Comma-separated strain list (e.g., -0.004,-0.002,0.0,0.002,0.004)")
parser.add_argument("--keepvolume", action = "store_true", help = "Adjust the unstretched directions to keep the volume unchanged")
parser.add_argument("--deleteold", action = "store_true", help="Whether to delete existing directory automatically.")
args = parser.parse_args()

# deformed 'x' => unstretched 'y' and 'z' => [1, 2]
axis_map = {"x": 0, "y": 1, "z": 2}
all_axes = set("xyz")
unstretched_axes = all_axes - set(args.type)
change_axis = [axis_map[ax] for ax in unstretched_axes]

print("Stretch type    :", args.type)
print("Initial strain  :", args.init)
print("Final   strain  :", args.final)
print("Strain  interval:", args.interval)
print("Strain  list    :", args.strains)
print("Keep volume     :", args.keepvolume)

myroot = os.getcwd()
os.makedirs("y_full_relax_temp", exist_ok=True)
# copy y_full_relax to y_full_relax_temp
cp_vaspfiles(oldir = "y_full_relax", newdir = "y_full_relax_temp")

# stretch list
if args.strains is not None:
    strain_list = [float(x) for x in args.strains.split(",")]
    
else:
    if args.init > args.final:
        raise ValueError("Error: --init should be less than or equal to --final")
    strain_list = []
    strain = args.init
    while strain <= args.final + 1e-8:
        strain_list.append(round(strain, 8))
        strain += args.interval

stretch_list = [1 + x for x in strain_list]
atoms, lc = my_read_vasp(os.path.join("y_full_relax_temp", "CONTCAR"))
atomsc = atoms.copy()
films_stretch = stretch_list_along_direction_to_cell(atomsc, stretch_factor_list=stretch_list, stretch_direction_list=[args.type])
print('=========================')
print('Stretch list:', stretch_list)
print(f"Total {len(films_stretch)} structures generated.")

workdir = os.path.join(myroot, "y_stretch")
# Check if the directory already exists
#rm_i(workdir = workdir)
if args.deleteold:
    if os.path.exists(workdir):
        print(f"The directory {workdir} already exists. Deleting it automatically as per --deleteold y.")
        shutil.rmtree(workdir)

os.makedirs(os.path.join(workdir, "y_dir"), exist_ok=True)

for film, stretch in zip(films_stretch, stretch_list):
    dist_path = os.path.join(workdir, f"y_dir/{stretch:.8f}/")
    filmc = film.copy()
    os.makedirs(dist_path, exist_ok=True)
    cp_vaspfiles(cpfiles=[ "KPOINTS", "INCAR", "POTCAR", "Y_CONSTR_LATT"], oldir = "y_full_relax_temp", newdir = dist_path)
    print("==========================")
    print(f"Stretch factor: {stretch:.8f}")
    print(f"Directory     : {dist_path}")

    lattice_before = np.array(atomsc.get_cell())
    lattice_after  = np.array(film.get_cell())
    temp, _ = adjust_lattice_for_volume_conservation(np.array(atomsc.get_cell()), np.array(film.get_cell()), change_axis = change_axis)
    lattice_after_adjusted = np.array(temp)
    compare_three_lattices(lattice_before, lattice_after, lattice_after_adjusted)
    if args.keepvolume:
        print("Adjusting the unstretched directions to keep the volume unchanged.")
        filmc.set_cell(lattice_after_adjusted, scale_atoms=True)

    my_write_vasp(os.path.join(dist_path, 'POSCAR'), filmc, lattice_scale_factor = lc * stretch, label = f"Stretched by {stretch:.8f}. ")

shutil.rmtree("y_full_relax_temp")

print("==========================")
print("Please check the input files and POSCARs in each directory.")
print("Then you can run VASP (pei_vasp_univ_sbatch_jobs) in each directory to get the energies.")
print("After that, you can use 'pei_vasp_univ_plot_all -stretch' to get the plot.")
print("==========================")
