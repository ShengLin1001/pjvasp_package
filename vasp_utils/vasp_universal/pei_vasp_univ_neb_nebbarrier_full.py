#!/public3/home/scg6928/mysoft/env/pyenv/dft/bin/python

import os
import glob
import subprocess
from ase.io import read

# get subfoler list
myroot = os.getcwd()
lfolder = glob.glob('[0-9][0-9]')
lfolder = [f for f in lfolder if os.path.isdir(f)]
lfolder.sort()
file_output = os.path.join(myroot, 'p_neb_full.dat')
print("# Directories found: " + " ".join(sorted(lfolder)))

# For full, list of list, each image, each frame
lforce_bs = []
# here is energies relative to initial state
lenergies = []
# absolute energies
lenergies_absolute = []
ldist_cums = []
# list of int
lframes = []

atoms_ref = read(os.path.join(myroot, lfolder[0], 'CONTCAR'), format='vasp')
nions = len(atoms_ref)
nions_more = nions + 1
print(f"# Number of ions   : {nions}")

dist_cum = 0.0

# get lenergy, ldist_cum, lforce_b
for id, folder in enumerate(lfolder):
    path = os.path.join(myroot, folder)
    file_outcar = os.path.join(path, 'OUTCAR')

    print("\n===========================")
    print(f"Processing folder  : {folder}")

    # read the full index of OUTCAR
    if not os.path.exists(file_outcar):
        raise FileNotFoundError(f"No OUTCAR in {path}!")
    else:
        outcar = read(file_outcar, format='vasp-out', index=':')

    nframes = len(outcar)
    print(f"Nframes            :", nframes)
    lframes.append(nframes)
    force_bs = []
    energies = []
    energies_absolute = []
    dist_cums = []
    
    # append energies
    if folder == '00':
        energy0 = outcar[-1].get_total_energy()

    if folder in [lfolder[0], lfolder[-1]]:
        energy = outcar[-1].get_total_energy()
        energies_absolute.append(energy)
        energy -= energy0
        energies.append(energy)
    else:
        for iframe in range(nframes):
            energy = outcar[iframe].get_total_energy()
            energies_absolute.append(energy)
            energy -= energy0
            energies.append(energy)
    lenergies.append(energies)
    lenergies_absolute.append(energies_absolute)
    print('Energies (eV)      :', energies)

    # append dist_cums
    # # 1. dist.pl file1 file2
    # # 2. grep form OUTCAR
    # # to check 02, 03, ..., -2 , the dist_next of id i is equal to dist_pre of id i+1
    #id_temp = 0
    dist_pre = 0
    dist_next = 0
    # First to '00' loop: [[0.0]]
    # Then to 
    if folder in '00':
        dist_pre = 0.0
        dist_cum += dist_pre
        dist_cums.append(dist_cum)
    elif folder == lfolder[-1]:
        for id_frame in range(len(ldist_cums[id-1])):
            dist_cum_old = ldist_cums[id-1][id_frame]
            dist_pre = ldist_next[id_frame]
            dist_cum_temp = dist_cum_old + dist_pre
            dist_cums.append(dist_cum_temp)
    else:
        cmd = f"grep 'NEB: distance' {file_outcar} "
        output_run = subprocess.run(cmd, shell=True,
                                        capture_output=True, 
                                        text=True,            # 必须，因为要处理字符串
                                        check=True)
        
        # ldist_next will transfer to next image lfolder[-1]
        lines = output_run.stdout.strip().split('\n')
        ldist_next = []
        # nframes = len(lines)
        for id_frame, line in enumerate(lines):
            if len(ldist_cums[id-1]) == 1:
                dist_cum_old = ldist_cums[id-1][0]
            else:
                dist_cum_old = ldist_cums[id-1][id_frame]

            dist_pre  = float(line.split()[-3])
            dist_next = float(line.split()[-2]) 
            ldist_next.append(dist_next)

            dist_cum_temp = dist_cum_old + dist_pre
            dist_cums.append(dist_cum_temp)

    ldist_cums.append(dist_cums)
    print('Dist cums (A)      :', dist_cums)

    # get force_bs
    if folder in [lfolder[0], lfolder[-1]]:
        force_b = 0.0
        force_bs.append(force_b)
    else:
        cmd = f"grep 'NEB: projections' {file_outcar} "
        output_run = subprocess.run(cmd, shell=True,
                                        capture_output=True,
                                        text=True,            # 必须，因为要处理字符串
                                        check=True)
        lines = output_run.stdout.strip().split('\n')
        for id_frame, line in enumerate(lines):
            line = lines[id_frame]
            force_b = float(line.split()[-1])
            force_bs.append(force_b)

    lforce_bs.append(force_bs)
    print('Force_bs (eV/A)    :', force_bs)


# adjust to same length
print("\nAdjusting to same length for all images...")
max_nframes = max(lframes)

if len(set(lframes[1:-1])) != 1:
    raise ValueError("Middle images have different number of frames!")

print("Lframes       :", lframes)
print("Max nframes   :", max_nframes)
# For lenergies
# In line 55-58, already maked sure the 0, and -1 index only have the last energies
lenergies[0] = lenergies[0] * max_nframes
lenergies[-1] = lenergies[-1] * max_nframes
lenergies_absolute[0] = lenergies_absolute[0] * max_nframes
lenergies_absolute[-1] = lenergies_absolute[-1] * max_nframes
# For ldist_cums
# In line 76-79, already maked sure the 0 index only have the [0.0] dist_cum
# In line 80-85, already maked sure the -1 index have the dist_cums according to max_nframes
ldist_cums[0] = ldist_cums[0] * max_nframes
# For lforce_bs
# In line 115-117, already maked sure the 0, and -1 index only have the [0.0] force_bs
lforce_bs[0] = lforce_bs[0] * max_nframes
lforce_bs[-1] = lforce_bs[-1] * max_nframes

# write to file
with open(file_output, 'w') as f:
    f.write("%3s %8s %16s %16s %16s %5s %16s\n" % ("id", "frame", "dist_cum(Å)", "energy_rel(eV)", "force_b(eV/Å)", "natoms", "energy_abs(eV)"))
    for frame in range(max_nframes):
        for id, folder in enumerate(lfolder):
            dist_cum = ldist_cums[id][frame]
            energy = lenergies[id][frame]
            energy_abs = lenergies_absolute[id][frame]
            force_b = lforce_bs[id][frame]
            f.write(f"{int(folder):3} {frame:8} {dist_cum:16.8f} {energy:16.8f} {force_b:16.8f} {int(nions): 5} {energy_abs:16.8f}\n")