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
file_output = os.path.join(myroot, 'p_neb_py.dat')
print("# Directories found: " + " ".join(sorted(lfolder)))

lforce_b = []
lenergy = []
ldist_cum = []
ldist_cum_outcar = []

atoms_ref = read(os.path.join(myroot, lfolder[0], 'CONTCAR'), format='vasp')
nions = len(atoms_ref)
nions_more = nions + 1
print(f"# Number of ions   : {nions}")


dist_cum = 0.0
dist_cum_outcar = 0.0

# get lenergy, ldist_cum, lforce_b
for id, folder in enumerate(lfolder):
    path = os.path.join(myroot, folder)
    file_outcar = os.path.join(path, 'OUTCAR')

    print("\n===========================")
    print(f"Processing folder  : {folder}")

    # default read the last index -1 of OUTCAR
    if not os.path.exists(file_outcar):
        raise FileNotFoundError(f"No OUTCAR in {path}!")
    else:
        outcar = read(file_outcar, format='vasp-out')
    print(f"Reading OUTCAR from: {path}")
    #print(outcar)
    
    # get energy
    energy = outcar.get_total_energy()
    if folder == '00':
        energy0 = energy
    energy -= energy0
    lenergy.append(energy)
    print('Energy (eV):',energy)

    # get dist
    # 1. dist.pl file1 file2
    if folder == '00':
        dist = 0.0
    else:
        output_run = subprocess.run(["dist.pl", f"{lfolder[id-1]}/CONTCAR", f"{folder}/CONTCAR"],
                                capture_output=True,
                                 universal_newlines=True,  # 返回文本而不是字节
                                check=True)  
        dist = float(output_run.stdout)
    dist_cum += dist
    ldist_cum.append(dist_cum)
    print("Distance (Å):", dist_cum)
    # 2. grep form OUTCAR
    # to check 02, 03, ..., -2 , the dist_next of id i is equal to dist_pre of id i+1
    id_temp = 0
    dist_pre = 0
    dist_next = 0
    if folder == '00':
        dist_pre = 0.0
    elif folder == lfolder[-1]:
        dist_pre = dist_next_old
    else:
        cmd = f"grep 'NEB: distance' {file_outcar} | tail -1 "
        output_run = subprocess.run(cmd, shell=True,
                                        capture_output=True, 
                                        text=True,            # 必须，因为要处理字符串
                                        check=True)
        
        # dist_nex in lfolder[-2] will be transferred to dist_pre in next image lfolder[-1]
        dist_pre  = float(output_run.stdout.split()[-3])
        dist_next = float(output_run.stdout.split()[-2]) 

        if id_temp > 0:
            if abs(dist_pre - dist_next_old) > 1e-5:
                raise ValueError(f"Distance mismatch between images: {folder} dist_pre={dist_pre}, previous dist_next={dist_next_old}")
        
        id_temp += 1
    dist_pre_old = dist_pre
    dist_next_old = dist_next
        
    dist_cum_outcar += dist_pre
    ldist_cum_outcar.append(dist_cum_outcar)
    print("Distance from OUTCAR (Å):", dist_cum_outcar)

    # Check dist and dist_cum consistency
    # if abs(dist_pre - dist) > 1e-6:
    #     raise ValueError(f"Distance mismatch in {folder}: dist.pl={dist}, OUTCAR dist_pre={dist_pre}")

    # if abs(dist_cum - dist_cum_outcar) > 1e-5:
    #     raise ValueError(f"Cumulative distance mismatch in {folder}: dist.pl={dist_cum}, OUTCAR={dist_cum_outcar}")


    # get force_b
    #$force = `grep 'NEB: projections' $directories[$i]/OUTCAR|tail -1`;
    if folder in [lfolder[0], lfolder[-1]]:
        force_b = 0.0
    else:
        cmd = f"grep 'NEB: projections' {file_outcar} | tail -1 "
        output_run = subprocess.run(cmd, shell=True,
                                        capture_output=True,
                                        text=True,            # 必须，因为要处理字符串
                                        check=True)
        force_b = float(output_run.stdout.split()[-1])
    lforce_b.append(force_b)
    print("Force_b (eV/Å):", force_b)
        

with open(file_output, 'w') as f:
    f.write("%3s %16s %16s %16s %8s\n" % ("id", "dist_cum(Å)", "energy_rel(eV)", "force_b(eV/Å)", "image"))
    for dist_cum, energy, force_b, folder in zip(ldist_cum, lenergy, lforce_b, lfolder):
        f.write(f"{int(folder):3} {dist_cum:16.8f} {energy:16.8f} {force_b:16.8f} {int(folder):3}\n")

print("\n\n"+"="*100)
os.system("cat -n " + file_output)
print("="*100)