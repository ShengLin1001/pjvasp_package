import numpy as np
from myvasp import vasp_func as vf
import copy, os, sys, shutil, time
import StructureForge.InputTools.VASPReader as vr

#folders contain all struct of testing set
data = "/store/tliu/dataset_nnp/0101/vasp/"
folders = [data + f + '/' for f in os.listdir(data) if os.path.isdir(data + f)]


energy = []
forces_cfg = []
folders.sort()
for s in folders:
    #for nnp and run lammps in 
    atoms = vf.my_read_vasp(s + 'CONTCAR')
    # read vasp        reference value
    ref = vr.reader(s)
    n = len(ref.atoms)
    # write forces
    forces = [] 
    for j in range(len(ref.atoms)):
        forces.append([ref.atoms[j].Force[0],ref.atoms[j].Force[1],ref.atoms[j].Force[2]]) 
    forces_cfg.append(np.array(forces)) 
    # write energy
    energy.append(ref.energy)
    


energy = np.array(energy)
forces = np.array(forces_cfg)

np.save('energy.npy',energy)
np.save('forces.npy',forces)


