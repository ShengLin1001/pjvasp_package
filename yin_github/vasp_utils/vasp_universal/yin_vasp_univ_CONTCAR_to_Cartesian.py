#!/home/yin/opt/bin/python3

import numpy as np
from ase.io.vasp import read_vasp, write_vasp
from ase.constraints import FixedLine
import sys, os


filename1 = 'CONTCAR'    # use ./CONTCAR leads to error sometimes


f = open(filename1,'r')
a0 = float( f.readlines()[1] )
f.close()
# print('==> a0 = {0}'.format(a0))

ASE_Atoms = read_vasp(filename1)
atoms_pos = ASE_Atoms.positions
natoms = atoms_pos.shape[0]
# print('==> natoms = {0}'.format(natoms))

ASE_Atoms.set_cell( ASE_Atoms.cell[:]/a0 )
ASE_Atoms.set_positions( atoms_pos/a0, apply_constraint=False )  # fuking important!


if len(sys.argv)  == 1:        
    filename2 = 'POSCAR_Cartesian'

elif len(sys.argv)  == 2:
    if sys.argv[1] == '-FFT':
        filename2 = 'POSCAR_Cartesian_FFT'
        constrained_atoms = np.arange(natoms)  # apply constraint to all atoms
        relax_direction = ASE_Atoms.cell[2,:]          # relax only in the a3 direction
        constraint_a3only = FixedLine(constrained_atoms, relax_direction)
        ASE_Atoms.set_constraint(constraint_a3only)


write_vasp(filename2, ASE_Atoms,
label='system_name', direct=False)


with open(filename2) as f:
    lines = f.readlines()

lines[1] = ' %.16f \n' % (a0)

with open('file_temp', "w") as f:
    f.writelines(lines)

os.replace('file_temp', filename2)



# check
ASE_atoms1 = read_vasp(filename1)
ASE_atoms2 = read_vasp(filename2)

temp = np.linalg.norm( ASE_atoms2.cell[:] - ASE_atoms1.cell[:] )
if temp>1e-10:
    sys.exit('==> Abort! wrong cell, with norm = {0}'.format(temp))


temp = np.linalg.norm( ASE_atoms2.positions - ASE_atoms1.positions )
if temp>1e-10:
    sys.exit('==> Abort! wrong positions, with norm = {0}'.format(temp))

