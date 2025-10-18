import numpy as np
import argparse
import os,shutil

os.system('rm -r work_dir/*')

# input.nn.no without SFs
SF_all = open("SFs_all.dat","a")
with open('input.nn.all','r') as r:
    lines = r.readlines()
with open('input.nn.no','w') as w:
    for l in lines:
        if l[0] == '#':
            w.write(l)
        elif 'symfunction_short' not in l:
            w.write(l)
        else:
            SF_all.write(l)   # add to SF_all 
SF_all.close()
            
# append two 0 value SF, ortherwise n2p2 will break
with open('input.nn.no','a') as w:
    w.write('symfunction_short Mg 3 Mg Mg 1  1.000 1.000 0.100\n')

# input.nn.all contian all SFs
# write input.nn.no without SFs
# SFs_all.dat 
i=0
with open('input.nn.all','r') as r:
    lines = r.readlines()
for l in lines:
    if l[0] == '#':
        pass
    elif 'symfunction_short' not in l:
        pass
    else:
        shutil.copyfile('input.nn.no','input.nn')   # sub
        with open('input.nn','a') as w:
            w.write(l) 
        os.system('cp input.data scaling/')
        os.system('mv input.nn scaling/')
        os.system('cp -r scaling work_dir/%04d'%i)
        os.chdir('work_dir/%04d'%i)
        os.system('qsub sub.pbs')
        os.chdir('../../')
        i += 1



