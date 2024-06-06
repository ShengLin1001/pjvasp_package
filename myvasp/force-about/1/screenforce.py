#!/usr/bin/python -u 
# -*- coding: utf-8 -*-
# Calculate the force of each atom. Screen the force not fullfil convergence criteria and the 
# corresponding atom number.
# Written by qmlearner at Homestay Uni. 


#import math 
#import numpy as np
#import re
import os

a=open('sep_forcevectors.dat','r+')

af='atomicforce.dat'

if os.path.exists(af):
    print("---------------------------------------------")
    print("%s file already exists. It will be deleted." %(af)) 
    print("---------------------------------------------")
    os.remove(af)

b=open(af,'a')


la=a.readlines()

###############################################
#Calculate the force performing on each atom###
###############################################
for i in range(len(la)):
    x=float(la[i].strip().split()[0])
    y=float(la[i].strip().split()[1])
    z=float(la[i].strip().split()[2])
   
    F=(x**2 + y**2 + z**2)**0.5
    b.write(str(F))
    b.write('\n')

a.close()
b.close()

################################################
#Output the force that does not fulfill        #
#convergence criteria.                         #
################################################

c=open('atomicforce.dat','r')

ce='lt_ediffg.dat'

if os.path.exists(ce):
    print("------------------------------------------")
    print("%s file already exists. It will be deleted." %(ce)) 
    print("------------------------------------------")
    os.remove(ce)

d=open(ce,'a')


lc=c.readlines()

##################################################
#Change the ediffg according to your own settings#
##################################################
ediffg=0.03

for i in range(len(lc)):
    F=float(lc[i].strip())       
    if F >= ediffg:
        d.write(str(i+1))
        d.write('\t')
        d.write(str(F))
        d.write('\n')

c.close()
d.close()
