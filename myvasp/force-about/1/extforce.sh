#!/bin/bash
# Written by qmlearner at Homestay Uni.

######################################################
#Define the number of atoms (N_atoms) of your system.#
######################################################
N_atoms=249
N_Lines_1=$(echo "${N_atoms}+1" | bc)
N_Lines_2=$(echo "${N_atoms}+2" | bc)

########################################
#Extract the force vectors of each atom# 
#for each ionic step.                  #
########################################
if [ ! -f "OUTCAR" ]; then
 echo "Does not find OUTCAR."
else
 sed -n "/TOTAL-FORCE/,+${N_Lines_1}p" OUTCAR > forcevectors.dat
fi

N_ionic_steps=$(grep -o 'TOTAL-FORCE' forcevectors.dat | wc -l)

cp forcevectors.dat forcevectors-back.dat

for ((i=1;i<=${N_ionic_steps};i++));
do
  sed -n "3,${N_Lines_2}p" forcevectors.dat | awk '{printf "%10s   %10s   %10s   \n",$4,$5,$6}' > sep_forcevectors.dat
  sed  "1,${N_Lines_2}d" forcevectors.dat > forcevectors_del.dat
  mv forcevectors_del.dat forcevectors.dat
  python screenforce.py
  mv lt_ediffg.dat lt_ediffg_${i}.dat
done

rm forcevectors.dat sep_forcevectors.dat atomicforce.dat
