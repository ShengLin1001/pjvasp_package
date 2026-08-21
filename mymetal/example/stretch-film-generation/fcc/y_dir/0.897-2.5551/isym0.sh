#!/bin/bash
bash cp_out.sh -o slurm*.out -t 1/ 
rm slurm*.out
rm POSCAR CONTCAR
cp 1/CONTCAR POSCAR
bash sed_incar.sh ISYM 0
bash sed_incar.sh ISTART 1
sbatch sub*
