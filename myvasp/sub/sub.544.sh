#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 12
source /public3/soft/modules/module.sh
module load mpi/intel/17.0.7-thc
export PATH=/public3/home/scg6928/mysoft/vasp/544-opt/vasp.5.4.4.pl2/bin:$PATH
srun vasp_std_544_opt
