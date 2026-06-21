#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 1        

# skip this jobid
srun pei_vasp_univ_monitor_error "$SLURM_JOB_ID"