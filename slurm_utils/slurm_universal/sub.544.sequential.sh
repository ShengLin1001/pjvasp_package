#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 128         
source /public3/soft/modules/module.sh
module load mpi/intel/17.0.7-thc
export PATH=/public3/home/scg6928/mysoft/vasp/vasp/544-yin/vasp.5.4.4.pl2/bin:$PATH

echo "==> sbatch all jobs under y_dir sequentially"
for dir in ./y_dir/*/; do
    cd $dir
    echo ================ $dir
    echo 'submit dir:' `pwd`
    #sbatch sub.*
    srun vasp_std
    cd - > /dev/null
done
echo "done!"


