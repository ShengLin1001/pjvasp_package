# vasp_utils 

This repo contains scripts for various VASP workflows, including job submission (with slurm) and post-analyzing (with python3). 



# Setup

1. clone the repos and add them to your PATH:     
   https://github.com/BinglunYin/vasp_utils    
   https://github.com/BinglunYin/slurm_utils    

1. create links to `python3` and `pip3` at `$HOME/opt/bin/` (do not recommended).

1. Recommended you use python virtual environments. Befor you do post job, source activate the environment.

1. pip install a package:     
   ```shell
   pip3 install --upgrade --user   git+https://github.com/BinglunYin/myalloy_package.git 
   ```
   or download the zip file and 
   ```shell
   pip3 install --user -e  ./   
   ```






# Usage 1: run workflows based on `y_full_relax` 

1. Calculate the reference state in the folder `y_full_relax`. This folder will serve as the basis for the following workflows.  
   ```shell
   project_folder
   ├── y_full_relax
   │   ├── INCAR
   │   ├── KPOINTS
   │   ├── POSCAR
   │   ├── POTCAR
   │   ├── sub.vasp
   ```



1. For example. If you would like to calculate the EOS of the reference state, run the following command at the level of `y_full_relax`:    
   ```shell
   yin_vasp_run_eos
   ```

   When all the jobs finish, run the command:    
   ```shell
   yin_vasp_plot_all  -eos 
   ```

   Then you will have the EOS result.   
   
1. All the workflow commands are with the name `yin_vasp_run_*`.
   



# Usage 2: run a batch of POSCAR files 


1. In your `project_folder`, create `y_src` to contain all the source files. These files will not be changed during calculation. 

   ```shell
   project_folder
   ├── y_src
   │   ├── INCAR
   │   ├── KPOINTS
   │   ├── POTCAR
   │   ├── sub.vasp
   │   ├── poscars2
   │   │   ├── POSCAR_001 
   │   │   ├── POSCAR_002 
   │   │   ├── POSCAR_003 
   ```

   
1. To submit jobs, run the following command at the level of `y_src`.   
    
   ```shell
   yin_vasp_univ_sub_poscars2
   ```

1. To monitor the job status, run

   ```shell
   yin_vasp_univ_post
   ```

1. When all the jobs finish, remove the failed jobs from `y_dir`, e.g., aborted, non-converged, etc, and then post-analyze.

   ```shell
   yin_vasp_univ_post  -v  
   ```

   Make sure all the remaining jobs in `y_dir` are successfully completed.


# Data management: clean, compress, and store

Inside the `project_folder`, run

   ```shell
   yin_vasp_univ_clean_up
   cd ..
   mv  project_folder  yyyymmdd_project_folder
   nohup  ~/yin_github/linux_config/mktar2  -xz  yyyymmdd_project_folder/  &
   ```

When the compression is successfully completed (check log), move them to `/store`.



