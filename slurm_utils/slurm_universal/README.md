# This `slurm_universal` directory is designed to put some well-defined sub.* script

1. Submit a single VASP job in a specified subdirectory.

        `sub.544.sh` → `pei_vasp_univ_sbatch`

2. Submit VASP jobs sequentially for all subdirectories under `y_dir`.

        `sub.544.sequential.sh` → `pei_vasp_univ_sbatch_sequential` → `pei_vasp_univ_sbatch`

3.  Submit VASP jobs in parallel for all subdirectories under `y_dir` with a specified number of parallel jobs.

        `sub.544.parallel.sh` → `pei_vasp_univ_sbatch_chunk` → `pei_vasp_univ_sbatch_sequential` →
        `pei_vasp_univ_sbatch`

4. Submit VASP jobs in parallel for all subdirectories under `y_dir`.

        `sub.544.parallel.sh` ← `pei_vasp_univ_sbatch` ← `pei_vasp_univ_sbatch_parallel`

Note: All preparatory steps, including running `pei_vasp_univ_clean_up_full`, modifying the `INCAR` settings, and updating the `KPOINTS` settings, must be completed before running the scripts listed above. These scripts are only responsible for submitting jobs to SLURM and will not modify any input files.
