# `slurm_universal` 目录用于存放一些定义明确的 `sub.*` 脚本

注意：在运行下面列出的脚本之前，必须先完成所有前处理步骤，包括运行 `pei_vasp_univ_clean_up_full`、修改 `INCAR` 设置以及更新 `KPOINTS` 设置。

只有 `pei_vasp_univ_sbatch` 会删除旧的输出文件，但会保留 `*.out` 文件。

```bash
rm -f CHG CONTCAR* DOSCAR* DYNMAT EIGENVAL IBZKPT OPTIC OSZICAR* OUTCAR* PROCAR* PCDAT W* XDATCAR* PARCHG* vasprun.xml SUMMARY.* REPORT
```

## 顺序提交工作流

这里包含两种顺序提交工作流。应根据你的需求选择：是希望 VASP 在同一个 Slurm 资源分配中依次运行，还是希望每个子目录分别提交自己的子 Slurm 作业。

### 所有子目录共用一个 Slurm 资源分配

当你希望只排队一次，并在同一个 Slurm 资源分配中逐个运行尚未完成的计算时，使用 `pei_slurm_univ_vasp_544_sequential_single_allocation.sh`。

在包含 `y_dir` 的目录中运行：

```bash
sbatch pei_slurm_univ_vasp_544_sequential_single_allocation.sh
```

该工作流为：

```bash
pei_slurm_univ_vasp_544_sequential_single_allocation.sh -> pei_vasp_univ_sbatch_sequential_single_allocation -> pei_vasp_univ_sbatch
```

该模式使用 `pei_slurm_univ_vasp_544_sequential_single_allocation.sh` 所申请的计算资源，通过 `srun` 直接运行 VASP，并在提交目录中写入一个父 Slurm 日志。当排队时间较长，并且希望用同一个资源分配处理所有选定子目录时，这种模式比较有用。坏处是一旦有一个任务报错，monitor 脚本会直接把整个父作业终止，无法继续处理剩余的子目录。

如果只想运行指定的子目录，可以传入一个由逗号分隔的 `lsubdir` 列表。列表中的名称应为 `y_dir` 下的子目录 basename：

```bash
sbatch pei_slurm_univ_vasp_544_sequential_single_allocation.sh "0.98,1.00,1.02"
```

### 每个子目录对应一个子 Slurm 作业

当你希望一个父 Slurm 作业依次遍历 `y_dir` 下的子目录，而每个计算仍然作为独立的子 Slurm 作业运行时，使用 `pei_slurm_univ_vasp_544_sequential_each_subdir_sbatch.sh`。

在包含 `y_dir` 的目录中运行：

```bash
sbatch pei_slurm_univ_vasp_544_sequential_each_subdir_sbatch.sh
```

该工作流为：

```bash
pei_slurm_univ_vasp_544_sequential_each_subdir_sbatch.sh -> pei_vasp_univ_sbatch_sequential_each_subdir_sbatch -> sbatch --wait pei_slurm_univ_vasp_544.sh
```

对于 `y_dir` 下的每个一级子目录，该辅助脚本会：

1. 进入该子目录。
2. 如果 `OUTCAR` 中已经包含收敛标记，则跳过该子目录。
3. 在该子目录内部运行 `sbatch --wait pei_slurm_univ_vasp_544.sh`。
4. 子 Slurm 作业结束后，再次检查 `OUTCAR` 是否收敛。
5. 将未收敛的目录记录到父 Slurm 日志中，然后继续处理下一个子目录。

这种模式会生成两级 Slurm 日志：

- 父脚本 `pei_slurm_univ_vasp_544_sequential_each_subdir_sbatch.sh` 的日志会写入提交顺序作业时所在的目录。它记录目录进入、作业提交状态以及收敛检查结果。
- 每个子脚本 `pei_slurm_univ_vasp_544.sh` 的日志会写入对应的计算子目录中，因为 `sbatch` 是在进入该子目录之后调用的。监控工作流可以检查这个子日志，并在检测到错误信息时终止对应的子作业。

如果只想运行指定的子目录，可以传入一个由逗号分隔的 `lsubdir` 列表。列表中的名称应为 `y_dir` 下的子目录 basename：

```bash
sbatch pei_slurm_univ_vasp_544_sequential_each_subdir_sbatch.sh "0.98,1.00,1.02"
```

## 其他工作流

1. 在一个计算目录中提交一个 VASP 作业。

```bash
pei_slurm_univ_vasp_544.sh -> pei_vasp_univ_sbatch
```

   在某个计算目录中运行：

```bash
sbatch pei_slurm_univ_vasp_544.sh
```

2. 将所有未完成的子目录作为相互独立的 Slurm 作业提交，不等待每个作业完成。

```bash
pei_vasp_univ_sbatch_parallel -> sbatch pei_slurm_univ_vasp_544.sh
```

   在包含 `y_dir` 的目录中运行：

```bash
pei_vasp_univ_sbatch_parallel ./y_dir pei_slurm_univ_vasp_544.sh
```

3. 通过顺序提交包装脚本提交 `y_dir` 中选定的分块。

```bash
pei_vasp_univ_sbatch_chunk -> pei_slurm_univ_vasp_544_sequential_single_allocation.sh 
pei_vasp_univ_sbatch_chunk -> pei_slurm_univ_vasp_544_sequential_each_subdir_sbatch.sh
```

   当你有意创建多个父顺序作业，并希望每个父作业负责 `y_dir` 的一个子集时，可以使用这种模式。提交脚本参数决定每个分块使用哪一种顺序提交模式。

```bash
pei_vasp_univ_sbatch_chunk <num_chunks> <submit_script> [--dry-run]
pei_vasp_univ_sbatch_chunk 10 pei_slurm_univ_vasp_544_sequential_single_allocation.sh
pei_vasp_univ_sbatch_chunk 10 pei_slurm_univ_vasp_544_sequential_each_subdir_sbatch.sh
```
