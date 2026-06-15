# `slurm_universal` 通用提交脚本

本目录存放定义明确的 `sub.*` / `pei_slurm_univ_*` 提交脚本。

> 注意：运行下面的脚本前，必须先完成所有前处理，包括 `pei_vasp_univ_clean_up_full`、修改 `INCAR`、更新 `KPOINTS`。
> 只有 `pei_vasp_univ_sbatch` 会删除旧的输出文件（保留 `*.out`）：
>
> ```bash
> rm -f CHG CONTCAR* DOSCAR* DYNMAT EIGENVAL IBZKPT OPTIC OSZICAR* OUTCAR* PROCAR* PCDAT W* XDATCAR* PARCHG* vasprun.xml SUMMARY.* REPORT
> ```

## 架构

整套批量提交围绕**一个软件无关的引擎** `pei_slurm_univ_submit`：它遍历 `root_dir` 下的一级子目录，对每个目录执行由 `--mode` 选择的动作，并打印汇总。VASP 相关脚本（父 `.sh` 与运行 helper）都构建在它之上，直接调用引擎。

```
pei_slurm_univ_submit --mode parallel|each-subdir|single-alloc [options]
    --root-dir DIR               # 作业根目录（默认 ./y_dir）
    --lsubdir a,b,c              # 只跑这些一级子目录（默认：全部）
    --submit-script SCRIPT       # parallel / each-subdir：每个目录里 sbatch 的脚本
    --run-cmd "CMD ..."          # single-alloc：在每个目录内直接运行的命令
    --skip-if-file-contains F M  # 当目录内文件 F 含标记 M 时跳过
    --chunks K --chunk-parent S  # 把目录切成 K 块，每块 sbatch 一个父脚本 S
    --dry-run                    # 只打印动作，不真正提交/运行
```

三种 `--mode`（每个目录的动作）：

| `--mode` | 每个目录做什么 | 运行位置 | 父作业占核 |
|---|---|---|---|
| `parallel` | `sbatch SCRIPT`，不等待 | 登录节点 | 子作业各自申请 |
| `each-subdir` | `sbatch --wait SCRIPT`，结束后复检标记 | `-n 1` 父分配内 | 1 核父 + 子作业 |
| `single-alloc` | 在目录内直接 `run-cmd`（退出码 0 完成 / 10 跳过 / 其他 失败） | `-n 128` 父分配内 | 父持有全部核 |

**chunk 是正交选项**：`--chunks K` 把目录切成 K 块、每块 `sbatch` 一个父作业，用来把总核数预算铺开成 **K 条并发流**。它只对两个 sequential 模式有意义（默认串行的工作借此并发提速）；对 `parallel` 没意义（parallel 本就最大并发）。

## VASP 工作流

### 1. 单目录提交

在某个计算目录中运行：

```bash
sbatch pei_slurm_univ_vasp_544.sh
```

调用链：`pei_slurm_univ_vasp_544.sh -> pei_vasp_univ_sbatch`

### 2. 顺序：同一个 Slurm 资源分配（single-alloc）

只排队一次，在同一分配里逐个 `srun` 运行未完成的计算。在包含 `y_dir` 的目录中：

```bash
sbatch pei_slurm_univ_vasp_544_sequential_single_allocation.sh
# 只跑指定子目录：
sbatch pei_slurm_univ_vasp_544_sequential_single_allocation.sh "0.98,1.00,1.02"
```

调用链：`...single_allocation.sh -> pei_slurm_univ_submit --mode single-alloc -> pei_vasp_univ_sbatch`

坏处：一旦某个任务报错，monitor 脚本会终止整个父作业，无法继续剩余子目录。

### 3. 顺序：每个子目录一个子 Slurm 作业（each-subdir）

父作业（仅 `-n 1`）依次遍历子目录，每个计算作为独立子作业 `sbatch --wait` 运行：

```bash
sbatch pei_slurm_univ_vasp_544_sequential_each_subdir_sbatch.sh
sbatch pei_slurm_univ_vasp_544_sequential_each_subdir_sbatch.sh "0.98,1.00,1.02"
```

调用链：`...each_subdir_sbatch.sh -> pei_slurm_univ_submit --mode each-subdir -> sbatch --wait pei_slurm_univ_vasp_544.sh`

这种模式生成两级日志：父脚本日志写在提交目录，记录目录进入/提交状态/收敛复检；每个子作业 `pei_slurm_univ_vasp_544.sh` 的日志写在对应计算子目录。

### 4. 并行：互相独立、不等待

把所有未完成子目录作为独立 Slurm 作业提交：

```bash
pei_slurm_univ_submit --mode parallel --root-dir ./y_dir \
    --submit-script pei_slurm_univ_vasp_544.sh \
    --skip-if-file-contains OUTCAR "reached required accuracy"
```

调用链：`pei_slurm_univ_submit --mode parallel -> sbatch pei_slurm_univ_vasp_544.sh`

用 `--skip-if-file-contains OUTCAR "reached required accuracy"` 跳过已完成目录；加 `--dry-run` 预览、`--lsubdir "a,b,c"` 只跑部分。

### 5. 分块（chunk）：把核预算铺开成多条并发的顺序流

把 `y_dir` 切成 N 块，每块提交一个父顺序作业。`--chunk-parent` 决定每块用哪种 sequential 模式：

```bash
pei_slurm_univ_submit --chunks 10 --root-dir ./y_dir \
    --chunk-parent pei_slurm_univ_vasp_544_sequential_single_allocation.sh
pei_slurm_univ_submit --chunks 10 --root-dir ./y_dir \
    --chunk-parent pei_slurm_univ_vasp_544_sequential_each_subdir_sbatch.sh
```

加 `--dry-run` 只打印切分结果而不提交。

## 软件无关的入口

引擎不绑定 VASP，可用于 n2p2 / LAMMPS 等任何 `sbatch` 脚本——只需把 `--submit-script` 换成对应的 Slurm 脚本：

```bash
pei_slurm_univ_submit --mode parallel --root-dir ./train_dirs \
    --submit-script pei_slurm_univ_n2p2_train.sh --dry-run
pei_slurm_univ_submit --mode parallel --root-dir ./y_dir \
    --submit-script pei_slurm_univ_vasp_544.sh \
    --skip-if-file-contains OUTCAR "reached required accuracy"
```

## 文件一览

| 文件 | 角色 |
|---|---|
| `pei_slurm_univ_lib.sh` | sourced 共享库：输出助手 + 目录解析 + 脚本路径解析 |
| `pei_slurm_univ_submit` | 通用引擎（`--mode` + 正交 `--chunks`），所有批量提交的唯一入口 |
| `pei_slurm_univ_vasp_544.sh` | VASP 单作业父脚本（`-n 128`），也作 each-subdir 的子作业 |
| `pei_slurm_univ_vasp_544_sequential_single_allocation.sh` | single-alloc 父脚本（`-n 128`） |
| `pei_slurm_univ_vasp_544_sequential_each_subdir_sbatch.sh` | each-subdir 父脚本（`-n 1`） |
| `vasp_utils/vasp_universal/pei_vasp_univ_load_env` | sourced：加载 module、把 `vasp_std` 放进 PATH、逐文件检查输入 |
| `vasp_utils/vasp_universal/pei_vasp_univ_sbatch` | 真正运行单个 VASP 计算的 runner（OUTCAR 判定 / CONTCAR→POSCAR / 清理 / `srun` / 退出码 0·10·1） |

之前的薄包装命令（`pei_vasp_univ_sbatch_parallel`、`pei_vasp_univ_sbatch_chunk`、`pei_vasp_univ_sbatch_sequential_*`、`pei_slurm_univ_sbatch_parallel`、`pei_slurm_run_vasp_544_*_chunk`）已删除，统一改为直接调用引擎 `pei_slurm_univ_submit`（用法见上文各节）。
