# `slurm_utils` 提交脚本

Slurm 提交脚本按用途分到三个子目录（均已加入 PATH，命令可在任意位置裸名调用）：

| 子目录 | 内容 |
|---|---|
| `slurm_universal/` | 软件无关的核心：Python CLI + launcher 重试 + 监控 + 命令速查 |
| `slurm_vasp/` | VASP 专属监控启动器 |
| `slurm_n2p2/` | n2p2 专属：训练 / scaling 作业脚本 |
| `docs/` | 架构与模式详解（见 `docs/submission-architecture-and-modes.md`） |

> 注意：运行下面的脚本前，必须先完成所有前处理，包括 `pei_vasp_univ_clean_up_full`、修改 `INCAR`、更新 `KPOINTS`。
> 只有 `pei_vasp_univ_sbatch` 会删除旧的输出文件（保留 `*.out`）：
>
> ```bash
> rm -f CHG CONTCAR* DOSCAR* DYNMAT EIGENVAL IBZKPT OPTIC OSZICAR* OUTCAR* PROCAR* PCDAT W* XDATCAR* PARCHG* vasprun.xml SUMMARY.* REPORT
> ```

## 架构

整套批量提交围绕**一个软件无关的 Python 引擎** `pei_slurm_univ_submit.py`：它遍历
`path_root/dir_root` 下的一级子目录，生成 Slurm 脚本，并按 `--mode` 选择提交方式。
VASP、LAMMPS 和 n2p2 的差异只体现在 module profile、launcher 和 `--cmd`。

```
pei_slurm_univ_submit.py --mode parallel|each-subdir|single-alloc [options]
    --path_root PATH             # 项目绝对根目录；默认当前目录
    --dir_root DIR               # 计算子目录的父目录；相对 path_root，默认 ./dir
    --lsubdir a b c              # 只处理这些一级子目录；默认自动扫描全部
    --chunks K                   # sequential 模式的并发调度流数量
    --chunk_parent_layout auto|shared|per-chunk
    --module_profile_type NAME   # module 环境配置
    --launcher_type srun|mpirun|none
    --cmd "CMD ..."              # 每个计算目录实际执行的命令
    --partition PART --nodes N --ncores C
    --if_sbatch [True|False]     # 默认 False：只生成脚本；裸写等价于 True
```

三种 `--mode`（每个目录的动作）：

| `--mode` | 每个目录做什么 | 父作业布局 |
|---|---|---|
| `parallel` | 生成并直接 `sbatch sub_slurm_univ.sh`，不等待 | 每目录一个独立计算作业 |
| `each-subdir` | chunk worker 逐个 `sbatch --wait sub_slurm_univ.sh` | 默认 1 个 `-n 1` shared 父作业管理 K 条 worker |
| `single-alloc` | chunk 父作业内直接执行 `cmd`，按退出码统计 | K 个持有完整计算资源的父作业 |

`--chunks K` 表示 **K 条并发调度流**，不再等同于 K 个父作业。默认
`--chunk_parent_layout auto` 的解析规则为：

- `each-subdir` → `shared`：只 `sbatch` 一个父作业，父作业内并发启动 K 个 Bash worker；
- `single-alloc` → `per-chunk`：仍 `sbatch` K 个父作业，每个 chunk 独占自己的计算资源；
- `parallel` → 忽略 chunks，因为本来就是每目录一个独立作业。

`each-subdir` 如需恢复历史 K 父作业行为，可显式传
`--chunk_parent_layout per-chunk`。`single-alloc` 不允许使用 `shared`。

## 常见用法

内置 preset 已提供常见环境和资源组合：

```bash
pei_slurm_univ_submit.py --list-presets
pei_slurm_univ_submit.py --show-preset zcm6-vasp-0

# each-subdir：5 条流、1 个 shared 父作业；不加 --if_sbatch 时只生成
pei_slurm_univ_submit.py --preset zcm6-vasp-0 --chunks 5 --if_sbatch

# each-subdir 历史布局：5 条流、5 个父作业
pei_slurm_univ_submit.py --preset zcm6-vasp-0 --chunks 5 \
    --chunk_parent_layout per-chunk --if_sbatch

# single-alloc：5 条流、5 个计算资源父作业
pei_slurm_univ_submit.py --preset zcm6-vasp-0 --mode single-alloc \
    --chunks 5 --if_sbatch
```

生成文件布局：

```text
<path_root>/
├─ <dir_root>/<case>/sub_slurm_univ.sh
└─ slurm/
   ├─ sub_slurm_each_subdir_chunk001.sh
   ├─ sub_slurm_each_subdir_chunk002.sh
   ├─ ...
   └─ sub_slurm_each_subdir_parent.sh      # shared 且 chunks > 1 时生成
```

chunk worker 文件名保持不变，可在失败恢复时单独 `sbatch`。shared parent 只引用本次
生成的明确 worker 列表，不会扫描并误执行旧运行残留的 chunk 文件。每个 worker 的详细
输出写到 `slurm/slurm-<parent-jobid>-sub_slurm_each_subdir_chunkNNN.out`，父作业 stdout
只记录 worker 启动信息和最终汇总。

## 软件无关的入口

引擎不绑定 VASP。选择 module profile、launcher 和 `cmd` 后，同一个入口可运行
n2p2、LAMMPS 或普通 shell 命令：

```bash
pei_slurm_univ_submit.py --preset zcm6-n2p2-train-0 --if_sbatch
pei_slurm_univ_submit.py --preset zcm6-lammps-0 --cmd "lmp -in lmp.in" --if_sbatch
```

## 文件一览（按子目录）

### `slurm_universal/`（软件无关核心）

| 文件 | 角色 |
|---|---|
| `pei_slurm_univ_submit.py` | Python CLI：preset / argparse / 参数回显，业务逻辑下沉到 `mymetal.slurm.submit` |
| `pei_slurm_univ_launch_retry` | 只对 Slurm/MPI 启动失败特征进行有限重试 |
| `pei_slurm_univ_monitor_error` | 扫描运行中作业 stdout，遇致命错误关键字时报警/处理 |
| `pei_slurm_univ_useful_command.sh` | 常用命令速查表（Slurm / VASP / …） |

### `slurm_vasp/`（VASP 专属）

| 文件 | 角色 |
|---|---|
| `pei_slurm_univ_vasp_monitor` | 便捷启动器：`sbatch` 监控作业 |

### `slurm_n2p2/`（n2p2 专属）

| 文件 | 角色 |
|---|---|
| `pei_slurm_univ_n2p2_train` | n2p2 训练辅助脚本 |
| `pei_slurm_univ_n2p2_scaling` | n2p2 scaling 辅助脚本 |

### 依赖（位于 `vasp_utils/vasp_universal/`，不在本目录）

| 文件 | 角色 |
|---|---|
| `pei_vasp_univ_sbatch` | 真正运行单个 VASP 计算的 runner（OUTCAR 判定 / CONTCAR→POSCAR / 清理 / `srun` / 退出码 0·10·1） |
| `pei_vasp_univ_monitor_error` | 扫描运行中作业 stdout，遇报错关键字即 `scancel` |

> CLI、runner 与 launcher 重试包装器都按裸名经 PATH 互相调用，因此对应工具目录必须在 PATH 上。
> 历史薄包装命令已经删除，统一使用 `pei_slurm_univ_submit.py`。
