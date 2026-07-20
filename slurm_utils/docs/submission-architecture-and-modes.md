# 提交脚本：架构、模式与流程详解

本文档讲清楚 `slurm_utils/`（`slurm_universal/` 通用核心、`slurm_vasp/` VASP 脚本）与 `vasp_utils/vasp_universal/` 这套批量提交脚本：
**每个文件干什么、三种模式各自的作用和流程、以及端到端怎么跑起来**。

> 速查手册见 [`../README.md`](../README.md)；本文是其背后的原理与流程说明。

---

## 一、整体心智模型

所有批量提交都从 `pei_slurm_univ_submit.py` 进入。CLI 只负责 preset 和 argparse；目录
切分、脚本生成与提交拓扑由 `mymetal.slurm.submit` 实现。父脚本和每目录 base 脚本均为
运行时自动生成文件，不再依赖历史固定父 wrapper。

目录发现以 `dir_root` 为递归搜索根：包含根本身在内，查找任意深度、名称严格等于
`y_dir` 的目录，再把每个 `y_dir` 的一级子目录汇总为计算目录。因此从包含多个
`mode/.../y_dir/<case>` 的顶层提交时，所有 mode 会进入同一次排序和 chunk 切分。
显式传入 `-dir_root` 不会退回一级扫描；preset 默认使用当前目录 `.`。

```
pei_slurm_univ_submit.py
  └─ mymetal.slurm.submit.pei_slurm_univ_submit
       ├─ parallel
       │    └─ 每目录 sbatch sub_slurm_univ.sh
       ├─ each_subdir -chunks K
       │    └─ 1 个 shared parent (-n1)
       │         ├─ chunk001 worker ── 逐个 sbatch --wait 子作业
       │         ├─ ...
       │         └─ chunk00K worker ── 逐个 sbatch --wait 子作业
       └─ single_alloc -chunks K
            └─ K 个计算资源父作业，各自在分配内顺序执行 cmd
```

**源码层：**

| 文件 | 角色 |
|---|---|
| `slurm_utils/slurm_universal/pei_slurm_univ_submit.py` | 薄 CLI、preset、module profile 注册表 |
| `mymetal/slurm/submit.py` | 软件无关的目录切分、脚本生成、父作业布局与提交逻辑 |
| `pei_vasp_univ_sbatch` | **真正跑一个 VASP 计算的 runner**：判收敛、断点续算、清理、`srun`、给退出码 |

**运行时生成层：**

| 文件 | 资源 | 角色 |
|---|---|---|
| `<case>/sub_slurm_univ.sh` | 用户指定 | parallel / each_subdir 的单目录计算子作业 |
| `slurm/sub_slurm_each_subdir_chunkNNN.sh` | `-n 1` | 可独立 sbatch 的 chunk worker，保留历史恢复接口 |
| `slurm/sub_slurm_each_subdir_parent.sh` | `-n 1` | 并发启动并等待本次 K 个 worker 的 shared parent |
| `slurm/sub_slurm_single_alloc_chunkNNN.sh` | 用户指定 | single_alloc 的每 chunk 计算资源父作业 |

---

## 二、最底层：runner `pei_vasp_univ_sbatch` 的流程

这是所有 VASP 计算的"原子操作"，**在一个已有的 Slurm 分配内**对单个目录执行。
理解它，上面所有模式才讲得通。

```
pei_vasp_univ_sbatch <dir> [exe] [完成标记] [EDIFF过大标记]
  │
  1. cd 进 dir            ── 进不去 → exit 1
  2. OUTCAR 含"reached required accuracy"？
        是 → 打印"已完成，跳过" → exit 10   ★关键安全闸：已完成就此返回，不删任何东西
  3. 断点续算判断：
        CONTCAR 非空 → cp CONTCAR POSCAR（从上次结构续算）
                       且若 OUTCAR 含"EDIFF 过大" → 自动改 EDIFF=1e-10、ALGO=Normal
        否则 POSCAR 在 → 从 POSCAR 全新开始
        否则           → exit 1（没有可用输入）
  4. 清理旧输出：rm -f CHG CONTCAR* ... OUTCAR* ... vasprun.xml ...   ★会删 OUTCAR，所以第2步必须在它之前
  5. srun vasp_std
  6. 跑完 OUTCAR 又含完成标记？
        是 → exit 0（本次跑完成）
        否 → exit 1（仍未收敛）
```

**退出码约定（整套系统的通用语言）：**

| 退出码 | 含义 |
|---|---|
| `0` | 本次运行后完成 |
| `10` | 早已完成，被跳过（不删除、不 srun） |
| `1` | 失败或缺输入 |

> 第 2 步是安全闸：已收敛的目录 runner 一进去就 `exit 10`，根本走不到第 4 步的删除。
> 因此对已完成的目录重复运行是安全的。

---

## 三、引擎的三种 `-mode`：作用与流程

引擎对每个目录做的动作由 `-mode` 决定。`nodes/ncores` 始终描述真正的计算资源；
`each_subdir` 的 worker/shared parent 只负责编排，会另外固定成 `-N 1 -n 1`。

| `-mode` | 每个目录的动作 | 父作业占核 | `chunks=K` 时的布局 |
|---|---|---|---|
| `parallel` | `sbatch sub_slurm_univ.sh`，不等 | 子作业各自申请 | chunks 被忽略 |
| `each_subdir` | worker 内 `sbatch --wait sub_slurm_univ.sh` | shared parent / worker 均 `-n 1` | 默认 1 父 + K worker |
| `single_alloc` | chunk 父作业内直接执行 cmd | 每父作业持有用户指定资源 | K 个资源父作业 |

引擎只做结构性检查和编排，不提前读取 OUTCAR 等业务文件。是否完成/跳过由 `cmd`
或下游 runner 自己判定，沿用退出码 `0/10/其它` 契约。

`parallel` / `each_subdir` 可用 `-child_wall_time 2-00:00:00` 给每个
`<case>/sub_slurm_univ.sh` 写入 `#SBATCH --time=2-00:00:00`。不传该选项时不生成
`--time` 行；它不限制 each_subdir 的 shared parent/chunk worker，`single_alloc` 也忽略它。

父脚本使用独立的 `-parent_wall_time 7-00:00:00`。它作用于 `each_subdir` 的 shared
parent / 可独立恢复的 chunk worker，以及 `single_alloc` 的 chunk 父脚本；`parallel`
没有父作业，会忽略该参数。shared parent 用普通 Bash 启动 worker 时，worker 内的
`#SBATCH` 指令只是注释，不会产生二次资源申请。

### 3.1 `parallel` —— 全部铺开，互不等待

```
你在登录节点跑引擎
  └─ 对每个 dir：
        生成 sub_slurm_univ.sh
        sbatch sub_slurm_univ.sh   （丢出去就走，不等）
  └─ 汇总提交成功/失败
       ↓（调度器各自排队、并行起跑）
   每个子作业按 module_profile + launcher + cmd 运行
```

- **作用**：吞吐最大化。N 个目录 = N 个独立作业同时排队。
- **完成判定下推**：引擎不 grep 业务文件；VASP 可让 `pei_vasp_univ_sbatch` 自行判断。
- **`failed` 含义**：这里只表示 `sbatch` 提交命令本身失败，不代表计算内容失败。

### 3.2 `each_subdir` —— 顺序，每个目录一个独立子作业

```
sbatch shared parent (-n 1)
  └─ 并发启动 K 个 chunk worker（普通后台 Bash 进程）
       └─ 每个 worker 对本组 dir 顺序执行：
            sbatch --wait sub_slurm_univ.sh
            根据子作业退出码统计完成/未收敛
  └─ wait 全部 worker，聚合成功/失败后 exit 0/1
```

- **作用**：每条 worker 内串行、worker 之间并发；每个计算仍是独立 Slurm 子作业。
- **三级日志**：shared parent stdout 记录 worker 启动/聚合；每个 worker 写独立日志；
  每个计算子作业日志仍在对应计算目录。
- **特点**：一个目录算崩不影响别的（子作业相互独立），父进程继续往下走。
- **历史恢复**：chunk worker 文件继续保留，可用 `sbatch chunkNNN.sh` 单独重跑。

### 3.3 `single_alloc` —— 顺序，复用同一个分配

```
sbatch sub_slurm_single_alloc_chunkNNN.sh（用户指定资源）
  └─ 对每个 dir：
        在 dir 内直接执行 launcher + cmd
        看退出码：0 完成 / 10 跳过 / 其它 失败（计数后继续下一个）
  └─ 汇总 completed/skipped/failed
```

- **作用**：每个 chunk 只排一次队；一次拿到该 chunk 的完整资源，在同一个分配里逐个
  `srun/mpirun` 跑完本组目录。
- **跳过靠 runner**：真正运行时，已收敛的目录 runner 返回 10，引擎记 skip。
  不传 `-if_sbatch` 时只生成脚本，不会调用 runner 或 launcher。
- **引擎层遇错不中止**：某目录未收敛只计入 `failed` 并继续下一个。
- **可选的外部中止**：若另外起一个 `pei_vasp_univ_monitor_error`（独立 `-n 1` 监控作业）
  在旁边轮询，它会扫描各作业 stdout，一旦发现致命错误关键词
  （`hermitian` / `sloshing` / `ERROR in subspace rotation PSSYEVX` / `highest band`）
  就 `scancel` 该作业——对 single_alloc 而言那个 job 就是一个 chunk 父作业，于是该组
  剩余目录无法继续，其它 chunk 父作业不受影响。

---

## 四、`-chunks` —— 正交叠加，把核预算铺成 K 条并发流

chunk 不是第四种模式。`-chunks K` 只负责把排序后的目录用 base/rem 算法均分成
K 条调度流；父作业数量由 `-chunk_parent_layout` 决定：

```
-chunk_parent_layout auto
  ├─ each_subdir  → shared：1 个父作业 + K 个后台 worker
  ├─ single_alloc → per_chunk：K 个计算资源父作业
  └─ parallel     → chunks 被忽略
```

- **each_subdir/shared**：K 个 worker 同时各自执行 `sbatch --wait`，所以仍可有最多 K 个
  计算子作业同时运行；调度器中只多出一个 `-n 1` 编排父作业。
- **each_subdir/per_chunk**：显式兼容开关；K 个 chunk 文件各自 `sbatch`，恢复历史 K 父作业。
- **single_alloc/per_chunk**：每条流需要真实计算资源，必须是 K 个独立父分配；不允许 shared。
- **shared 汇总**：父脚本等待全部 worker，任一 worker 失败则父作业最终退出 1。worker
  独立日志避免并发输出交错。
- **残留文件安全**：shared parent 嵌入本次生成的 worker 绝对路径列表，不使用 glob，旧的
  `chunk006.sh` 等文件即使保留也不会被误执行。

---

## 五、端到端调用链

| # | 场景 | 你敲的命令 | 调用链 |
|---|---|---|---|
| 1 | 并行全铺开 | `pei_slurm_univ_submit.py -preset zcm6_vasp_0 -mode parallel -if_sbatch` | 引擎 → N×`sbatch sub_slurm_univ.sh` |
| 2 | each_subdir shared | `pei_slurm_univ_submit.py -preset zcm6_vasp_0 -chunks 4 -if_sbatch` | 1 父 → 4 worker → 逐目录子作业 |
| 3 | each_subdir 历史布局 | 上一命令加 `-chunk_parent_layout per_chunk` | 4 父 → 各自逐目录子作业 |
| 4 | single_alloc | 上一命令改 `-mode single_alloc` | 4 个资源父作业 → 各自在分配内逐目录执行 |

**通用旋钮：**

- `-dir_root PATH`：从该位置递归发现所有 `y_dir/<case>`；默认当前提交位置。
- `-lsubdir a b c`：按 basename 过滤递归结果；多个 `y_dir` 中的同名目录会全部选中。
- `-child_wall_time D-HH:MM:SS`：可选的计算子作业 wall time；不设置则不额外限制。
- `-parent_wall_time D-HH:MM:SS`：可选的父作业 wall time；不设置则不额外限制。
- 不传 `-if_sbatch`：只生成/覆盖脚本，不提交作业。
- `-if_sbatch False`：与不传相同；裸写 `-if_sbatch` 等价于 True。
- 业务完成检查放在 `-cmd` 或 runner 内，不在通用编排层提前读取结果文件。

---

## 六、软件无关性

引擎不绑定 VASP。切换 preset 或显式设置 module profile、launcher、cmd 即可服务
n2p2 / LAMMPS 等工作流：

```bash
pei_slurm_univ_submit.py -preset zcm6_n2p2_train_0 -if_sbatch
pei_slurm_univ_submit.py -preset zcm6_lammps_0 \
    -cmd "lmp -in lmp.in" -if_sbatch
```

VASP 相关的 `pei_vasp_univ_sbatch` 只是建在通用引擎之上的业务 runner；父脚本仍由
同一个 `mymetal.slurm.submit` 动态生成。
