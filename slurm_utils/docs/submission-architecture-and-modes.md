# 提交脚本：架构、模式与流程详解

本文档讲清楚 `slurm_utils/slurm_universal/` 与 `vasp_utils/vasp_universal/` 这套批量提交脚本：
**每个文件干什么、三种模式各自的作用和流程、以及端到端怎么跑起来**。

> 速查手册见 [`../slurm_universal/README.md`](../slurm_universal/README.md)；本文是其背后的原理与流程说明。

---

## 一、整体心智模型

整套东西就 **4 个"活"文件 + 3 个父作业文件**，分四层。所有批量提交最终都汇聚到**一个引擎**，
引擎再把"每个目录具体做什么"交给底层 runner 或子作业。

```
你敲的命令
   │
   ├─ sbatch 父.sh ───────────┐         （场景 1/3/4：需要 #SBATCH 资源声明）
   └─ pei_slurm_univ_submit ──┤         （场景 2/5：登录节点直接调引擎）
                              ▼
                ┌─────────────────────────────┐
                │  pei_slurm_univ_submit       │  引擎：遍历目录 + 按 --mode 派活 + 汇总
                │  （source lib）              │
                └──────────────┬──────────────┘
            parallel │ each-subdir │ single-alloc │ --chunks
                     ▼           ▼              ▼          ▼
              sbatch 子作业  sbatch --wait   在分配内直接   sbatch K 个
              （不等）       子作业（等）     调 runner      父作业
                                 │              │
                                 ▼              ▼
                       pei_slurm_univ_vasp_544.sh   pei_vasp_univ_sbatch
                       （-n128 子作业，内部又调 →）  （真正跑 VASP 的 runner）
```

**地基层（删任何一个全崩）：**

| 文件 | 角色 |
|---|---|
| `pei_slurm_univ_lib.sh` | 共享库，被 source。提供 `ps_ok/ps_warn/ps_fail`、preflight 回显、**目录解析** `ps_resolve_job_dirs`、脚本路径解析、失败汇总 |
| `pei_slurm_univ_submit` | **通用引擎**。软件无关，只懂"遍历目录 → 派活 → 汇总" |
| `pei_vasp_univ_load_env` | VASP 专属，被父 `.sh` source：加载 module、把 `vasp_std` 放进 PATH、逐文件 ✅ 检查输入 |
| `pei_vasp_univ_sbatch` | **真正跑一个 VASP 计算的 runner**：判收敛、断点续算、清理、`srun`、给退出码 |

**父作业层（必须是文件，因为 `#SBATCH` 指令只能写在文件里）：**

| 文件 | 资源 | 角色 |
|---|---|---|
| `pei_slurm_univ_vasp_544.sh` | `-n 128` | 单作业父脚本；也作 each-subdir 模式的子作业 |
| `pei_slurm_univ_vasp_544_sequential_single_allocation.sh` | `-n 128` | single-alloc 父脚本 |
| `pei_slurm_univ_vasp_544_sequential_each_subdir_sbatch.sh` | `-n 1` | each-subdir 父脚本（只调度，不跑 VASP） |

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

## 三、引擎的三种 `--mode`：作用与流程

引擎对每个目录做的"动作"由 `--mode` 决定。三种模式最大的区别有两个：
① 每个目录干什么；② **引擎本身在哪运行**（登录节点 vs. 分配内）。

| `--mode` | 每个目录的动作 | 引擎本身在哪跑 | 父作业占核 | 并发度 |
|---|---|---|---|---|
| `parallel` | `sbatch 子脚本`，不等 | **登录节点** | 子作业各自申请 128 | 全部目录同时排队 |
| `each-subdir` | `sbatch --wait 子脚本`，完后复检 | **-n 1 父分配内** | 1 核父 + 逐个 128 子 | 一次一个子作业 |
| `single-alloc` | 分配内直接调 runner | **-n 128 父分配内** | 父独占 128 | 一次一个，复用同一分配 |

引擎对每个目录共通的前置动作（来自 lib）：
解析目录列表（`--root-dir` 全部 / `--lsubdir a,b,c` 子集）→ 逐目录进入 →
可选 `--skip-if-file-contains` 预过滤 → 末尾打印 `submitted/converged/skipped/failed` 汇总 + 失败目录清单。

### 3.1 `parallel` —— 全部铺开，互不等待

```
你在登录节点跑引擎
  └─ 对每个 dir：
        OUTCAR 已含完成标记？→ skip（连提交都省了）
        否则 → sbatch pei_slurm_univ_vasp_544.sh   （丢出去就走，不等）
  └─ 汇总 submitted/skipped/failed
       ↓（调度器各自排队、并行起跑）
   每个子作业 = pei_slurm_univ_vasp_544.sh (-n128)
        内部：preflight → load_env → 检查输入 → pei_vasp_univ_sbatch . vasp_std
```

- **作用**：吞吐最大化。N 个目录 = N 个独立作业同时排队。
- **双重跳过**：引擎先 grep OUTCAR 预过滤；子作业里的 runner 还会再判一次（exit 10）。
- **`failed` 含义**：这里只表示 `sbatch` 提交命令本身失败，不代表 VASP 算错（算得怎样要看各子作业日志）。

### 3.2 `each-subdir` —— 顺序，每个目录一个独立子作业

```
sbatch 父.sh (-n 1)  →  父分配内 exec 引擎
  └─ 对每个 dir：
        已含完成标记？→ skip
        否则 → sbatch --wait pei_slurm_univ_vasp_544.sh   （阻塞，等子作业结束）
               子作业结束后 → 复检 OUTCAR 标记 → converged / 未收敛
  └─ 汇总 submitted/converged/skipped/failed
```

- **作用**：串行但隔离。父作业只调度（所以只要 1 核、不加载 VASP module），
  每个计算是**独立子作业**（各自 128 核、各自日志）。
- **两级日志**：父日志（进入哪个目录、提交、复检结果）在提交目录；
  每个子作业日志在对应计算目录。
- **特点**：一个目录算崩不影响别的（子作业相互独立），父进程继续往下走。

### 3.3 `single-alloc` —— 顺序，复用同一个分配

```
sbatch 父.sh (-n 128)  →  父分配内 exec 引擎
  └─ 对每个 dir：
        在 dir 内直接 eval "pei_vasp_univ_sbatch . vasp_std ..."   （就在这 128 核里 srun）
        看退出码：0 完成 / 10 跳过 / 其它 失败（计数后继续下一个）
  └─ 汇总 completed/skipped/failed
```

- **作用**：**只排一次队**。一次拿到 128 核，在同一个分配里逐个 `srun` 跑完所有目录，
  省去反复排队的等待。
- **跳过靠 runner**：真正运行时，已收敛的目录 runner 返回 10，引擎记 skip。
  传入 `--dry-run` 时，引擎只打印每个目录将执行的 runner 命令，不会调用 runner 或 `srun`。
- **引擎层遇错不中止**：`set +e; eval; set -e`，某目录未收敛只计入 `failed` 并继续下一个。
- **可选的外部中止**：若另外起一个 `pei_vasp_univ_monitor_error`（独立 `-n 1` 监控作业）
  在旁边轮询，它会扫描各作业 stdout，一旦发现致命错误关键词
  （`hermitian` / `sloshing` / `ERROR in subspace rotation PSSYEVX` / `highest band`）
  就 `scancel` 该作业——对 single-alloc 而言那个 job 就是整个父作业，于是整批顺序流被中止、
  无法继续剩余目录。这就是 README 中 "monitor 脚本会终止整个父作业" 的来源。

---

## 四、`--chunks` —— 正交叠加，把核预算铺成 K 条并发流

chunk 不是第四种模式，而是**叠在 sequential 模式上的一个开关**：

```
你在登录节点跑：
  pei_slurm_univ_submit --chunks K --chunk-parent <某个sequential父.sh> --root-dir ./y_dir
  └─ 把一级子目录均分成 K 组（base/rem 算法，余数摊到前几组）
  └─ 对每组 → sbatch <父.sh> "$root_dir" "$group"
       ↓
   K 个父作业并发起跑，每个父作业内部就是 §3.2 或 §3.3 的顺序流
```

- **意义**：单个 sequential 父作业是串行的；切 K 块 = K 个父作业并发，
  等于把总核预算（K×128）铺成 **K 条并发的顺序流**，
  在"全并行会瞬间占满队列"和"全串行太慢"之间取平衡。
- **只对 sequential 有意义**：`parallel` 本来就最大并发，再 chunk 没用。
- `--chunk-parent` 用哪个 `.sh`，就决定每条流内部是 `each-subdir` 还是 `single-alloc`。

---

## 五、五个场景的端到端调用链

| # | 场景 | 你敲的命令 | 调用链 |
|---|---|---|---|
| 1 | 单目录 | `sbatch pei_slurm_univ_vasp_544.sh` | 父.sh → runner |
| 2 | 并行全铺开 | `pei_slurm_univ_submit --mode parallel --root-dir ./y_dir --submit-script pei_slurm_univ_vasp_544.sh --skip-if-file-contains OUTCAR "reached required accuracy"` | 引擎 → N×(sbatch 父.sh → runner) |
| 3 | 顺序·同一分配 | `sbatch pei_slurm_univ_vasp_544_sequential_single_allocation.sh [root] [lsubdir]` | 父.sh(-n128) → 引擎 single-alloc → 逐个 runner |
| 4 | 顺序·每目录子作业 | `sbatch pei_slurm_univ_vasp_544_sequential_each_subdir_sbatch.sh [root] [lsubdir]` | 父.sh(-n1) → 引擎 each-subdir → 逐个 sbatch --wait 父.sh → runner |
| 5 | 分块并发流 | `pei_slurm_univ_submit --chunks 4 --root-dir ./y_dir --chunk-parent <场景3或4的父.sh>` | 引擎 chunk → K×(sbatch 父.sh → 场景3/4 流程) |

**通用旋钮：**

- `--lsubdir "a,b,c"`：只跑这几个一级子目录（默认全部）。
- `--dry-run`：只打印动作、不提交/不运行。除 `single-alloc` 外都支持。
- `--skip-if-file-contains FILE MARKER`：当目录内 `FILE` 含 `MARKER` 时跳过该目录。

---

## 六、软件无关性

引擎 `pei_slurm_univ_submit` 不绑定 VASP——把 `--submit-script` 换成别的 Slurm 脚本即可用于
n2p2 / LAMMPS 等任何 `sbatch` 工作流：

```bash
pei_slurm_univ_submit --mode parallel --root-dir ./train_dirs \
    --submit-script pei_slurm_univ_n2p2_train.sh --dry-run
```

VASP 相关的部分（`pei_vasp_univ_load_env` 环境加载、`pei_vasp_univ_sbatch` 计算 runner、
三个 `pei_slurm_univ_vasp_544*.sh` 父脚本）是建在通用引擎之上的 VASP 适配层。
