slurm_utils 提交脚本
=====================

``slurm_utils`` 是软件无关的 Slurm 作业生成与提交工具集。核心引擎
``pei_slurm_univ_submit.py`` 递归发现 ``y_dir`` 计算目录，生成 Slurm 脚本，
并按 ``-mode`` 选择提交拓扑。VASP、LAMMPS 和 n2p2 的差异只体现在 module
profile、launcher 和 ``-cmd``。

.. note::

   这些脚本不注册 ``console_scripts``。请将 ``slurm_utils/slurm_universal/``
   等目录加入 ``PATH``，或使用完整路径调用。

.. contents:: 本页内容
   :local:
   :depth: 2

架构
----

所有批量提交都从 ``pei_slurm_univ_submit.py`` 进入。CLI 只负责 preset 和
argparse；目录切分、脚本生成与提交拓扑由 ``mymetal.slurm.submit`` 实现。

.. code-block:: text

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

三种 ``-mode`` （每个目录的动作）：

.. list-table::
   :header-rows: 1
   :widths: 20 45 35

   * - ``-mode``
     - 每个 case
     - 调度方式
   * - ``parallel``
     - 生成并直接提交独立计算作业
     - 每个 case 一个作业，不等待
   * - ``each_subdir``
     - 父/worker 依次 ``sbatch --wait`` 子作业
     - 默认一个 shared parent 管理 K 个 chunk
   * - ``single_alloc``
     - 在父 allocation 中直接执行 ``-cmd``
     - K 个 chunk 各持有计算资源

``-chunks K`` 表示 K 条并发调度流。默认 ``-chunk_parent_layout auto`` 的
解析规则：

- ``each_subdir`` → ``shared``：只 ``sbatch`` 一个父作业，父作业内并发启动
  K 个 Bash worker
- ``single_alloc`` → ``per_chunk``：仍 ``sbatch`` K 个父作业，每个 chunk
  独占自己的计算资源
- ``parallel`` → 忽略 chunks

slurm_universal
---------------

pei_slurm_univ_submit.py
~~~~~~~~~~~~~~~~~~~~~~~~

软件无关的作业生成/提交入口。

.. code-block:: bash

   pei_slurm_univ_submit.py -list_presets
   pei_slurm_univ_submit.py -show_preset zcm6_vasp_0

   # each_subdir：5 条流、1 个 shared 父作业；不加 -if_sbatch 时只生成
   pei_slurm_univ_submit.py -preset zcm6_vasp_0 -chunks 5 \
       -child_wall_time 2-00:00:00 -parent_wall_time 7-00:00:00 -if_sbatch

   # single_alloc：5 条流、5 个计算资源父作业
   pei_slurm_univ_submit.py -preset zcm6_vasp_0 -mode single_alloc \
       -chunks 5 -if_sbatch

参数说明：

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - 参数
     - 说明
   * - ``-path_root PATH``
     - 项目绝对根目录；默认当前目录
   * - ``-dir_root DIR``
     - ``y_dir`` 递归搜索根；相对 ``path_root``，默认当前目录
   * - ``-lsubdir a b c``
     - 按 basename 过滤全部 ``y_dir``；默认处理发现的全部
   * - ``-chunks K``
     - 并发调度流数量
   * - ``-chunk_parent_layout``
     - ``auto`` / ``shared`` / ``per_chunk``
   * - ``-module_profile_type``
     - module 环境配置（preset 或 ``none``）
   * - ``-launcher_type``
     - ``srun`` / ``mpirun`` / ``none``
   * - ``-cmd "CMD ..."``
     - 每个计算目录实际执行的命令
   * - ``-partition`` / ``-nodes`` / ``-ncores``
     - Slurm 资源参数
   * - ``-child_wall_time``
     - 计算子作业时限（仅 parallel / each_subdir）
   * - ``-parent_wall_time``
     - 父作业时限（仅 each_subdir / single_alloc）
   * - ``-if_sbatch``
     - 默认 ``False``：只生成脚本；裸写等价于 ``True``

目录发现以 ``dir_root`` 为递归搜索根：包含根本身在内，查找任意深度、名称
严格等于 ``y_dir`` 的目录，再把每个 ``y_dir`` 的一级子目录汇总为计算目录。

生成文件布局：

.. code-block:: text

   <path_root>/
   ├─ <dir_root>/**/y_dir/<case>/sub_slurm_univ.sh
   └─ slurm/
      ├─ sub_slurm_each_subdir_chunk001.sh
      ├─ sub_slurm_each_subdir_chunk002.sh
      ├─ ...
      └─ sub_slurm_each_subdir_parent.sh      # shared 且 chunks > 1 时生成

preset 注册表：

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - preset
     - module profile
     - 用途
   * - ``none``
     - 无 module
     - dry-run / 本地测试
   * - ``zcm6_vasp_0``
     - mpi/intel/17.0.7-thc + VASP 5.4.4
     - VASP 计算
   * - ``zcm6_lammps_0``
     - gcc/12.2 + mpi + fftw + LAMMPS
     - LAMMPS 计算
   * - ``zcm6_lammps_1``
     - 同上但 LAMMPS-NC 版本
     - LAMMPS NC 计算
   * - ``zcm6_n2p2_0``
     - eigen + gsl + openmpi + n2p2
     - n2p2 训练/scaling

pei_slurm_univ_launch_retry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

srun/mpirun 启动失败重试包装器。只重试 "launcher 没起来" 的暂态错误
（调度器抖动、作业步创建被禁用），**不重试** 程序本身的退出码 （VASP 不收敛、
段错误、OOM 被杀等重试多少次都是白烧机时）。

.. code-block:: bash

   pei_slurm_univ_launch_retry srun -n 128 vasp_std
   pei_slurm_univ_launch_retry mpirun -np 24 nnp-train

   # 经 MY_LAUNCHER 传给下游
   export MY_LAUNCHER="pei_slurm_univ_launch_retry srun -n 128"

可调参数（环境变量）：

- ``PEI_LAUNCH_MAX_RETRY``：最大重试次数（默认 99）
- ``PEI_LAUNCH_SLEEP``：每次重试间隔秒数（默认 10）

pei_slurm_univ_sbatch_retry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sbatch 提交失败重试包装器。与 ``launch_retry`` 分工互补：

- ``launch_retry`` 包 ``srun``/``mpirun``，处理"作业内启动器起不来"
- ``sbatch_retry`` 包 ``sbatch``，处理"提交这一步没被 slurmctld 接住"

提交成功判据（命中其一即算成功，绝不重试）：

1. ``sbatch`` 退出码为 0
2. 输出含 ``"Submitted batch job"``

.. code-block:: bash

   pei_slurm_univ_sbatch_retry --wait sub_slurm_univ.sh
   pei_slurm_univ_sbatch_retry sub.544.sh
   pei_slurm_univ_sbatch_retry sub.*

.. note::

   提交失败重试是安全的——没有作业被创建，不会产生重复作业。这与
   ``launch_retry`` 不同（后者重试时握着整个分配空转）。

pei_slurm_univ_monitor_error
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

通用作业监控器提交入口。按 ``-<engine>`` 选出对应的监控脚本，生成一个
1 节点 1 核的看守作业并 ``sbatch`` 出去。

.. code-block:: bash

   pei_slurm_univ_monitor_error -vasp [-skip_ljobid JOBID,...] [-phase_check]

目前支持的引擎：

- ``-vasp`` → ``pei_vasp_univ_monitor_error``：VASP 错误关键字检测、
  迭代步数监控、可选 CNA 相变检查

可调参数（环境变量）：

- ``PEI_MONITOR_PARTITION``：看守作业分区（默认 ``amd_512``）
- ``PEI_MONITOR_JOBNAME``：作业名（默认 ``pei_monitor``）
- ``PEI_MONITOR_DIR``：脚本和输出落脚点（默认 ``./slurm_monitor``）

pei_slurm_univ_useful_command.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

命令速查表 + 网络探测工具。

.. code-block:: bash

   pei_slurm_univ_useful_command.sh                  # 速查表
   pei_slurm_univ_useful_command.sh -claude          # 探测 Anthropic/Claude
   pei_slurm_univ_useful_command.sh -openai -proxy   # 探测 OpenAI + 显示代理
   pei_slurm_univ_useful_command.sh -net -ip         # 探测全部 + 公网 IP
   pei_slurm_univ_useful_command.sh -monitor          # 汇总运行中 Slurm 作业
   pei_slurm_univ_useful_command.sh -all -timeout 3

slurm_vasp
----------

pei_slurm_univ_vasp_monitor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

VASP 监控作业的便捷启动器。内部执行
``sbatch pei_slurm_univ_vasp_monitor.sh``。

slurm_n2p2
----------

pei_slurm_univ_n2p2_train
~~~~~~~~~~~~~~~~~~~~~~~~~~

n2p2 训练作业脚本。固定配置：``amd_512`` 分区，1 节点 24 核，加载 eigen/gsl/
openmpi module，执行 ``srun nnp-train``。

pei_slurm_univ_n2p2_scaling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

n2p2 scaling 作业脚本。固定配置：``amd_512`` 分区，1 节点 1 核，执行
``srun nnp-scaling 10000``。

.. note::

   这两个 n2p2 脚本是固定配置的 ``#SBATCH`` 脚本，不是通用入口。如需自定义
   分区或核数，直接编辑脚本或在 ``pei_slurm_univ_submit.py`` 中用
   ``-preset zcm6_n2p2_0`` 配合 ``-cmd "nnp-train"`` 使用。

最小 dry-run
------------

假设 ``/abs/project/y_dir/001`` 与 ``002`` 已存在：

.. code-block:: bash

   python slurm_utils/slurm_universal/pei_slurm_univ_submit.py \
       -path_root /abs/project -dir_root . \
       -mode each_subdir -chunks 2 \
       -module_profile_type none -launcher_type none \
       -cmd "echo dry-run" \
       -partition debug -nodes 1 -ncores 1

命令不含 ``-if_sbatch``，应生成：

.. code-block:: text

   /abs/project/
   |-- y_dir/
   |   |-- 001/sub_slurm_univ.sh
   |   `-- 002/sub_slurm_univ.sh
   `-- slurm/
       |-- sub_slurm_each_subdir_chunk001.sh
       |-- sub_slurm_each_subdir_chunk002.sh
       `-- sub_slurm_each_subdir_parent.sh

检查后再提交：

.. code-block:: bash

   bash -n /abs/project/y_dir/001/sub_slurm_univ.sh
   bash -n /abs/project/slurm/sub_slurm_each_subdir_parent.sh
   grep -R "^#SBATCH" /abs/project/slurm /abs/project/y_dir

示例脚本
--------

.. literalinclude:: ../../examples/slurm_utils_overview.py
   :language: python
   :caption: slurm_utils 子包功能总览示例

.. seealso::

   - :doc:`slurm` — Slurm 提交与 dry-run （简版）
   - :doc:`vasp_universal` — VASP 单目录工具集
   - :doc:`vasp_workflow_bulk` — VASP bulk workflow
   - :doc:`workflows` — Workflow 总览
   - :doc:`../reference/scripts` — 完整脚本参考
   - ``slurm_utils/docs/submission-architecture-and-modes.md`` — 架构详解
