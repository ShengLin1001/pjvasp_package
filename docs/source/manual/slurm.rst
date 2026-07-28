Slurm 提交与 dry-run
====================

统一入口
--------

``slurm_utils/slurm_universal/pei_slurm_univ_submit.py`` 是软件无关入口。
它递归发现 ``y_dir/<case>``，生成计算脚本和可选父脚本。默认
``-if_sbatch=False``，因此只生成、不提交。

三种模式
--------

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
     - 在父 allocation 中直接运行 ``-cmd``
     - K 个 chunk 各持有计算资源

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

检查后再提交
------------

.. code-block:: bash

   bash -n /abs/project/y_dir/001/sub_slurm_univ.sh
   bash -n /abs/project/slurm/sub_slurm_each_subdir_parent.sh
   grep -R "^#SBATCH" /abs/project/slurm /abs/project/y_dir

确认 partition、account、node/core、wall time、module block、launcher、实际
``-cmd`` 和日志路径。只有明确接受生成内容后，才在同一命令末尾添加
``-if_sbatch``。

状态与恢复
----------

* dry-run：只有文件，没有 job ID。
* submitted/pending：使用 ``squeue``/``scontrol`` 查看原因。
* running：查看 parent、worker 和 child stdout，区分编排错误与计算错误。
* failed/timeout：先确认 child 是否仍在队列，再单独恢复明确的 chunk，避免重复作业。
* completed：调度器成功仍需检查程序退出码和科学收敛。

``-chunks K`` 表示 K 条并发调度流；在 ``each_subdir`` 默认 shared parent
布局下，并不表示 K 个父作业。完整架构和 preset 以
``slurm_utils/docs/submission-architecture-and-modes.md`` 与 CLI ``-h`` 为准。
