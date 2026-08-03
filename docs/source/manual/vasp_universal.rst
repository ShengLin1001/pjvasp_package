vasp_utils/vasp_universal 单目录工具集
========================================

``vasp_utils/vasp_universal`` 是 VASP 计算的单目录管理工具集。这些脚本不是
Python package 的一部分，而是独立的 bash/python 命令行工具，用于在 HPC 上管理
VASP 计算目录的运行、续算、后处理、清理和监控。

.. note::

   这些脚本不注册 ``console_scripts``。请将 ``vasp_utils/vasp_universal/`` 加入
   ``PATH``，或使用完整路径调用。首次使用应先 ``-h`` 或阅读源码头部注释。

.. contents:: 本页内容
   :local:
   :depth: 2

核心 runner：pei_vasp_univ_sbatch
----------------------------------

``pei_vasp_univ_sbatch`` 是所有 VASP 计算的"原子操作"——它在**一个已有的 Slurm
分配内**对单个目录执行 VASP 计算。理解它，上面所有提交模式才讲得通。

功能
~~~~

1. **进入目录**：``cd <job_dir>``，进不去则 ``exit 1``
2. **完成检查**：OUTCAR 含完成标记 → ``exit 10`` （已完成，跳过，不删任何东西）
3. **断点续算**：CONTCAR 非空 → ``cp CONTCAR POSCAR``；若 OUTCAR 含 "EDIFF 过大"
   → 自动改 ``EDIFF`` = ``1e-10``、``ALGO`` = ``Normal``
4. **清理旧输出**：删除 CHG、CONTCAR、OUTCAR 等（第 2 步必须在它之前）
5. **运行 VASP**：``srun vasp_std``
6. **退出码判定**：OUTCAR 又含完成标记 → ``exit 0``；否则 ``exit 1``

退出码约定
~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - 退出码
     - 含义
   * - ``0``
     - 本次运行完成计算
   * - ``10``
     - 计算已完成，跳过（安全闸：不删任何东西）
   * - ``1``
     - 计算失败或所需文件缺失

完成标记自动检测
~~~~~~~~~~~~~~~~

脚本根据 INCAR 的 ``NSW`` 和 ``IBRION`` 自动判断是**离子弛豫**还是**静态计算**：

- 离子弛豫（``NSW≠0`` 且 ``IBRION≠-1``）：标记 ``reached required accuracy``
- 静态计算（``NSW=0`` 或 ``IBRION=-1``）：标记 ``aborting loop because EDIFF is reached``

可通过第三个参数手动覆盖。

用法
~~~~

.. code-block:: bash

   pei_vasp_univ_sbatch <job_dir> [vasp_executable] [completion_marker] [smaller_ediff_marker]

   # 默认：vasp_std + 自动检测标记
   pei_vasp_univ_sbatch /abs/project/y_dir/001

   # 指定可执行文件
   pei_vasp_univ_sbatch /abs/project/y_dir/001 vasp_gam

   # 手动指定完成标记
   pei_vasp_univ_sbatch /abs/project/y_dir/001 vasp_std "reached required accuracy"

流程图
~~~~~~

.. figure:: /_static/images/generated/vasp_universal_overview.png
   :alt: vasp_universal 子包核心功能总览

   ``vasp_utils/vasp_universal`` 子包功能总览：runner 流程、目录扫描、清理对比、
   INCAR 操作和退出码约定。

.. seealso::

   提交入口 :doc:`slurm` 和 :doc:`workflows` 说明了 ``pei_vasp_univ_sbatch``
   如何被 ``pei_slurm_univ_submit.py`` 的三种 ``-mode`` 编排。

状态汇总：pei_vasp_univ_post
-----------------------------

``pei_vasp_univ_post`` 递归发现 ``y_dir``，汇总每个 ``y_dir`` 一级子目录的计算
状态，写出多个汇总文件。

产出文件
~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - 文件
     - 内容
   * - ``y_post_data.txt``
     - 每个目录的能量、体积等数据
   * - ``y_post_data_2.txt``
     - 扩展数据
   * - ``y_post_time.txt``
     - 计算时间统计
   * - ``y_post_param.txt``
     - INCAR 参数汇总
   * - ``y_post_param_2.txt``
     - 扩展参数
   * - ``y_post_diff.txt``
     - 差异信息
   * - ``y_post_warning.txt``
     - 警告信息

用法
~~~~

.. code-block:: bash

   # 基本用法：在包含 y_dir 的目录运行
   pei_vasp_univ_post

   # verbose 模式：额外调用绘图脚本
   pei_vasp_univ_post -v

清理脚本
--------

``vasp_universal`` 提供四个清理脚本，各有不同的清理范围：

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - 脚本
     - 清理范围
   * - ``pei_vasp_univ_clean_up_full``
     - 全部输出 + slurm 日志 + ``y_`` 文件 + ML 临时文件（最彻底）
   * - ``pei_vasp_univ_clean_up_small``
     - VASP 输出 + slurm 日志（保留 ``y_`` 文件和 ML 临时文件）
   * - ``pei_vasp_univ_clean_old_slurm``
     - 仅保留每个 y_dir 子目录中最新的 slurm 文件，删除其余
   * - ``pei_vasp_univ_clean_outcar``
     - 清理 OUTCAR 末尾的 null 字符（不删 OUTCAR）

.. warning::

   ``pei_vasp_univ_clean_up_full`` 和 ``clean_up_small`` 会删除 ``OUTCAR*``
   和 ``CONTCAR*``。运行前确认已完成的数据已提取。先在非生产目录验证。

INCAR 操作
----------

pei_vasp_univ_find_and_change
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

批量修改 INCAR 参数。支持 comment/rewrite 两种模式，支持多参数对。

.. code-block:: bash

   # 修改 NBANDS 和 ALGO
   pei_vasp_univ_find_and_change -nbands 18 -algo Normal

   # 注释掉某个 tag（不改值）
   pei_vasp_univ_find_and_change -isif comment

   # 注释并改值
   pei_vasp_univ_find_and_change -isif comment:3

参数规则
^^^^^^^^

- ``-tag value``：将 ``tag`` 设为 ``value``
- ``-tag comment``：注释掉 ``tag`` 行
- ``-tag comment:value``：改值并保持注释状态
- 参数验证通过后才统一应用（从左到右）
- INCAR tag 自动转大写

脚本会查阅同目录的 ``incar_tag_help_vaspwiki.tsv`` 获取每个 tag 的 VASP Wiki
帮助说明。

pei_vasp_univ_increase_nbands
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

从 OUTCAR 读取 ``NELECT``、``NIONS`` 和当前 ``NBANDS``，按比例增加 NBANDS
并调用 ``pei_vasp_univ_find_and_change`` 写入 INCAR。

.. code-block:: bash

   # 默认读 OUTCAR
   pei_vasp_univ_increase_nbands

   # 指定 OUTCAR 路径
   pei_vasp_univ_increase_nbands /path/to/OUTCAR

POSCAR 操作
-----------

pei_vasp_univ_transfer_normal_to_selective
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

将 POSCAR 从普通模式转为选择性动力学模式，在第 8 行插入 ``Selective dynamics``，
并为每个原子行添加 ``F F T`` （z 方向自由，xy 固定）。

pei_vasp_univ_transfer_selective_to_normal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

反向操作：删除 ``Selective dynamics`` 行和所有 ``F``/``T`` 标记。

pei_vasp_univ_cp_contcar_cartesian_poscar
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

将 ``CONTCAR`` 复制为 ``POSCAR``，保持笛卡尔坐标格式。读取 CONTCAR 第二行
的缩放因子 ``a0``，重新缩放坐标。

收敛数据提取与绘图
------------------

pei_vasp_univ_extract_convergence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

从 ``y_dir/*/OUTCAR`` 提取每个离子步的收敛数据（``energy(sigma->0)`` 和
``max force norm``），写入 ``y_post_convergence/y_post_convergence_<subdir>.txt``。

支持 ``OUTCAR`` 和 ``OUTCAR.gz``。header 包含 ``natoms``、``EDIFFG``、``ISIF``
元数据。

.. code-block:: bash

   # 在包含 y_dir 的目录运行
   pei_vasp_univ_extract_convergence

pei_vasp_plot_convergence.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

读取 ``pei_vasp_univ_extract_convergence`` 产出的数据文件，绘制每个 ``y_dir``
子目录的弛豫收敛曲线（能量/原子 + 最大力）。

.. code-block:: bash

   python pei_vasp_plot_convergence.py

能量分量提取
------------

pei_vasp_univ_extract_energy_components
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

从 ``y_dir/*/OUTCAR`` 提取 VASP 能量分量（TOTEN、PSCENC、TEWEN、DENC、EXHF、
XCENC、Double counting、EENTRO、EBANDS、EATOM、Ediel_sol），写入
``p_post_energy_components.txt``。

.. code-block:: bash

   pei_vasp_univ_extract_energy_components

监控脚本
--------

pei_vasp_univ_monitor_error
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

扫描运行中作业的 stdout，检测致命错误关键字（``hermitian``、``sloshing``、
``DAV: 100``、``RMM: 100``），触发 ``scancel``。可选 ``-phase_check`` 模式
对每个目录跑 CNA 相变检查。

.. code-block:: bash

   pei_vasp_univ_monitor_error [-skip_ljobid JOBID,...] [-phase_check] [-cancel_on_phase_change]

pei_vasp_univ_monitor_slurm_state
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

提取 ``y_dir`` 所有子目录中 slurm 文件的最后 5 行，快速查看作业状态。

续算脚本
--------

pei_vasp_univ_resubmit_isif3
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

将 ``y_dir`` 中未收敛目录的 INCAR 改为 ``ISIF=3`` 并重新提交。

pei_vasp_univ_resubmit_isym0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

将 ``y_dir`` 中目录的 INCAR 改为 ``ISYM=0`` 并重新提交。

环境与结构信息
--------------

pei_vasp_univ_load_env
~~~~~~~~~~~~~~~~~~~~~~

VASP 运行环境加载脚本（**source 而非执行**）。加载 module、将 ``vasp_std``
加入 ``PATH``。可通过环境变量覆盖默认路径。

.. code-block:: bash

   source pei_vasp_univ_load_env

pei_vasp_univ_get_struct_infos
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Python 脚本，提取 ``y_dir`` 每个子目录的结构信息（晶格、厚度、原子数、CNA
相分布）和最终能量，写出两个固定宽度汇总表。

.. code-block:: bash

   python pei_vasp_univ_get_struct_infos [-root .] [-no_recursive]

pei_vasp_univ_check_phase_transition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Python CLI，对单个 VASP 目录跑 CNA 扫描，比较 POSCAR 和 CONTCAR 的相分布，
报告是否发生相变。设计为被 ``pei_vasp_univ_monitor_error`` 逐目录调用。

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - 退出码
     - 含义
   * - ``0``
     - 已检查，无相变
   * - ``1``
     - 检测到相变
   * - ``2``
     - 跳过（POSCAR/CONTCAR 缺失或为空）
   * - ``3``
     - 检查失败（结构不可读或 CNA 错误）

pei_vasp_univ_get_size_by_distance.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Python CLI，按 k 点距离 ``R_k`` 扫描计算 k 点网格 size。读取结构文件（默认
``CONTCAR``），调用 ``mymetal.calculate.calqm.kpoints.get_size_by_distance``。

.. code-block:: bash

   python pei_vasp_univ_get_size_by_distance.py \
       -input CONTCAR -output p_get_size_by_distance.txt \
       -rmin 20 -rmax 100 -rstep 1

.. seealso::

   :doc:`../tutorials/kpoints_sampling` 详细说明 ``R_k`` 与 k 点网格的关系。

示例脚本
--------

完整的示例脚本见：

.. literalinclude:: ../../examples/vasp_universal_overview.py
   :language: python
   :caption: vasp_universal 子包功能总览示例

运行方式：

.. code-block:: bash

   python docs/examples/vasp_universal_overview.py --output docs/_build/example-vasp-universal

安全检查清单
------------

1. 先 ``bash -n`` 检查脚本语法
2. 在非生产目录 dry-run
3. 清理脚本会删 ``OUTCAR*``——确认数据已提取
4. ``pei_vasp_univ_sbatch`` 会 ``cp CONTCAR POSCAR``——确认续算意图
5. 监控脚本的 ``-cancel_on_phase_change`` 会 ``scancel``——确认可接受

.. seealso::

   - :doc:`vasp` — VASP workflow 总览
   - :doc:`slurm` — Slurm 提交与 dry-run
   - :doc:`workflows` — Workflow 总览
   - :doc:`../reference/scripts` — 完整脚本参考
