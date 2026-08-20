Workflow 总览
=============

本页是 VASP、LAMMPS、n2p2 三个 workflow 的导航总览。每个 workflow 的详细
说明（应用场景、输入文件、执行流程、真实结果展示）见对应的子页面。

.. contents:: 本页内容
   :local:
   :depth: 2

组件边界
--------

``mymetal``
   可导入的 Python package：结构、I/O、计算、后处理与 ML 数据工具。

``vasp_utils/``、``slurm_utils/``、``lmp_utils/``、``n2p2_utils/``
   目录型 workflow 与 HPC 脚本；不是 Python package 的 console entry point。

VASP、LAMMPS、n2p2、Slurm、MPI 与 module
   集群提供的外部运行环境，本仓库不安装它们。

``mymetal/example/``
   用于演示和解析验证的 tracked fixture；不是生产计算模板。

推荐生命周期
------------

.. code-block:: text

   structure
      -> source inputs
      -> generated case directories
      -> reviewed submit scripts
      -> queued/running jobs
      -> completed or failed outputs
      -> parsed tables and figures

每一箭头都应有一个可观察检查：原子数、输入 diff、``bash -n``、Slurm state、
``OUTCAR`` convergence 或结果断言。不要把“脚本已生成”写成“作业已提交”，也不要
把“作业结束”自动等同于“电子或离子收敛”。

常见目录
--------

.. code-block:: text

   project/
   |-- y_full_relax/
   |-- y_src/
   |   |-- INCAR
   |   |-- KPOINTS
   |   |-- POTCAR
   |   `-- poscars2/
   `-- y_dir/
       |-- 0.997/
       |-- 1.000/
       `-- 1.003/

``y_src`` 保存源输入，``y_dir/<case>`` 保存独立计算。脚本可能假定这一命名；
运行前以实际源码和 ``-h`` 为准。

三个 workflow 的协作
--------------------

VASP、LAMMPS 和 n2p2 在一个典型的材料势函数研究中形成流水线：

.. code-block:: text

   ┌─────────────┐     ┌──────────────┐     ┌─────────────┐
   │  VASP       │     │  n2p2 训练    │     │  LAMMPS     │
   │  (DFT 参考) │ ──> │  (NNP 拟合)  │ ──> │  (性质验证) │
   └─────────────┘     └──────────────┘     └─────────────┘
        │                     │                     │
        ▼                     ▼                     ▼
   OUTCAR → input.data   nnp-train → weights   pair_style hdnnp
   (能量/力/结构)         (learning-curve)      (stretch/Cij/GSFE)

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - 阶段
     - 工具
     - 产出
   * - DFT 参考
     - VASP workflow (:doc:`vasp`)
     - ``OUTCAR`` （能量、力、应力、收敛状态）
   * - 数据准备
     - ``mymetal.ml.n2p2.dataset``
     - ``input.data`` （n2p2 格式的结构 + 能量 + 力）
   * - NNP 训练
     - n2p2 workflow (:doc:`n2p2`)
     - ``weights.*.out`` 、 ``scaling.data`` 、 ``learning-curve.out``
   * - 性质验证
     - LAMMPS workflow (:doc:`lammps`)
     - ``p_post_stretch.*`` 、 ``y_post_cij_energy.*`` 、 ``y_post_gsfe.*``
   * - 提交与调度
     - Slurm workflow (:doc:`slurm`)
     - ``sub_slurm_univ.sh`` （统一入口，三 workflow 共用）

真实计算案例
------------

本仓库文档中的真实结果图表基于 Au HCP/FCC 的集群计算，数据源如下：

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - workflow
     - 数据路径
     - 关键产出
   * - VASP
     - ``construct_dataset/calculate/A11-1/``
     - stretch: a₀=2.85925 Å, poly coeffs
   * - VASP
     - ``construct_dataset/calculate/A11-2/``
     - Cij (5 常数), HOEC (2/3/4 阶)
   * - VASP
     - ``construct_dataset/calculate/A11-3/``
     - cohesive: E_coh=-3.867 eV/atom, B=188.17 GPa
   * - VASP
     - ``construct_dataset/calculate/decohesion/``
     - γ∞ ≈ 1969 mJ/m²
   * - LAMMPS
     - ``pj-test-properties-gold/y_stretch/fcc/``
     - FCC a₀=4.15852 Å (UNNEP 势)
   * - LAMMPS
     - ``pj-test-properties-gold/y_Cij_energy/fcc/``
     - C11=154.72, C44=32.69 GPa
   * - LAMMPS
     - ``pj-test-properties-gold/y_gsfe/fcc/FCC_111/``
     - γ_USFE ≈ 106.5 mJ/m²
   * - n2p2
     - ``train/stage0/0/``
     - E RMSE 6.5 meV/atom, F RMSE 31 meV/Å (1500 epochs)

按任务选择路径
--------------

* 只构建结构：从 :doc:`../getting_started/au111_slab` 开始。
* 计算已有能量的 surface energy：使用 :doc:`../tutorials/surface_energy`。
* 汇总现有 ``OUTCAR``：使用 :doc:`../tutorials/outcar_batch`。
* 生成或提交 HPC 作业：先读 :doc:`slurm`，再读 :doc:`vasp`。
* 查参数和单位：使用 :doc:`../api`；查脚本位置：使用
  :doc:`../reference/scripts`。

.. seealso::

   - :doc:`vasp` — VASP workflow 图文总览（含真实结果展示）
   - :doc:`lammps` — LAMMPS workflow 图文总览（含真实结果展示）
   - :doc:`n2p2` — n2p2 workflow 图文总览（含真实结果展示）
   - :doc:`slurm` — Slurm 提交与 dry-run
