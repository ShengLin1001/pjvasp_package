.. _tutorial-post-cij-comparison:

二阶弹性常数 Cij: VASP vs LAMMPS 对比
=====================================

:Audience: 想对比 DFT 和 NNP 势的弹性常数结果的用户
:Time: 6 分钟
:Requires: Python、numpy、pjvasp_package
:Runs VASP: No

目标与边界
----------

演示 :func:`mymetal.post.Cij_energy.post_lammps_Cij_energy` 和
:func:`mymetal.post.Cij_energy.fit_cij_energy` 的分析方法。使用 Au HCP (VASP)
和 FCC (LAMMPS UNNEP) 的真实计算结果，对比二阶弹性常数。

.. note::

   VASP 数据来自 zcm6 集群 A11-2 工作流 (Au HCP)，
   LAMMPS 数据来自 ``pj-test-properties-gold/y_Cij_energy/fcc``。

真实数据结果
------------

.. list-table:: Cij 对比 (GPa)
   :header-rows: 1
   :widths: 15 25 25 35

   * - 常数
     - VASP Au HCP
     - LAMMPS FCC (UNNEP)
     - 说明
   * - C11
     - 259.12
     - 154.72
     - 面内刚度
   * - C12
     - 181.82
     - 118.32
     - 面内 Poisson 耦合
   * - C13
     - 142.92
     - 118.32
     - 面内-面外耦合
   * - C33
     - 242.75
     - 154.72
     - c 轴刚度
   * - C44
     - 20.50
     - 32.69
     - 剪切模量

.. list-table:: 导出弹性性质
   :header-rows: 1
   :widths: 20 25 25 30

   * - 性质
     - VASP Au HCP
     - LAMMPS FCC
     - 公式
   * - E_x (GPa)
     - 120.45
     - 52.17
     - Young's modulus (x)
   * - E_z (GPa)
     - 150.10
     - 52.17
     - Young's modulus (z)
   * - nu_xy
     - 0.558
     - 0.433
     - Poisson ratio
   * - mu_xz (GPa)
     - 20.50
     - 32.69
     - = C44

.. figure:: /_static/images/generated/post_cij_comparison.png
   :alt: Cij bar chart comparing VASP Au HCP and LAMMPS FCC

   VASP Au HCP vs LAMMPS FCC (UNNEP) 二阶弹性常数对比柱状图。
   HCP 的 C11/C33 大于 FCC，反映 Au HCP 的更高对称性约束。

结果含义
--------

* Au HCP 的 ``C12/C11 = 0.70``，Poisson 比高达 0.558，说明 HCP Au 的
  面内横向耦合很强。
* LAMMPS FCC (UNNEP) 的 ``C12/C11 = 0.77``，Poisson 比 0.433，在 FCC
  金属典型范围 (0.3-0.4) 内。
* Au HCP 的 ``C44 = 20.50 GPa`` 小于 FCC 的 32.69 GPa，说明 HCP 基面
  滑移方向更容易剪切。

相关 API
--------

* :func:`mymetal.post.Cij_energy.post_lammps_Cij_energy`
* :func:`mymetal.post.Cij_energy.fit_cij_energy`
* :func:`mymetal.post.Cij_energy.write_cij_energy`
* :func:`mymetal.post.Cij_energy.read_cij_energy`
* :doc:`cij_energy_fitting` （合成数据演示）
* :doc:`post_hoec_energy` （高阶弹性常数）
