.. _tutorial-post-kpar-ncore:

KPAR/NCORE 并行性能基准后处理
=============================

:Audience: 想理解 KPAR/NCORE 并行基准后处理的用户
:Time: 6 分钟
:Requires: Python、numpy、pjvasp_package
:Runs VASP: No

目标与边界
----------

演示 :func:`mymetal.post.kpar_ncore.post_kpar_ncore` 的分析方法。展示
KPAR/NCORE 参数对 VASP 计算时间和能量精度的影响，帮助选择最优并行配置。

背景
--------

VASP 的并行通过两个参数控制：

* **KPAR**：k 点并行组数。增大 KPAR 减少 MPI 通信，但每组可用核数减少。
* **NCORE**：每核处理的能带数。增大 NCORE 减少 MPI 通信，但 OpenMP
  线程数减少。

两者乘积不能超过总核数。本后处理扫描 KPAR * NCORE 组合，记录：

* **elapsed time** (s)：壁钟时间
* **energy(sigma->0)** (eV)：能量精度

合成数据演示
------------

.. figure:: /_static/images/generated/post_kpar_ncore_demo.png
   :alt: KPAR/NCORE benchmark: time and energy delta vs NCORE grouped by KPAR

   上：elapsed time vs NCORE，按 KPAR 分组。下：delta energy vs NCORE
   (meV/atom)。最优组合在时间最短且能量偏差可接受处。

参数说明
--------

.. list-table:: ``post_kpar_ncore`` 参数
   :header-rows: 1
   :widths: 25 15 60

   * - 参数
     - 默认
     - 含义
   * - ``path_workflow``
     - ``None``
     - 包含 ``y_dir`` 的工作流目录
   * - ``run_post``
     - ``True``
     - 是否先运行 ``pei_vasp_univ_post``
   * - ``save_fig_path``
     - ``None``
     - 图表保存路径
   * - ``save_txt_path``
     - ``None``
     - 结果表保存路径

KPAR 扫描范围
--------------

默认扫描的 KPAR/NCORE 组合：

.. list-table:: 默认基准网格
   :header-rows: 1
   :widths: 15 85

   * - KPAR
     - NCORE 值
   * - 128
     - 1
   * - 64
     - 1, 2
   * - 32
     - 1, 2, 4
   * - 16
     - 1, 2, 4, 8
   * - 8
     - 1, 2, 4, 8, 16
   * - 4
     - 1, 2, 4, 8, 16, 32

相关 API
--------

* :func:`mymetal.post.kpar_ncore.post_kpar_ncore`
* :func:`mymetal.post.kpar_ncore.read_kpar_ncore_times`
* :func:`mymetal.post.kpar_ncore.read_kpar_ncore_energies`
* :func:`mymetal.post.kpar_ncore.write_kpar_ncore_times`
