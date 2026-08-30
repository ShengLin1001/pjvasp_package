.. _tutorial-post-convergence:

ENCUT 与 KPOINTS 收敛测试后处理
===============================

:Audience: 想理解 convergence 后处理如何提取收敛曲线的用户
:Time: 8 分钟
:Requires: Python、numpy、pjvasp_package
:Runs VASP: No

目标与边界
----------

演示 :func:`mymetal.post.convergence.post_convergence` 的分析方法。使用
Au FCC 的真实 VASP 收敛测试结果，展示 ENCUT 和 KPOINTS 收敛曲线的生成、
排序、能量差计算和可接受的收敛判据。

.. note::

   本教程的数据来自 Au FCC 的 VASP 收敛计算（zcm6 集群，20230907_au_dft）。
   ``post_convergence`` 函数需要完整的 ``y_convergence/`` 目录结构。

背景公式
--------

收敛测试的目标是找到能量收敛到目标精度所需的最小参数：

* **ENCUT 收敛**：增大截断能直到 ``|E - E_conv| < tol`` (meV/atom)
* **KPOINTS 收敛**：增大 k 点网格直到能量收敛

能量差定义为：

.. math::

   \\Delta E(n) = E(n) - E_{\\text{conv}}

其中 ``E_conv`` 取最大参数下的能量（或外推值）。

真实数据结果
------------

.. list-table:: Au FCC ENCUT 收敛结果
   :header-rows: 1
   :widths: 20 25 25 30

   * - ENCUT (eV)
     - Etot (eV)
     - Delta E (meV)
     - 状态
   * - 230
     - -3.9197
     - 4.12
     - 未收敛
   * - 250
     - -3.9359
     - -12.07
     - 过收敛（振荡）
   * - 350
     - -3.9243
     - -0.37
     - 接近收敛
   * - 550
     - -3.9238
     - 0.13
     - 收敛基准

.. figure:: /_static/images/generated/post_convergence_encuts.png
   :alt: Au FCC ENCUT convergence: energy and delta-from-converged vs ENCUT

   Au FCC ENCUT 收敛曲线。上：总能量 vs ENCUT；下：与收敛基准的能量差
   (meV/atom)。橙色阴影标记 ``|Delta E| > 1 meV/atom`` 的未收敛区域。

.. figure:: /_static/images/generated/post_convergence_kpoints.png
   :alt: Au FCC KPOINTS convergence: energy vs kpoint grid size

   Au FCC KPOINTS 收敛曲线。能量 vs k 点网格 (NxNxN)。k 点收敛振荡
   较大，需要 15x15x15 以上才能稳定在 1 meV 以内。

结果含义
--------

* ENCUT 在 350 eV 后能量变化小于 1 meV/atom，可视为收敛。
* KPOINTS 收敛较慢，在 15x15x15 后振荡幅度 < 1 meV。
* 250 eV 处的"过收敛"是 ENCUT 不够时的数值振荡，不是物理效应。

相关 API
--------

* :func:`mymetal.post.convergence.post_convergence`
* :func:`mymetal.post.convergence.my_write_convergence`
* :func:`mymetal.post.convergence.my_read_convergence`
* :func:`mymetal.post.general.my_sort`

下一步
------

* :doc:`post_stretch_analysis` — 类似的后处理流程，用于拉伸计算；
* :doc:`post_cij_comparison` — 比较不同方法得到的弹性常数；
* :doc:`../api/post` — 完整的后处理 API 参考。
