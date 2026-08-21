.. _tutorial-post-relax-convergence:

离子弛豫收敛轨迹后处理
=======================

:Audience: 想理解 ionic relaxation 收敛轨迹如何提取和可视化的用户
:Time: 6 分钟
:Requires: Python、numpy、pjvasp_package
:Runs VASP: No

目标与边界
----------

演示 :func:`mymetal.post.relax_convergence.my_univ_post_convergence` 的
分析方法。使用真实 VASP 弛豫数据，展示 ionic step 的能量和力收敛轨迹。

.. note::

   本教程的数据来自 VASP 弛豫计算（zcm6 集群，GSFE slab 弛豫）。
   ``my_univ_post_convergence`` 读取 ``pei_vasp_univ_extract_convergence``
   生成的数据文件。

背景
--------

VASP 的 ionic relaxation 逐步降低能量和力，直到满足 EDIFFG 判据：

* **EDIFF < 0**：能量收敛判据 (eV)
* **EDIFFG < 0**：力收敛判据 (eV/Ang)
* **EDIFFG > 0**：能量收敛判据 (eV)

本后处理追踪每个 ionic step 的：

* ``energy(sigma->0)`` 相对于最后一帧的差值 (meV/atom)
* 最大力 norm (eV/Ang)，与 ``|EDIFFG|`` 力判据对比

真实数据结果
------------

.. figure:: /_static/images/generated/post_relax_convergence.png
   :alt: Ionic relaxation convergence: energy and force vs step on log scale

   上：``|E - E_final|`` (meV/atom) vs ionic step，对数纵轴。
   下：``|F_max|`` (eV/Ang) vs ionic step，对数纵轴。红色虚线为 ``|EDIFFG|``
   力判据。两个面板都取绝对值，避免大值遮蔽近收敛细节。

参数说明
--------

.. list-table:: ``my_univ_post_convergence`` 参数
   :header-rows: 1
   :widths: 25 15 60

   * - 参数
     - 默认
     - 含义
   * - ``workflow``
     - ``'y_post_convergence'``
     - 数据文件所在目录
   * - ``yscale``
     - ``'log'``
     - 纵轴刻度（``'log'`` 或 ``'linear'``）

数据文件格式
------------

输入数据文件由 bash 脚本 ``pei_vasp_univ_extract_convergence`` 生成，格式：

.. code-block:: text

   # convergence data extracted from OUTCAR
   # natoms 12
   # ediffg -.1E-02
   # isif 3
      frame   energy_sigma0(eV)   max_force_norm(eV/A)
       1      -45.64932601          0.00000000
       2      -45.64935928          0.00000000

相关 API
--------

* :func:`mymetal.post.relax_convergence.my_univ_post_convergence`
* :func:`mymetal.post.relax_convergence.read_univ_post_convergence`
