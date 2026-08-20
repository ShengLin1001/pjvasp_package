.. _tutorial-reciprocal-lattice:

倒格子矢量与 RK 自动 k 点网格
=============================

:Audience: 需要计算倒格子或根据 RK/KSPACING 生成 k 点网格的用户
:Time: 6 分钟
:Requires: Python、ASE、numpy、pjvasp_package
:Runs VASP: No

目标与背景
----------

演示 :mod:`mymetal.calculate.calqm.kpoints` 中两个常被忽略但很实用的函数：

1. :func:`mymetal.calculate.calqm.kpoints.cal_reciprocal_matrix` —
   叉积法计算倒格子矢量；
2. :func:`mymetal.calculate.calqm.kpoints.cal_reciprocal_matrix2` —
   矩阵求逆法计算倒格子矢量；
3. :func:`mymetal.calculate.calqm.kpoints.get_size_by_distance` —
   根据 RK 乘积或 KSPACING 计算 VASP 自动 k 点网格。

本教程与 :doc:`kpoints_sampling` 互补：后者展示 MP/Gamma k 点位置和 RK 扫描，
本教程聚焦倒格子矢量本身和 RK → mesh 的映射。

.. note::

   ``cal_reciprocal_matrix`` 和 ``cal_reciprocal_matrix2`` 两种方法在数学上
   等价，结果在数值精度内一致。叉积法更直观，矩阵求逆法更简洁。

运行完整脚本
------------

.. code-block:: console

   $ python docs/examples/reciprocal_lattice.py --output docs/_build/example-reciprocal

.. literalinclude:: ../../examples/reciprocal_lattice.py
   :language: python
   :linenos:

预期输出
--------

.. code-block:: text

   label    formula  |b1|     |b2|     |b3|     old_mesh        new_mesh
   FCC Cu   Cu4      1.7405  1.7405  1.7405  [22, 22, 22]  [23, 23, 23]
   BCC Fe   Fe2      2.1893  2.1893  2.1893  [28, 28, 28]  [28, 28, 28]
   BCC Fe   Fe2      2.1893  2.1893  2.1893  [28, 28, 28]  [28, 28, 28]
   HCP Mg   Mg2      2.2602  2.2602  1.2060  [29, 29, 15]  [29, 29, 16]

结果图
------

.. figure:: /_static/images/generated/reciprocal_lattice.png
   :alt: Reciprocal lattice vector lengths and RK-based k-point meshes for FCC/BCC/HCP

   左图：FCC Cu、BCC Fe、HCP Mg 三种结构的 ``|b1|``、``|b2|``、``|b3|``
   倒格子矢量长度对比（Å⁻¹）。右图：RK=80 时 ``get_size_by_distance``
   给出的 old（round）和 new（ceil）k 点网格对比。

结果含义
--------

* **倒格子矢量长度**：FCC Cu ``a=3.61 Å`` → ``|b| = 2π/a × √3 ≈ 1.74 Å⁻¹``；
  HCP Mg 的 ``|b3|`` 明显小于 ``|b1|``/``|b2|``，因为 ``c > a``。
* **RK → mesh**：RK=80 时，FCC Cu 立方 cell 给出 ``[22,22,22]``（old）/
  ``[23,23,23]``（new）；HCP Mg 因 ``c/a ≈ 1.62``，z 方向网格更稀
  （``[29,29,15]``）。old 法用 ``round``，new 法用 ``ceil``，所以 new
  网格总是 ≥ old 网格。

验证方法
--------

脚本执行以下断言：

* ``cal_reciprocal_matrix`` 和 ``cal_reciprocal_matrix2`` 结果 ``allclose``；
* k 点网格为正整数；
* 图片非空白。

相关 API
--------

* :func:`mymetal.calculate.calqm.kpoints.cal_reciprocal_matrix`
* :func:`mymetal.calculate.calqm.kpoints.cal_reciprocal_matrix2`
* :func:`mymetal.calculate.calqm.kpoints.get_size_by_distance`
* :func:`mymetal.calculate.calqm.kpoints.get_kpoints_by_size`
* :doc:`kpoints_sampling` （MP/Gamma k 点位置与 RK 扫描）
