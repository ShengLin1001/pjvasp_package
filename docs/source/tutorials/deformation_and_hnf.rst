.. _tutorial-deformation-and-hnf:

变形矩阵与 Hermite Normal Form
===============================

:Audience: 需要计算变形梯度或分析超胞变换矩阵的用户
:Time: 6 分钟
:Requires: Python、ASE、numpy、pjvasp_package
:Runs VASP: No

目标与背景
----------

演示两个底层矩阵工具：

1. :func:`mymetal.calculate.calmechanics.deformation.cal_deform_matrix` —
   从初始 cell 和变形后 cell 计算变形梯度 ``F = B @ A⁻¹``；
2. :func:`mymetal.calculate.calmath.matrix.hermite_normal_form` —
   计算整数矩阵的 Hermite Normal Form (HNF)，用于判断两个超胞是否
   commensurate、寻找最小公倍超胞。

变形梯度是应变分析（ :doc:`strain_deformation` ）和弹性常数拟合
（Cij/HOEC）的基础。HNF 在异质结构匹配中用于把两个不同超胞映射到
共同 cell。

运行完整脚本
------------

.. code-block:: console

   $ python docs/examples/deformation_and_hnf.py --output docs/_build/example-deform-hnf

.. literalinclude:: ../../examples/deformation_and_hnf.py
   :language: python
   :linenos:

预期输出
--------

.. code-block:: text

   == cal_deform_matrix ==
   initial_cell (Å):
   [[3.61 0.   0.  ]
    [0.   3.61 0.  ]
    [0.   0.   3.61]]
   deformed_cell (Å):
   [[3.7905 0.     0.    ]
    [0.     3.61   0.    ]
    [0.     0.     3.61  ]]
   deformation gradient F = deformed_cell @ inv(initial_cell):
   [[1.05 0.   0.  ]
    [0.   1.   0.  ]
    [0.   0.   1.  ]]

   == hermite_normal_form ==
   matrix                          det   original -> HNF
   M_a (upper-triangular)           12   [[2 1 0] [0 2 1] [0 0 3]]  ->  [[2 1 0] [0 2 1] [0 0 3]]
   M_b (unit upper-triangular)       1   [[1 2 3] [0 1 4] [0 0 1]]  ->  [[1 0 0] [0 1 0] [0 0 1]]
   M_c (general dense)              -2   [[3 1 2] [1 2 1] [2 1 1]]  ->  [[-1 -1 0] [-1 1 0] [7 0 -12]]

结果图
------

.. figure:: /_static/images/generated/deformation_and_hnf.png
   :alt: Deformation gradient schematic and Hermite Normal Form heatmaps

   左图：FCC Cu 初始 cell（蓝）与 x 方向 5% 拉伸后 cell（红）叠加示意。
   中图/右图：一般稠密整数矩阵 ``M_c`` 变换前后的热图，每个格子标注数值。

结果含义
--------

* **变形梯度**：``F[0,0] = 1.05``，其余对角元为 1.0，非对角元为 0，
  正好对应 x 方向 5% 单轴拉伸。
* **HNF M_a**：上三角矩阵的 HNF 是自身（已满足 HNF 条件）。
* **HNF M_b**：单位上三角矩阵的 HNF 为单位矩阵（det=1，可逆整数矩阵的
  HNF 是单位阵）。
* **HNF M_c**：一般稠密矩阵的 HNF 把非对角元消为零、对角元化为特定形式，
  det 不变（-2）。

验证方法
--------

脚本执行以下断言：

* ``F[0,0] ≈ 1.05``，``F[1,1] = F[2,2] = 1.0``，非对角元为 0；
* HNF 的行列式等于原矩阵行列式；
* 图片非空白。

相关 API
--------

* :func:`mymetal.calculate.calmechanics.deformation.cal_deform_matrix`
* :func:`mymetal.calculate.calmath.matrix.hermite_normal_form`
* :doc:`strain_deformation` （变形 → Lagrangian 应变的完整流程）
