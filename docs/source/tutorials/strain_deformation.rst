.. _tutorial-strain-deformation:

计算形变矩阵与 Lagrangian / Euler 应变张量
==========================================

:Audience: 想用 ``mymetal`` 把 cell 变形转换为 deformation gradient
   ``F``、Lagrangian 应变 ``E``、Euler 应变 ``e`` 和主应变的用户
:Time: 10 分钟
:Requires: Python、ASE、numpy、pjvasp_package
:Runs VASP: No

目标与边界
----------

演示 :mod:`mymetal.calculate.calmechanics` 中三个纯 numpy 工具：

1. :func:`cal_deform_matrix` -- 给定初始 cell ``A`` 与变形 cell ``B``
   （列矢量），返回 ``F = B · A^{-1}``。
2. :func:`cal_strain_matrix` -- 由 ``F`` 计算 Lagrangian 应变
   ``E = 1/2 (F^T F - I)`` 和 Euler 应变 ``e = 1/2 (I - (F F^T)^{-1})``。
3. :func:`cal_principal_and_shear_strain` -- 对 ``E`` 做特征分解，得到
   主应变（特征值）、主方向（特征向量）和偏应变（deviatoric）矩阵。

.. note::

   ``mymetal.calculate.calmechanics.strain`` 顶层 ``import hetbuilder``
   仅用于异质结专用 helper。本教程只调用上述三个 numpy 函数，因此
   ``strain_deformation.py`` 顶部用 ``unittest.mock`` 把 ``hetbuilder``
   替换为占位对象，使脚本在最小安装下也能运行。

公式与单位
----------

形变梯度 ``F`` 把初始 cell 列矢量 ``A`` 映射到变形 cell ``B``：

.. math::

   \mathbf{B} = \mathbf{F} \mathbf{A}, \quad
   \mathbf{F} = \mathbf{B} \mathbf{A}^{-1}.

Lagrangian 应变（小应变与大应变通用）：

.. math::

   \mathbf{E} = \frac{1}{2} \left(\mathbf{F}^\mathrm{T}\mathbf{F}
                                     - \mathbf{I}\right).

Euler 应变（大应变）：

.. math::

   \mathbf{e} = \frac{1}{2} \left(\mathbf{I}
                                  - \left(\mathbf{F}\mathbf{F}^\mathrm{T}\right)^{-1}\right).

主应变与主方向：

.. math::

   \mathbf{E} \mathbf{v}_i = \lambda_i \mathbf{v}_i, \quad
   \mathrm{diag}(\lambda_1, \lambda_2, \lambda_3).

本教程对每种变形只显示 Lagrangian ``E``；``E`` 与 ``e`` 在小应变下
近似相等。

输入与 provenance
-----------------

无外部文件依赖。Reference cell 由 ``ase.build.bulk('Cu', 'fcc', a=3.61,
cubic=True)`` 在运行时构建。三种变形分别为：

* **Uniaxial x**: ``F = diag(1+s, 1, 1)``，``s = 0.05``；
* **Biaxial xy**: ``F = diag(1+s, 1+s, 1)``，``s = 0.05``；
* **Simple shear xy**: ``F = I + γ · e_x ⊗ e_y``，``γ = 0.10``。

运行完整脚本
------------

.. code-block:: console

   $ python docs/examples/strain_deformation.py --output docs/_build/example-strain

.. literalinclude:: ../../examples/strain_deformation.py
   :language: python
   :linenos:

预期输出
--------

.. code-block:: text

   reference cell (A):
     [3.610000, 0.000000, 0.000000]
     [0.000000, 3.610000, 0.000000]
     [0.000000, 0.000000, 3.610000]

   Uniaxial x
   strain = +0.05
     deformation gradient F:
       [+1.050000, +0.000000, +0.000000]
       [+0.000000, +1.000000, +0.000000]
       [+0.000000, +0.000000, +1.000000]
     Lagrangian strain E:
       [+0.051250, +0.000000, +0.000000]
       [+0.000000, +0.000000, +0.000000]
       [+0.000000, +0.000000, +0.000000]
     principal strains (eigenvalues of E):
       [0.05125 0.      0.     ]

   Biaxial xy
   strain = +0.05
     deformation gradient F:
       [+1.050000, +0.000000, +0.000000]
       [+0.000000, +1.050000, +0.000000]
       [+0.000000, +0.000000, +1.000000]
     Lagrangian strain E:
       [+0.051250, +0.000000, +0.000000]
       [+0.000000, +0.051250, +0.000000]
       [+0.000000, +0.000000, +0.000000]
     principal strains (eigenvalues of E):
       [0.05125 0.05125 0.     ]

   Simple shear xy
   gamma = +0.10
     deformation gradient F:
       [+1.000000, +0.100000, +0.000000]
       [+0.000000, +1.000000, +0.000000]
       [+0.000000, +0.000000, +1.000000]
     Lagrangian strain E:
       [+0.000000, +0.050000, +0.000000]
       [+0.050000, +0.005000, +0.000000]
       [+0.000000, +0.000000, +0.000000]
     principal strains (eigenvalues of E):
       [-0.04756246  0.05256246  0.        ]
   wrote: .../docs/_build/example-strain/strain_deformation.png

.. figure:: /_static/images/generated/strain_deformation.png
   :alt: Lagrangian strain tensors for uniaxial x, biaxial xy and simple shear xy

   三个子图分别展示三种参考变形对应的 Lagrangian 应变张量 ``E``。
   颜色编码：红色为正应变（拉伸），蓝色为负应变（压缩）。对 Uniaxial x
   和 Biaxial xy，对角元 ``E_ii = 0.0513`` 等于 ``0.5*((1.05)^2 - 1)``，
   非对角元为零。对 Simple shear xy，非对角元 ``E_12 = E_21 = γ/2 = 0.05``，
   对角元 ``E_11 = γ^2/2 = 0.005`` 反映二阶几何效应。

结果含义
--------

* Uniaxial x 应变只在对角 ``E_11`` 上有非零值 ``0.0513``，对应工程
  应变 ``s = 0.05``；差异来自 Lagrangian 应变的二阶项
  ``0.5 * s^2 = 0.00125``。
* Biaxial xy 的 ``E_11 = E_22 = 0.0513``，``E_33 = 0``；如果加上
  Poisson 耦合，``E_33`` 会变为负值，本教程故意忽略。
* Simple shear 的 ``E_12 = γ/2 = 0.05`` 是教材定义；主应变为
  ``±sqrt((γ/2)^2 + (γ^2/2)^2) ≈ ±0.0526``，对应于沿 ±45° 方向的
  拉伸/压缩。
* 所有 ``E`` 矩阵都对称，符合应变张量定义。

参数说明
--------

.. list-table:: :func:`cal_deform_matrix` 关键参数
   :header-rows: 1
   :widths: 22 14 64

   * - 参数
     - 默认
     - 含义
   * - ``initial_cell``
     - 必填
     - 初始 cell，``3x3`` numpy 数组，列矢量（每列是一个 lattice
       vector）。
   * - ``deformed_cell``
     - 必填
     - 变形 cell，``3x3`` numpy 数组，列矢量。
   * - 返回值
     - -
     - ``F``，``3x3`` numpy 数组。``B = F · A``。

.. list-table:: :func:`cal_strain_matrix` 关键参数
   :header-rows: 1
   :widths: 22 14 64

   * - 参数
     - 默认
     - 含义
   * - ``deformation_matrix``
     - 必填
     - ``F``，``3x3`` numpy 数组。
   * - 返回值
     - -
     - ``[E_lagrangian, E_eulerian]``，两个 ``3x3`` 矩阵。

.. list-table:: :func:`cal_principal_and_shear_strain` 关键参数
   :header-rows: 1
   :widths: 22 14 64

   * - 参数
     - 默认
     - 含义
   * - ``strain_matrix``
     - 必填
     - ``3x3`` numpy 数组，通常是 ``E_lagrangian``。
   * - 返回值
     - -
     - ``(principal_matrix, shear_matrix, eigvals, eigvecs)``；
       ``principal_matrix`` 是 ``diag(eigvals)``，``shear_matrix`` 是
       偏应变（去对角部分）。

验证方法
--------

* Uniaxial x: ``E[0,0] = 0.5*((1+s)^2 - 1) = 0.05125``；其他对角元
  和全部非对角元为零。
* Biaxial xy: ``E[0,0] == E[1,1]``，且等于 Uniaxial 的 ``E[0,0]``。
* Simple shear: ``E[0,1] = E[1,0] = γ/2``；``E[0,0] = γ^2/2``。
* 三个 ``E`` 矩阵都对称（``1e-12`` 容差）。

常见错误
--------

``ValueError: Input matrix cannot be None``
   ``cal_strain_matrix`` 必须接收 ``F``，而不是 ``atoms`` 对象。先用
   ``cal_deform_matrix`` 计算 ``F``。

``ModuleNotFoundError: hetbuilder``
   ``mymetal.calculate.calmechanics.strain`` 顶层 ``import hetbuilder``。
   本教程用 ``unittest.mock`` 把 ``hetbuilder`` 替换为占位对象，仅当
   你不调用 ``cal_strain_matrix_root`` 等异质结 helper 时安全。

``E`` 不对称
   检查 ``F`` 是否正确。``F = B · A^{-1}`` 而不是 ``A^{-1} · B``；
   cell 必须以列矢量传入。

下一步
------

* :doc:`biaxial_stretch` 看如何对 slab 批量施加应变系列；
* :doc:`schmid_factor` 看如何把应变方向与 FCC 滑移系对应；
* :doc:`../manual/vasp` 了解 ``vasp_utils`` 中的 ``pei_vasp_run_stretch``
  workflow。

Related API
-----------

* :func:`mymetal.calculate.calmechanics.deformation.cal_deform_matrix`
* :func:`mymetal.calculate.calmechanics.strain.cal_strain_matrix`
* :func:`mymetal.calculate.calmechanics.strain.cal_principal_and_shear_strain`
* :func:`numpy.linalg.eig`
