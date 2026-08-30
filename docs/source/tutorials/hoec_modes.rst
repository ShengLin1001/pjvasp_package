.. _tut-hoec-modes:

Higher-order elastic constants: mode table and deformation gradient
====================================================================

目标
----

演示 ``mymetal.calculate.calmechanics.hoec`` 模块的核心功能：
运行内置 self-test、查看立方和六方对称性的独立弹性常数集合和变形模式表，
并用变形梯度验证 Lagrangian 应变与 Voigt 方向的关系。

前置条件
--------

* 已安装 ``mymetal`` 和 ``numpy``；
* 了解 Voigt 记号、Green-Lagrange 应变和 Brugger 弹性常数。

物理问题
--------

高阶弹性常数（HOEC）用 energy-strain 方法拟合：对每个变形模式
``eta = xi * d``，应变能密度为多项式
``[U(xi)-U(0)]/V0 = P2*xi² + P3*xi³ + P4*xi⁴``。点群对称性把完全对称
的 Brugger 常数缩减为独立集合（立方：3/6/11，六方：5/10/19），每个
``P_n`` 是独立常数的固定线性组合。

此模块基于 Wang & Li (PRB 79, 224102, 2009)，推广到任意点群。

核心代码
--------

.. code-block:: python

   from mymetal.calculate.calmechanics.hoec import (
       get_model, get_hoec_modes, get_strain_list,
       get_deformation_gradient, check_symmetry, selftest_hoec,
   )

   # 1. Self-test: verify cubic against Wang-Li Table I
   selftest_hoec()

   # 2. Cubic mode table
   model = get_model("cubic")
   modes = get_hoec_modes("cubic")
   print(model.names(2))  # SOEC: ['11', '12', '44']
   print(model.names(3))  # TOEC: ['111', '112', '123', '144', '155', '456']
   print(model.names(4))  # FOEC: 11 constants

   # 3. Strain list and deformation gradient
   xi_list = get_strain_list(emax=0.06, de=0.02)
   d_A = modes["A"]  # uniaxial x: d = (1, 0, 0, 0, 0, 0)
   F = get_deformation_gradient(d_A, xi=0.02)

Convention
----------

* ``eta`` 是 engineering-Voigt 应变（``eta4 = 2*eps_yz`` 等，即工程剪应变
  带因子 2）；
* ``d`` 是 6 分量 Voigt 方向，与 ``eta`` 同记号；
* ``F`` 是对称变形梯度，满足 ``E = 1/2(F^T F - I) = xi * d_tensor``，
  其中 ``d_tensor`` 把 Voigt 方向转换为对称张量（剪切项除以 2）；
* ``get_strain_list`` 返回对称的 ``xi`` 列表（含 0），用于多项式拟合。

完整可运行脚本
--------------

.. literalinclude:: ../../examples/hoec_modes.py
   :language: python

运行：

.. code-block:: shell

   python docs/examples/hoec_modes.py --output docs/_build/example-hoec

输出检查
--------

* self-test 报告 cubic 和 hex 的 rank 均满秩（3/3, 6/6, 11/11, 5/5, 10/10,
  19/19）；
* Cubic model 有 3 SOEC、6 TOEC、11 FOEC；
* Hex model 有 5 SOEC、10 TOEC、19 FOEC；
* Cubic 模式表有 11 个模式（A..K）；
* Hex 模式表有 23 个模式（M01..M23）；
* Mode A（单轴 x）的 ``E`` 等于 ``xi * diag(1, 0, 0)``，在 1e-10 内；
* Mode F（纯剪切）的 ``E`` 剪切项等于 ``xi * d / 2 = 0.03``，在 1e-10 内；
* ``check_symmetry`` 从 Cu cell 正确检测出 cubic。

相关 API
--------

* :mod:`mymetal.calculate.calmechanics.hoec`
* :doc:`../api/calculate` — HOEC API reference

下一步
------

* :doc:`strain_deformation` — 变形梯度与 Lagrangian/Euler 应变的基础概念；
* :doc:`post_hoec_energy` — 用真实 VASP 能量数据拟合 HOEC；
* :doc:`cij_energy_fitting` — 二阶弹性常数的 energy-strain 拟合。
