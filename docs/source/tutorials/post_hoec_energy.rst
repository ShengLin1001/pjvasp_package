.. _tutorial-post-hoec-energy:

高阶弹性常数 (HOEC) 拟合
========================

:Audience: 想理解 energy-strain 方法如何提取 2/3/4 阶弹性常数的用户
:Time: 12 分钟
:Requires: Python、numpy、pjvasp_package
:Runs VASP: No

目标与边界
----------

演示 :func:`mymetal.post.hoec_energy.post_hoec_energy` 背后的高阶弹性常数
拟合数学。使用 Au HCP 的真实 DFT 计算结果（来自 zcm6 集群 A11-2 工作流），
展示 19 种变形模式的能量-应变曲线、多项式拟合、以及 2/3/4 阶常数的求解。

.. note::

   本教程的数据来自真实的 Au HCP VASP 计算。``post_hoec_energy`` 函数需要
   完整的 ``y_hoec_energy/`` 目录结构，这里用已提取的结果数据复现相同的
   分析流程。脚本明确标注数据来源。

背景公式
--------

对六方晶系，energy-strain 方法对每种变形模式 ``d`` 施加应变 ``xi``，
记录应变能密度 ``u(xi) = [U(xi) - U(0)] / V0`` (GPa)，拟合多项式：

.. math::

   u(\\xi) = P_2\\,\\xi^2 + P_3\\,\\xi^3 + P_4\\,\\xi^4 + \\ldots

各阶系数 ``P_n`` 组合后构成线性系统：

.. math::

   A_2\\,@\\,\\text{SOEC} = P_2, \\quad A_3\\,@\\,\\text{TOEC} = P_3, \\quad
   A_4\\,@\\,\\text{FOEC} = P_4

六方晶系有 5 个 SOEC (C11, C12, C13, C33, C44)、10 个 TOEC、19 个 FOEC。
系统通过最小二乘法求解，使用精选择的典型模式子集。

运行完整脚本
------------

.. code-block:: console

   $ python docs/examples/post_hoec_energy_demo.py --output docs/_build/example-hoec

真实数据结果
------------

以下结果来自 Au HCP 的 VASP 计算（zcm6 集群，A11-2 工作流）：

.. list-table:: Au HCP 高阶弹性常数 (GPa)
   :header-rows: 1
   :widths: 20 20 20 20 20

   * - 阶数
     - C11
     - C12
     - C13
     - C33
   * - 2nd (SOEC)
     - 229.31
     - 173.32
     - 142.27
     - 232.96
   * -
     -
     -
     -
     -
   * - 阶数
     - C111
     - C112
     - C113
     - C123
   * - 3rd (TOEC)
     - -4804.20
     - -461.51
     - -814.76
     - -530.46

.. figure:: /_static/images/generated/post_hoec_energy_au.png
   :alt: Au HCP HOEC 2nd/3rd/4th order elastic constants bar chart

   Au HCP 高阶弹性常数。左：2 阶 (SOEC) 5 个常数；中：3 阶 (TOEC) 10 个常数；
   右：4 阶 (FOEC) 19 个常数。4 阶常数数值较大，使用对数刻度展示。

模式拟合与交叉验证
------------------

HOEC 工作流使用精选择的典型模式子集求解各阶常数，剩余模式作为独立验证：

.. list-table:: 各阶求解模式与交叉验证残差
   :header-rows: 1
   :widths: 15 20 20 45

   * - 阶数
     - 求解模式数
     - 验证模式数
     - 验证残差 RMS (GPa)
   * - 2nd
     - 5 (M01-M03,M05,M07)
     - 15
     - 10.75
   * - 3rd
     - 10
     - 10
     - 284.32
   * - 4th
     - 19
     - 1 (M06)
     - 200.32

.. figure:: /_static/images/generated/post_hoec_mode_fits.png
   :alt: Per-mode energy-strain curves with polynomial fits for Au HCP HOEC

   各变形模式的能量-应变曲线（数据点 + 多项式拟合）。每个模式在
   ``|xi| <= 0.12`` 窗口内做 4 阶多项式拟合。

结果含义
--------

* Au HCP 的 ``C44 = 22.68 GPa`` 远小于 ``C11 = 229.31 GPa``，反映基面
  滑移方向的低剪切阻力。
* 3 阶常数 ``C111 = -4804 GPa`` 量级远大于 2 阶，反映了非简谐效应的
  强烈非对称性。
* 4 阶常数 ``C1111 = 77087 GPa`` 进一步约束了大应变下的能量响应，
  对塑性变形模拟至关重要。
* 交叉验证残差在 2 阶最小 (10.75 GPa)，说明 SOEC 求解最稳定；高阶
  残差增大反映了模式间响应的波动。

参数说明
--------

.. list-table:: ``post_hoec_energy`` 关键参数
   :header-rows: 1
   :widths: 25 15 60

   * - 参数
     - 默认
     - 含义
   * - ``fitmax``
     - ``0.12``
     - 拟合窗口上限 ``|xi| <= fitmax``
   * - ``fitdeg``
     - ``4``
     - 多项式拟合阶数（Wang-Li 方法用 4 阶）
   * - ``maxorder``
     - ``4``
     - 求解的最高弹性常数阶数
   * - ``fix_soec``
     - ``None``
     - 从 Cij-energy 导入的 2 阶常数（固定 P2）
   * - ``solve_modes``
     - ``None``
     - 显式指定求解模式子集

验证方法
--------

* 拟合多项式 ``R²`` 接近 1；
* 交叉验证残差在合理范围内；
* 图片非空白。

相关 API
--------

* :func:`mymetal.post.hoec_energy.post_hoec_energy`
* :func:`mymetal.post.hoec_energy.fit_P`
* :func:`mymetal.post.hoec_energy.solve_constants`
* :func:`mymetal.post.hoec_energy.write_hoec_energy`
* :func:`mymetal.post.hoec_energy.read_hoec_energy`
* :doc:`cij_energy_fitting` （2 阶弹性常数基础）
