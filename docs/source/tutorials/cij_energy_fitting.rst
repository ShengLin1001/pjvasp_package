.. _tutorial-cij-energy-fitting:

用合成数据拟合二阶弹性常数 Cij
===============================

:Audience: 想理解 energy-strain 方法如何提取 Cij 的用户
:Time: 8 分钟
:Requires: Python、ASE、numpy、pjvasp_package
:Runs VASP: No

目标与边界
----------

演示 :func:`mymetal.post.Cij_energy.fit_cij_energy` 背后的 energy-strain
拟合数学。用合成 Cu-like 数据（C11=168, C12=121, C44=76 GPa）生成三组
应变-能量数据，做二次多项式拟合，提取 Cij 并与输入对比。

.. note::

   本教程的数据是合成的，不是真实 DFT 结果。``fit_cij_energy`` 函数本身
   需要真实 VASP 输出目录结构（``y_dir/<strain>``），这里用合成数据复现
   相同的拟合数学。脚本明确标注数据来源。

背景公式
--------

对立方晶系，energy-strain 方法对每种变形模式施加应变 ``eta``，记录能量
``U``，拟合：

.. math::

   \frac{U}{V_0} = \frac{1}{2} C_{ij}\, \eta_i\, \eta_j

三种独立变形模式：

* **uniaxial x**：``eta_xx`` 扫描，``U/V0 = 0.5*C11*eta²``
* **y-z shear**：``eta_yz`` 扫描，``U/V0 = 0.5*C44*eta²``
* **x-y biaxial**： ``eta_xx = eta_yy = eta`` 扫描，
  ``U/V0 = (C11+C12)*eta²`` （两个正应变各贡献 ``0.5*C11*eta²``，
  加上两个交叉项 ``0.5*C12*eta²``，总计 ``(C11+C12)*eta²``）

运行完整脚本
------------

.. code-block:: console

   $ python docs/examples/cij_energy_fitting.py --output docs/_build/example-cij

.. literalinclude:: ../../examples/cij_energy_fitting.py
   :language: python
   :linenos:

预期输出
--------

.. code-block:: text

   ================================================================
   Cij energy-strain fitting (synthetic cubic Cu-like data)
   reference volume V0 = 11.70 A^3
   strain grid = [-0.02, -0.01, 0.0, 0.01, 0.02]
   noise sigma = 1.0e-05 GPa (seed=42)
   ----------------------------------------------------------------
   mode                 fitted C (GPa)
   uniaxial_x           167.933
   yz_shear              75.946
   xy_biaxial           288.831  (= C11 + C12)
   ----------------------------------------------------------------
   extracted constants:
          input    fitted    rel.err(%)
          C11      168.000       167.933       0.0399
          C12      121.000       121.072       0.0593
          C44       76.000        75.946       0.0712

结果图
------

.. figure:: /_static/images/generated/cij_energy_fitting.png
   :alt: Strain-energy curves and input vs fitted Cij bar chart for synthetic Cu

   左图：三种变形模式的应变-能量曲线（数据点 + 二次拟合）。右图：输入
   vs 拟合的 C11/C12/C44 柱状图对比。

结果含义
--------

* 三种模式的二次拟合 ``R² ≈ 1``，拟合 Cij 与输入偏差 < 0.08%。
* ``xy_biaxial`` 模式拟合给出 ``C11+C12 = 288.8 GPa``，结合
  ``uniaxial_x`` 的 ``C11 = 167.9 GPa`` 可解出 ``C12 = 120.9 GPa``。
* 噪声 ``σ=1e-5 GPa`` 对拟合结果影响可忽略。

验证方法
--------

* 拟合 Cij 与输入偏差 < 5%；
* 图片非空白。

相关 API
--------

* :func:`mymetal.post.Cij_energy.fit_cij_energy`
* :func:`mymetal.post.Cij_energy.read_cij_energy`
* :func:`mymetal.post.Cij_energy.write_cij_energy`
* :doc:`strain_deformation` （变形 → Lagrangian 应变）
* :doc:`deformation_and_hnf` （变形矩阵）
