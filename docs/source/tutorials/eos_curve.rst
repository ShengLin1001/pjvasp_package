.. _tutorial-eos-curve:

拟合 Murnaghan / Birch-Murnaghan 状态方程
==========================================

:Audience: 已有一组 (V, E) 数据点，想拟合体弹模量 ``B0``、平衡体积 ``V0``
   和压力导数 ``B0'`` 的用户
:Time: 10 分钟
:Requires: Python、ASE、scipy、pjvasp_package
:Runs VASP: No

目标与边界
----------

演示一个完整的 EOS 拟合工作流：从一组 ``(volume, energy)`` 数据点出发，
用 ASE :class:`ase.eos.EquationOfState` 引擎拟合 Murnaghan 与
Birch-Murnaghan 两种形式，并用 :func:`scipy.optimize.curve_fit` 在
Murnaghan 形式上额外提取压力导数 ``B0'``。最后用项目通用绘图样式渲染
双面板图（拟合曲线 + 残差）。

.. note::

   本教程的数据是合成的：用一组 Cu-like 参数（``V0 = 11.7 Å³``、
   ``B0 = 140 GPa``、``B0' = 5.2``）从 Murnaghan 方程反算能量，再加
   ``σ = 0.0015 eV`` 的高斯噪声。这样可以在没有真实 VASP 结果的情况下
   验证拟合是否正确恢复输入参数。脚本不会把合成数据伪装成真实 DFT
   benchmark。

公式与单位
----------

Murnaghan EOS（PRB 28, 5480, 1983）：

.. math::

   E(V) = E_0
   + \frac{B_0 V}{B_0'} \left[
     \frac{(V_0/V)^{B_0'}}{B_0' - 1} + 1
   \right]
   - \frac{B_0 V_0}{B_0' - 1}.

输入能量为 eV，体积为 Å³，``B0`` 在拟合内部以 eV/Å³ 表示，输出时转换为
GPa（``1 eV/Å³ = 160.21766 GPa``）。Birch-Murnaghan 形式由 ASE 内部实现，
本教程只调用 ``eos="birchmurnaghan"``。

输入与 provenance
-----------------

无外部文件依赖。所有数据由
:func:`docs.examples.eos_curve.make_synthetic_ev` 在运行时生成，种子固定为
``20260729``，因此终端输出和图片都是确定性的。

运行完整脚本
------------

.. code-block:: console

   $ python docs/examples/eos_curve.py --output docs/_build/example-eos

.. literalinclude:: ../../examples/eos_curve.py
   :language: python
   :linenos:

预期输出
--------

.. code-block:: text

   true:    V0=11.7000 Å^3  B0=140.00 GPa  B0'=5.200  E0=-8.200000 eV
   murn(ASE):    V0=11.6961 Å^3  B0=148.30 GPa            E0=-8.200076 eV
   bm  (ASE):    V0=11.6962 Å^3  B0=148.47 GPa            E0=-8.200079 eV
   murn(scipy):  V0=11.6961 Å^3  B0=148.30 GPa  B0'=3.224  E0=-8.200076 eV
   wrote: .../docs/_build/example-eos/eos_curve.png

.. figure:: /_static/images/generated/eos_curve.png
   :alt: Murnaghan and Birch-Murnaghan EOS fits to synthetic Cu-like data

   左图：合成数据点（黑点）、Murnaghan 拟合（蓝实线）、Birch-Murnaghan
   拟合（红虚线）和真实 ``V0`` （灰点线）。右图：两种拟合的残差。

结果含义
--------

* 两种 EOS 拟合的 ``V0`` 都在 ``11.696 Å³`` 左右，与真实 ``11.7 Å³`` 偏差
  小于 ``0.005 Å³``。这是 EOS 拟合最稳定的输出。
* ``B0`` 偏差约 ``8 GPa`` （≈ 6%）。这是 ``B0' = 5.2`` 远大于 ``4`` 时
  Murnaghan/BM 形式对噪声敏感的典型表现：在 ``B0' ≈ 4`` 附近两种形式接近
  重合，``B0'`` 越大两个形式对数据细节越敏感。
* ``B0'`` 的 scipy 拟合给出 ``3.224``，与真实 ``5.2`` 偏差较大。这是预期
  行为：``B0'`` 是 EOS 三阶项，对噪声最敏感，需要比 ``V0``/``B0`` 更密的
  数据点和更小的噪声才能稳定恢复。
* 残差图的纵轴量级为 ``1 meV``，与合成噪声 ``σ = 1.5 meV`` 一致，说明拟合
  没有引入额外偏差。

参数说明
--------

.. list-table:: :func:`make_synthetic_ev` 关键参数
   :header-rows: 1
   :widths: 22 14 64

   * - 参数
     - 默认
     - 含义
   * - ``e0``
     - ``-8.2``
     - 平衡能量 eV。
   * - ``v0``
     - ``11.7``
     - 平衡体积 Å³。
   * - ``b0_gpa``
     - ``140.0``
     - 体弹模量 GPa。
   * - ``b0p``
     - ``5.2``
     - 压力导数 ``B0'``。
   * - ``n_points``
     - ``13``
     - (V, E) 数据点数量。
   * - ``strain_amp``
     - ``0.06``
     - 体积应变半范围（``V/V0 - 1``）。
   * - ``noise_std_ev``
     - ``0.0015``
     - 高斯噪声标准差 eV。
   * - ``seed``
     - ``20260729``
     - RNG 种子，固定以保证可复现。

.. list-table:: :func:`fit_with_ase` / :func:`fit_murnaghan_with_b0p` 关键参数
   :header-rows: 1
   :widths: 22 14 64

   * - 参数
     - 默认
     - 含义
   * - ``eos`` (fit_with_ase)
     - ``"murnaghan"``
     - 传给 ASE 的 EOS 名字，可选 ``"murnaghan"``、``"birchmurnaghan"``、
       ``"sjeos"``、``"vinet"`` 等。
   * - ``p0`` (fit_murnaghan_with_b0p)
     - 自动
     - 初值由数据最小值点 + 教科书 ``B0=140 GPa``、``B0'=4`` 构造。

验证方法
--------

* ``abs(fit_murn["v0_A3"] - meta["v0_A3"]) < 0.05``：ASE Murnaghan 拟合
  的 ``V0`` 必须恢复输入值；
* ``abs(fit_bm["v0_A3"] - meta["v0_A3"]) < 0.05``：ASE Birch-Murnaghan
  拟合同样要求；
* ``abs(fit_murn["b0_gpa"] - meta["b0_gpa"]) < 25.0``：体弹模量允许较宽
  容差，因为 ``B0' = 5.2`` 时两种形式对噪声敏感；
* ``abs(fit_murn_b0p["b0p"] - meta["b0p"]) < 3.0``：scipy 拟合的 ``B0'``
  允许较大偏差，因为 ``B0'`` 本身难以稳定恢复。

常见错误
--------

``ValueError: No minimum!`` (ASE sjeos)
   数据的最小值落在端点而非中间。检查 ``strain_amp`` 是否过小、``n_points``
   是否过少，或是否 ``V0`` 落在数据范围外。

``RuntimeError: Optimal parameters not found`` (scipy)
   ``B0'`` 的初值离真值太远。``fit_murnaghan_with_b0p`` 默认用 ``B0'=4``，
   对 ``B0'`` 在 ``3-7`` 范围内的金属通常足够；若数据范围过窄，可手动给
   ``p0=[E_min, V_min, 100, 4]`` 重试。

``B0'`` 拟合结果为负或极大
   数据点太少或噪声太大。``B0'`` 是 EOS 三阶项，需要至少 9-11 个点且
   ``strain_amp`` 不小于 ``0.04`` 才能稳定。

下一步
------

* :doc:`biaxial_stretch` 看如何对结构施加单轴/双轴应变并合成能量；
* :doc:`surface_energy` 看如何从已有 VASP 结果计算表面能；
* :doc:`../manual/vasp` 了解 EOS workflow 脚本
  ``vasp_utils/vasp_workflow_bulk/pei_vasp_run_eos.py`` 的位置。

Related API
-----------

* :class:`ase.eos.EquationOfState`
* :func:`scipy.optimize.curve_fit`
* :func:`mymetal.universal.plot.general.general_set_all_rcParams`
