.. _tutorial-biaxial-stretch:

对 Cu(111) slab 施加单轴/双轴应变并合成能量
============================================

:Audience: 想用 ``stretch_list_along_direction_to_cell`` 批量生成应变系列、
   并理解不同方向应变对应不同能量响应的用户
:Time: 10 分钟
:Requires: Python、ASE、numpy、pjvasp_package
:Runs VASP: No

目标与边界
----------

演示 :func:`mymetal.build.film.stretch.stretch_list_along_direction_to_cell`
的三个典型用法：沿 ``x`` 单轴、沿 ``z`` 单轴、沿 ``xy`` 双轴对同一个
4-layer Cu(111) slab 施加 ``[-3%, +3%]`` 范围的对称应变。每个方向得到
7 个 strained cell，再用一组 Cu-like 立方弹性常数合成对应的应变能曲线。

.. note::

   能量数据是合成的。模型把 slab 当作立方弹性固体，按 Voigt 简化应变能
   密度 ``u = 1/2 * C_IJ * eta_I * eta_J`` 计算每个方向的能量响应，并
   加上 ``σ = 0.5 meV/atom`` 的高斯噪声。模型故意忽略 Poisson 耦合，让
   三条曲线斜率明显不同；它不能替代真实 VASP 弛豫计算。

公式与单位
----------

立方晶系的二阶弹性应变能密度：

.. math::

   u(\eta) = \frac{1}{2} C_{11} (\eta_1^2 + \eta_2^2 + \eta_3^2)
   + C_{12} (\eta_1 \eta_2 + \eta_2 \eta_3 + \eta_3 \eta_1)
   + \frac{1}{2} C_{44} (\eta_4^2 + \eta_5^2 + \eta_6^2).

本教程只激活法向应变分量：

* 沿 ``x`` 单轴：``eta_1 = strain``，其余为 0，``u = 1/2 * C11 * strain^2``；
* 沿 ``z`` 单轴：``eta_3 = strain``，其余为 0，``u = 1/2 * C11 * strain^2``
  （立方对称下 ``C33 = C11``）；
* 沿 ``xy`` 双轴：``eta_1 = eta_2 = strain``，其余为 0，
  ``u = (C11 + C12) * strain^2``。

每原子能量 ``E = u * V0/N``，其中 ``V0`` 是 reference slab 体积，``N`` 是
原子数。``C11 = 168 GPa``、``C12 = 121 GPa`` 为 Cu 文献值（仅用于教学）。

输入与 provenance
-----------------

无外部文件依赖。Reference slab 由
:func:`mymetal.build.film.stretch.generate_film` 在运行时构建，应变系列由
:func:`stretch_list_along_direction_to_cell` 生成，能量由
:func:`docs.examples.biaxial_stretch.synthesize_energy` 合成，种子固定。

运行完整脚本
------------

.. code-block:: console

   $ python docs/examples/biaxial_stretch.py --output docs/_build/example-biaxial

.. literalinclude:: ../../examples/biaxial_stretch.py
   :language: python
   :linenos:

预期输出
--------

.. code-block:: text

   direction  strains(%)  E-E0 (meV/atom)
   x          [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]   20.003    8.984    2.799   -0.605    2.405    9.233   20.045
   z          [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]   20.814    9.384    2.625   -0.324    2.857    8.783   20.263
   xy         [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]   69.356   31.221    7.425   -0.512    7.076   30.904   69.915
   wrote: .../docs/_build/example-biaxial/biaxial_stretch.png

.. figure:: /_static/images/generated/biaxial_stretch.png
   :width: 960px
   :alt: Side views of Cu(111) slabs under -3% strain along x/z/xy, and the three strain-energy curves

   上排：三个方向在最大压缩 (``strain = -3%``) 下的侧视图（面内 3 倍重复）。
   下排：三个方向的应变能曲线。``xy`` 双轴的斜率约为 ``x``/``z`` 单轴的
   3.5 倍，因为 ``C11 + C12`` 是 ``C11/2`` 的约 3.5 倍（Cu 的
   ``C12/C11 ≈ 0.72``）。

结果含义
--------

* ``x`` 和 ``z`` 单轴曲线斜率接近（都由 ``C11`` 主导），但 ``z`` 方向
  因为真空层存在，物理上 slab 沿 z 拉伸几乎不会改变原子间斥力，真实
  VASP 计算中能量响应会远小于本教程的合成值。本教程的 ``z`` 曲线只用于
  展示 ``stretch_direction='z'`` 的几何效果。
* ``xy`` 双轴曲线在 ``strain = ±3%`` 处的能量约为 ``70 meV/atom``，是
  单轴的 3.5 倍。这是因为双轴同时激活 ``C11`` 和 ``C12``：两个方向都
  拉伸时，横向耦合项 ``C12 * eta1 * eta2`` 加倍贡献正能量。
* ``strain = 0`` 处的能量不为精确 0，因为有 ``σ = 0.5 meV/atom`` 噪声。
  断言允许 ``±1 meV`` 的偏差。

参数说明
--------

.. list-table:: :func:`stretch_list_along_direction_to_cell` 关键参数
   :header-rows: 1
   :widths: 28 14 58

   * - 参数
     - 默认
     - 含义
   * - ``atoms``
     - ``None``
     - 待拉伸的 ``ase.Atoms``。本教程传入 4-layer Cu(111) slab。
   * - ``stretch_factor_list``
     - ``[0.997, ..., 1.003]``
     - 拉伸因子列表。本教程设为 ``[0.97, 0.98, ..., 1.03]``。
   * - ``stretch_direction_list``
     - 12 个方向
     - 与 ``stretch_factor_list`` 一一对应的方向字符串。本教程每次只传
       一个方向，让函数内部把它广播到所有因子。
   * - ``stretch_3_direction_lists``
     - ``None``
     - 三轴独立拉伸的因子列表，形如
       ``[[fx, fy, fz], ...]``。提供时覆盖 ``stretch_factor_list``。
   * - ``my_scale_atoms``
     - ``True``
     - 是否同时缩放原子位置（保持 fractional 坐标不变）。

.. list-table:: 支持的方向字符串（节选）
   :header-rows: 1
   :widths: 30 70

   * - 方向
     - 拉伸的 cell 矢量
   * - ``"x"`` / ``"y"`` / ``"z"``
     - 单个笛卡尔方向
   * - ``"xy"`` / ``"xz"`` / ``"yz"``
     - 两个笛卡尔方向同时同比例
   * - ``"xyz"``
     - 三向同性拉伸（体积应变）
   * - ``"1"`` / ``"2"`` / ``"3"``
     - 单个晶格矢量 ``a1`` / ``a2`` / ``a3`` （与 ``x``/``y``/``z``
       在非正交 cell 下不同）

验证方法
--------

* 每个方向返回的 ``Atoms`` 数量必须等于 ``len(STRETCH_FACTORS)``；
* ``strain = 0`` 处的能量必须小于两端能量（加 ``1 meV`` 容差），保证
  抛物线开口向上；
* ``xy`` 方向在最大应变处的能量必须大于 ``x`` 方向，因为
  ``C11 + C12 > C11``。

常见错误
--------

``ValueError: \`atoms\` must be an instance of ASE Atoms.``
   ``atoms`` 未传入或传入 None。``stretch_list_along_direction_to_cell``
   不会自动构造参考结构。

``ValueError: \`stretch_direction\` must be one of ...``
   方向字符串不在支持的枚举里。注意 ``"XY"`` 大写不行，必须用小写
   ``"xy"``；``"1"``/``"2"``/``"3"`` 是晶格矢量下标，不是笛卡尔方向。

拉伸后 slab 沿 z 方向变薄但真空也变小
   ``stretch_along_direction_to_cell`` 缩放整个 cell，包括真空层。如果想
   保持真空厚度恒定，需要在拉伸后重新调用 ``cell.center(vacuum=...)``
   调整。

下一步
------

* :doc:`eos_curve` 看如何从 (V, E) 数据拟合 ``B0`` 和 ``B0'``；
* :doc:`../getting_started/au111_slab` 回到单 slab 的构建与 ``POSCAR``
  round-trip；
* :doc:`../manual/vasp` 了解 stretch workflow 脚本
  ``vasp_utils/vasp_workflow_bulk/pei_vasp_run_stretch`` 的位置。

Related API
-----------

* :func:`mymetal.build.film.stretch.stretch_list_along_direction_to_cell`
* :func:`mymetal.build.film.stretch.stretch_along_direction_to_cell`
* :func:`mymetal.build.film.stretch.generate_film`
