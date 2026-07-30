.. _tutorial-schmid-factor:

计算 FCC 滑移系的 Schmid 因子
=============================

:Audience: 想用 ``mymetal`` 计算单晶 FCC 在任意加载方向下的 Schmid 因子，
   并找出主导滑移系的用户
:Time: 10 分钟
:Requires: Python、pandas、numpy、pjvasp_package
:Runs VASP: No

目标与边界
----------

演示 :func:`mymetal.calculate.material_science.schmid.cal_fcc_schmid_factors`
的典型用法：

1. 给定加载方向 ``[u, v, w]``，自动枚举 FCC ``{111}<110>`` 全部 12 个
   滑移系，计算每个系的 Schmid 因子 ``m = |cos(λ) cos(φ)|``。
2. 对 10 个标准拉伸方向做扫描，给出最大 Schmid 因子与激活系数量。
3. 在 (001) 立方极图上做密集扫描，可视化最大 Schmid 因子的分布。

.. note::

   Schmid 因子是纯几何量，不依赖材料常数、温度或 VASP 计算。本教程的
   所有数值都是确定性的。

公式与单位
----------

FCC 完美位错柏氏矢量为 ``a/2<110>``，滑移面为 ``{111}``。对加载方向
``L = [u, v, w]``、滑移面法向 ``n`` 和滑移方向 ``b``，Schmid 因子为：

.. math::

   m = \left| \cos\lambda \cdot \cos\phi \right|
     = \left| \frac{\mathbf{L} \cdot \mathbf{b}}
                       {|\mathbf{L}|\,|\mathbf{b}|}
              \cdot \frac{\mathbf{L} \cdot \mathbf{n}}
                       {|\mathbf{L}|\,|\mathbf{n}|}
       \right|.

只保留 ``b`` 落在 ``n`` 所定义平面内的组合（即 ``b · n = 0``），FCC 共
12 个滑移系。其中 ``[1, 1, 1]`` 加载方向对 ``(1, 1, 1)`` 面本身 ``m = 0``
（因为 ``cos(φ) = 1`` 且 ``cos(λ) = 0``），但对其他 ``{111}`` 面仍可能
激活。

输入与 provenance
-----------------

无外部文件依赖。FCC 滑移系硬编码在
:func:`cal_fcc_schmid_factors` 中；加载方向列表与极图扫描网格在
``docs/examples/schmid_factor.py`` 中定义，种子固定。

运行完整脚本
------------

.. code-block:: console

   $ python docs/examples/schmid_factor.py --output docs/_build/example-schmid

.. literalinclude:: ../../examples/schmid_factor.py
   :language: python
   :linenos:

预期输出
--------

.. code-block:: text

   reference direction: [1, 1, 6]
     total slip systems: 12
     active (m > 1e-6):   10
     max Schmid factor:  0.451222

   scan summary:
   direction       max_schmid  n_active
     [1 0 0]       0.408248    8
     [1 1 0]       0.408248    4
     [1 1 1]       0.272166    6
     [1 1 2]       0.408248    8
     [1 1 3]       0.445362    10
     [1 1 6]       0.451222    10
     [1 2 3]       0.466569    9
     [1 3 4]       0.471056    9
     [2 1 3]       0.466569    9
     [3 2 1]       0.466569    9
   wrote: .../docs/_build/example-schmid/schmid_factor.png

.. figure:: /_static/images/generated/schmid_factor.png
   :width: 960px
   :alt: Polar Schmid-factor map and per-slip-system bar chart for FCC [1, 1, 6] loading

   左图：(001) 立方极图上的最大 Schmid 因子分布。颜色越亮代表最大
   Schmid 因子越大；接近 ``[1, 1, 1]`` 角的区域 ``m`` 较小（接近
   ``0.27``），靠近 ``[1, 1, 2]``--``[1, 1, 6]`` 边界的区域 ``m`` 较大
   （接近 ``0.47``）。右图：``[1, 1, 6]`` 加载方向下，全部 12 个 FCC
   滑移系的 Schmid 因子水平条形图；红色条为激活系（``m > 1e-6``），
   灰色条为休眠系。

结果含义
--------

* ``[1, 1, 6]`` 加载方向激活 10 个滑移系，最大 Schmid 因子为
  ``0.4512``；这是 FCC 单晶拉伸实验中典型的多滑移方向。
* ``[1, 0, 0]`` 与 ``[1, 1, 2]`` 的最大 Schmid 因子相同（都是
  ``0.4082``），因为它们对称等价。
* ``[1, 1, 1]`` 加载方向最大 Schmid 因子只有 ``0.2722``，对应于非
  ``(1, 1, 1)`` 面的滑移。这是 FCC 单晶 ``[1, 1, 1]`` 轴向被视为
  "硬" 方向的几何原因。
* ``[1, 3, 4]`` 给出本扫描的最大值 ``0.4711``，但仍小于理论上限
  ``0.5``。

参数说明
--------

.. list-table:: :func:`cal_fcc_schmid_factors` 关键参数
   :header-rows: 1
   :widths: 22 14 64

   * - 参数
     - 默认
     - 含义
   * - ``normal_orientation``
     - ``[1, 1, 6]``
     - 加载方向 ``[u, v, w]``。整数三元组；不需要归一化。
   * - 返回值
     - -
     - ``pandas.DataFrame``，列包括 ``slip_plane``、``dislocation``、
       ``force_direction``、``schmid_factor``、``cos_lambda``、``cos_phi``。

.. list-table:: 本教程使用的扫描方向
   :header-rows: 1
   :widths: 30 70

   * - 方向
     - 含义
   * - ``[1, 0, 0]``、``[1, 1, 0]``、``[1, 1, 1]``
     - 立方胞三个高对称轴。
   * - ``[1, 1, 2]``、``[1, 1, 3]``、``[1, 1, 6]``
     - ``[1, 1, 1]``--``[1, 1, 0]`` 边界上的标准方向。
   * - ``[1, 2, 3]``、``[1, 3, 4]`` 等
     - 立方极图内部的常见拉伸方向。

验证方法
--------

* ``[1, 1, 6]`` 必须激活至少一个滑移系；
* 所有 Schmid 因子必须落在 ``[0, 0.5]``；
* ``[1, 1, 1]`` 加载方向对 ``(1, 1, 1)`` 滑移面本身的 ``m`` 必须为 0
  （加载方向平行于面法向），但其他 ``{111}`` 面可被激活；
* ``[1, 1, 6]`` 的最大 Schmid 因子必须严格大于 ``[1, 1, 1]``。

常见错误
--------

``KeyError: 'schmid_factor'``
   旧版 ``mymetal`` 返回 ``dict`` 而非 ``DataFrame``。运行
   ``python -m pip install -e .`` 重装包后再试。

返回的 Schmid 因子全为 0
   检查加载方向 ``[0, 0, 0]`` 是否被误传；方向必须非零。

``[1, 1, 1]`` 给出非零 Schmid 因子
   这是预期行为：``[1, 1, 1]`` 不激活 ``(1, 1, 1)`` 面，但激活
   ``(-1, 1, 1)``、``(1, -1, 1)``、``(1, 1, -1)`` 三个面。

下一步
------

* :doc:`strain_deformation` 看如何把应变张量与滑移方向对应；
* :doc:`../manual/vasp` 了解 ``vasp_utils`` 中的拉伸/剪切 workflow；
* :doc:`../getting_started/bulk_structures` 回到 FCC bulk 构建。

Related API
-----------

* :func:`mymetal.calculate.material_science.schmid.cal_fcc_schmid_factors`
* :class:`pandas.DataFrame`
