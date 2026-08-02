.. _tutorial-periodic-table-and-arkel:

周期表热图与 van Arkel 三角图
=============================

:Audience: 需要可视化元素性质分布或化学键特征的用户
:Time: 5 分钟
:Requires: Python、matplotlib、pjvasp_package
:Runs VASP: No

目标
----

演示 :mod:`mymetal.universal.plot.plotting` 中两个特色可视化函数：

1. :func:`mymetal.universal.plot.plotting.periodic_table_heatmap` —
   在周期表上按元素值着色；
2. :func:`mymetal.universal.plot.plotting.van_arkel_triangle` —
   把二元材料按电负性绘制在 van Arkel-Ketelaar 三角图上。

背景
----

**周期表热图** 适合展示元素性质的周期性趋势（如体弹模量、结合能、
电负性）。

**van Arkel 三角图** 用电负性差（纵轴 Δχ）和电负性均值（横轴 χ̄）定位
二元化合物，三个顶角分别代表离子键（大 Δχ）、金属键（低 χ̄）和共价键
（高 χ̄，小 Δχ）。

运行完整脚本
------------

.. code-block:: console

   $ python docs/examples/periodic_table_and_arkel.py --output docs/_build/example-periodic-arkel

.. literalinclude:: ../../examples/periodic_table_and_arkel.py
   :language: python
   :linenos:

结果图
------

.. figure:: /_static/images/generated/periodic_table_heatmap.png
   :width: 960px
   :alt: Periodic table heatmap of bulk moduli for selected pure elements

   15 种纯元素的体弹模量（GPa）周期表热图。过渡金属（W、Fe、Cu）体弹模量
   较高，碱金属（Na、K）较低。

.. figure:: /_static/images/generated/van_arkel_triangle.png
   :width: 720px
   :alt: Van Arkel-Ketelaar triangle for NaCl, MgO, Al2O3, SiO2, GaAs, ZnS, CuBr, InP

   8 种二元化合物的 van Arkel 三角图。NaCl、MgO 位于离子区（大 Δχ），
   GaAs、InP 位于共价区（小 Δχ），CuBr 位于中间。

结果含义
--------

* **周期表热图**：W ``310 GPa`` 最高，K ``3.5 GPa`` 最低，反映过渡金属
  d 键强度。
* **van Arkel 三角图**：NaCl ``Δχ=2.23`` 最离子，GaAs ``Δχ=0.37`` 最共价。

验证方法
--------

* 输入 dict 非空，材料列表非空；
* 两张图片均非空白。

相关 API
--------

* :func:`mymetal.universal.plot.plotting.periodic_table_heatmap`
* :func:`mymetal.universal.plot.plotting.van_arkel_triangle`
* :func:`mymetal.universal.plot.plotting.format_formula`
