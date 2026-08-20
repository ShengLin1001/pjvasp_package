.. _tutorial-plot-gallery-plotting:

plotting 模块 —— pymatgen 风格绘图
==================================

:模块: :mod:`mymetal.universal.plot.plotting`
:函数数: 10

模块概述
--------

``plotting`` 模块源自 ``pymatgen.util.plotting``，提供出版级绘图工具：
:func:`pretty_plot` 系列函数、周期表热图、van Arkel 三角图，以及一批
Figure/Axes 工厂助手。

示例图
------

.. figure:: /_static/images/generated/plot_gallery_pretty.png
   :width: 960px
   :alt: pretty_plot and pretty_polyfit_plot demonstration

   左：``pretty_plot`` 绘制散点；右：``pretty_polyfit_plot`` 绘制数据 +
   一次多项式拟合趋势线。

.. figure:: /_static/images/generated/plot_gallery_periodic.png
   :width: 960px
   :alt: periodic_table_heatmap demonstration

   8 种纯元素体弹模量的周期表热图（matplotlib 后端）。

.. figure:: /_static/images/generated/plot_gallery_arkel.png
   :width: 640px
   :alt: van_arkel_triangle demonstration

   5 种二元化合物的 van Arkel-Ketelaar 三角图。

公开函数
--------

pretty_plot
~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.plotting.pretty_plot
   :no-index:

用途
   获取出版级 Axes，自动设置字号、黄金比例画布、palettable 色循环。

pretty_plot_two_axis
~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.plotting.pretty_plot_two_axis
   :no-index:

用途
   双 y 轴出版级图（左右不同量纲）。

pretty_polyfit_plot
~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.plotting.pretty_polyfit_plot
   :no-index:

用途
   绘制数据点 + 多项式拟合趋势线。

periodic_table_heatmap
~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.plotting.periodic_table_heatmap
   :no-index:

.. note::

   ``pymatviz`` 默认为 ``True``（返回 plotly figure，无法 ``savefig`` 到
   PNG）。需要 matplotlib 后端时传 ``pymatviz=False``。

format_formula
~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.plotting.format_formula
   :no-index:

用途
   将化学式字符串转为 LaTeX 格式，用于图例标注。

van_arkel_triangle
~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.plotting.van_arkel_triangle
   :no-index:

.. note::

   该函数 **不创建自己的 figure**，而是在 ``plt.gca()`` 上绘制。调用前需
   先 ``fig, ax = plt.subplots(figsize=(8, 7))``。输入是元素符号对的列
   表（如 ``[["Na", "Cl"], ["Ga", "As"]]``），不是化学式字符串。

get_ax_fig / get_ax3d_fig / get_axarray_fig_plt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.plotting.get_ax_fig
   :no-index:

.. autofunction:: mymetal.universal.plot.plotting.get_ax3d_fig
   :no-index:

.. autofunction:: mymetal.universal.plot.plotting.get_axarray_fig_plt
   :no-index:

用途
   工厂助手：接受可选 Axes 参数，为 None 时自动创建 figure + Axes。

add_fig_kwargs
~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.plotting.add_fig_kwargs
   :no-index:

用途
   装饰器：为返回 matplotlib figure 的函数添加 ``**fig_kwargs`` 透传。

相关 API
--------

* :mod:`mymetal.universal.plot.plotting`
* :doc:`../tutorials/periodic_table_and_arkel` （专题教程）
* :doc:`plot_gallery` （总览）
