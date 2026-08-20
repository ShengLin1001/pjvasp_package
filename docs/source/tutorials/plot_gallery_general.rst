.. _tutorial-plot-gallery-general:

general 模块 —— 通用样式与标注助手
==================================

:模块: :mod:`mymetal.universal.plot.general`
:函数数: 18

模块概述
--------

``general`` 模块是整个 ``mymetal.universal.plot`` 包的样式基础。它定义了
全局 rcParams、字体、刻度、边距、图例、线条颜色等默认值，并提供一批标注
助手（色带、编号圆圈、箭头、自动避让文字、辅助线）。所有入口函数
（:func:`~mymetal.universal.plot.plot.my_plot` 等）都依赖这里的
:func:`general_font` / :func:`general_axes` / :func:`general_subplots_adjust`。

设计理念：把\"跨图一致\"的默认值集中在一处，用户不要从裸
``plt.subplots()`` 手调样式，而应通过 :func:`general_set_all_rcParams` 或
:func:`~mymetal.universal.plot.plot.my_plot` 获取统一风格。

示例图
------

.. figure:: /_static/images/generated/plot_gallery_general.png
   :width: 960px
   :alt: general module helpers demonstration

   四个子图分别演示 ``add_color_band`` + ``generate_gradient_colors``、
   ``add_circle_number`` + ``add_arrow``、``general_add_vlines_hlines`` +
   ``general_modify_line``、``general_adjust_text`` + ``general_margin_bin``。

公开函数
--------

general_set_all_rcParams
~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.general.general_set_all_rcParams
   :no-index:

用途
   一次性配置全部 matplotlib rcParams，返回一个 ``_general_modify_legend``
   闭包，画完图例后调用它即可统一图例样式。适合需要自己 ``plt.subplots()``
   的场景。

最小示例
   .. code-block:: python

      from mymetal.universal.plot.general import general_set_all_rcParams
      import matplotlib.pyplot as plt

      lg = general_set_all_rcParams(figure_subp=(1, 1))
      fig, ax = plt.subplots()
      ax.plot([0, 1], [0, 1], label="line")
      lg(ax.legend())

general_modify_legend
~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.general.general_modify_legend
   :no-index:

用途
   修饰图例边框：Square 样式、黑边白底、不透明。每个 ``ax.legend()`` 后
   都要调用。

add_color_band
~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.general.add_color_band
   :no-index:

用途
   在已有 Axes 上叠加半透明渐变色带，常用于标记相变区间或参考区域。

add_circle_number
~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.general.add_circle_number
   :no-index:

用途
   在数据坐标处绘制带编号的空心椭圆，用于在图上标注关键点序号。

add_arrow
~~~~~~~~~

.. autofunction:: mymetal.universal.plot.general.add_arrow
   :no-index:

用途
   在两点之间添加箭头（可带文字），支持 ``data`` 和 ``axes`` 坐标系。

generate_gradient_colors
~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.general.generate_gradient_colors
   :no-index:

用途
   生成两种颜色之间的渐变色列表，或从 colormap 采样；支持 reshape 为
   ``imshow`` 可用的形状。

general_modify_line
~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.general.general_modify_line
   :no-index:

用途
   统一修改 Axes 上所有线条的样式：移除虚线、按色表重新着色或渐变着色。

general_add_vlines_hlines
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.general.general_add_vlines_hlines
   :no-index:

用途
   批量添加垂直 / 水平辅助线（默认灰色虚线）。

general_adjust_text
~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.general.general_adjust_text
   :no-index:

用途
   自动调整文字位置以避免与标记点和坐标区边界重叠（adjustText 的封装）。

general_margin_bin
~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.general.general_margin_bin
   :no-index:

用途
   调整坐标轴边距和刻度分箱（MaxNLocator + prune）。

general_font
~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.general.general_font
   :no-index:

用途
   配置全局字体、网格、线宽、标记、图例等 rcParams（被 ``my_plot`` 内部
   调用）。

general_axes
~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.general.general_axes
   :no-index:

用途
   修饰单个 Axes 的刻度方向、标签、边距、顶右刻度开关。

general_subplots_adjust
~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.general.general_subplots_adjust
   :no-index:

用途
   按绝对英寸精算子图位置和间距（被 ``my_plot`` 内部调用）。

check_font_size / check_axes_size / check_all_rcParams
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.general.check_font_size
   :no-index:

.. autofunction:: mymetal.universal.plot.general.check_axes_size
   :no-index:

.. autofunction:: mymetal.universal.plot.general.check_all_rcParams
   :no-index:

用途
   自省助手：打印当前 Axes 的字号 / 尺寸 / 全部 rcParams，用于调试样式
   不一致问题。

get_ploted_figure
~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.general.get_ploted_figure
   :no-index:

用途
   返回当前 pyplot 状态中的 figure 对象。

get_points_on_markers_boundary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.general.get_points_on_markers_boundary
   :no-index:

用途
   生成围绕标记点边界的采样点，供 ``general_adjust_text`` 避让使用。

general_modify_band_plot
~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.general.general_modify_band_plot
   :no-index:

用途
   修饰带状图（band plot）的显示样式。

相关 API
--------

* :mod:`mymetal.universal.plot.general`
* :doc:`plot_gallery_plot` （画布创建器）
* :doc:`plot_gallery` （总览）
