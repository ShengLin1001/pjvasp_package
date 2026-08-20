.. _tutorial-plot-gallery-plot:

plot 模块 —— 画布创建器
========================

:模块: :mod:`mymetal.universal.plot.plot`
:函数数: 5

模块概述
--------

``plot`` 模块提供画布创建器，是用户最常调用的入口。:func:`my_plot` 封装了
:func:`~mymetal.universal.plot.general.general_font` +
:func:`~mymetal.universal.plot.general.general_subplots_adjust` +
:func:`~mymetal.universal.plot.general.general_axes`，返回已配好样式的
``(fig, axes)``，用户只需画数据、设标签、修饰图例、保存。

其他创建器：:func:`my_plot_brokenaxed` （断轴）、
:func:`my_plot_colorbar` （带 colorbar）、
:func:`my_plot_modify_ploted_figure` （改造已画好的图）。

示例图
------

.. figure:: /_static/images/generated/plot_gallery_plot.png
   :width: 960px
   :alt: my_plot two-panel demonstration

   用 ``my_plot(fig_subp=[1, 2])`` 创建的双子图，分别绘制 sin 和 cos 曲线，
   每个子图画完图例后调用 ``general_modify_legend``。

.. figure:: /_static/images/generated/plot_gallery_colorbar.png
   :width: 640px
   :alt: my_plot_colorbar demonstration

   用 ``my_plot_colorbar`` 创建的主坐标区 + 右侧 colorbar 坐标区，展示
   随机矩阵的 imshow。

公开函数
--------

my_plot
~~~~~~~

.. autofunction:: mymetal.universal.plot.plot.my_plot
   :no-index:

用途
   标准画布创建器。返回 ``(fig, axes)``，``axes`` 在单子图时为单个 Axes，
   多子图时为 ndarray。``fig_sharex=True`` 时建议用
   ``fig.subplots_adjust(hspace=...)`` 压缩行间距。

最小示例
   .. code-block:: python

      from mymetal.universal.plot.plot import my_plot
      from mymetal.universal.plot.general import general_modify_legend

      fig, ax = my_plot()
      ax.plot([0, 1, 2], [0, 1, 4], marker='o', label='data')
      ax.set_xlabel('x'); ax.set_ylabel('y')
      general_modify_legend(ax.legend())
      fig.savefig('out.png', bbox_inches='tight')

my_plot_brokenaxed
~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.plot.my_plot_brokenaxed
   :no-index:

用途
   创建断轴画布（依赖 ``brokenaxes`` 包），用于 y 轴有断裂区间的数据。

my_plot_modify_ploted_figure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.plot.my_plot_modify_ploted_figure
   :no-index:

用途
   改造已画好的图：设置标签、刻度、边距、线条颜色、图例、辅助线、保存
   等。常用于 DOS 等需要先画数据再统一修饰的场景。

my_plot_colorbar
~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.plot.my_plot_colorbar
   :no-index:

用途
   创建带右侧 colorbar 坐标区的画布。返回 ``(fig, (ax_main, ax_cbar))``。

.. note::

   返回值的第二个元素是 **元组** ``(ax_main, ax_cbar)``，不是单个 Axes。
   用 ``ax_main, ax_cbar = axes[0], axes[1]`` 解包。

general_modify_ploted_figure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.plot.general_modify_ploted_figure
   :no-index:

用途
   ``my_plot_modify_ploted_figure`` 的废弃别名，保留向后兼容。

相关 API
--------

* :mod:`mymetal.universal.plot.plot`
* :doc:`plot_gallery_general` （样式助手）
* :doc:`plot_gallery` （总览）
