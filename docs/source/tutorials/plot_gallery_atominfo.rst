.. _tutorial-plot-gallery-atominfo:

atominfo 模块 —— 原子结构信息
==============================

:模块: :mod:`mymetal.universal.plot.atominfo`
:函数数: 3

模块概述
--------

``atominfo`` 模块从 :class:`ase.Atoms` 对象提取结构信息并绘图：层间距、
z 位置、径向分布函数（RDF）。所有函数内部调用
:func:`~mymetal.universal.plot.plot.my_plot` 获取统一样式。

示例图
------

.. figure:: /_static/images/generated/plot_gallery_interlayer.png
   :width: 640px
   :alt: my_plot_interlayer_distance demonstration

   合成 Cu 薄板的层间距相对偏差图，以中间层为参考。

.. figure:: /_static/images/generated/plot_gallery_zpositions.png
   :width: 640px
   :alt: my_plot_zpositions demonstration

   原子 z 位置排序图，右上角标注材料厚度。

.. figure:: /_static/images/generated/plot_gallery_rdf.png
   :width: 640px
   :alt: my_plot_rdf demonstration

   合成 Cu 结构的径向分布函数直方图。

公开函数
--------

my_plot_interlayer_distance
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.atominfo.my_plot_interlayer_distance
   :no-index:

用途
   计算并绘制层间距。以中间层为参考（``d_ref``），绘制
   ``d_i / d_ref - 1``。可选保存图片和 txt 数据。

my_plot_zpositions
~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.atominfo.my_plot_zpositions
   :no-index:

用途
   绘制原子 z 位置的排序图，并标注材料厚度（Å）。

my_plot_rdf
~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.atominfo.my_plot_rdf
   :no-index:

用途
   计算并绘制径向分布函数（RDF）直方图。依赖
   :func:`mymetal.universal.atom.neighbor.get_neighbor_distances`。

相关 API
--------

* :mod:`mymetal.universal.plot.atominfo`
* :mod:`mymetal.universal.atom.neighbor`
* :doc:`plot_gallery` （总览）
