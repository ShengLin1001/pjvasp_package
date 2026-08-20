.. _tutorial-plot-gallery-oldplotdos:

oldplotdos 模块 —— DOS 绘图
===========================

:模块: :mod:`mymetal.universal.plot.oldplotdos`
:函数数: 5

模块概述
--------

``oldplotdos`` 模块用于绘制 VASP 态密度（DOS）和积分态密度（IDOS），
支持总 DOS、元素投影 s/p/d DOS、自旋极化、堆叠填充、Fermi 能级归零等。
:func:`my_plot_complete_dos` 和 :func:`my_plot_idos` 从 ``vasprun.xml``
文件读取数据；:func:`my_plot_horizontal_vertical` 和
:func:`my_plot_orientation` 是底层画布和方向助手。

.. note::

   ``my_plot_complete_dos`` / ``my_plot_idos`` /
   ``my_plot_element_spd_dos`` 需要真实的 ``vasprun.xml`` 文件，超出本
   VASP-free 画廊的范围。下图用合成高斯 DOS 数据演示绘图风格。

示例图
------

.. figure:: /_static/images/generated/plot_gallery_dos.png
   :width: 960px
   :alt: DOS plotting style demonstration

   用合成高斯数据绘制的总 DOS（左）和 s/p 投影 DOS（右），Fermi 能级
   归零，虚线标 Fermi 能级。

公开函数
--------

my_plot_horizontal_vertical
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.oldplotdos.my_plot_horizontal_vertical
   :no-index:

用途
   创建水平或垂直方向的 DOS 画布，设置坐标范围和标签。

my_plot_orientation
~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.oldplotdos.my_plot_orientation
   :no-index:

用途
   按指定方向（水平 / 垂直）绘制 DOS 曲线，支持填充、能带 / 电子数
   参考线。

my_plot_complete_dos
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.oldplotdos.my_plot_complete_dos
   :no-index:

用途
   从 ``vasprun.xml`` 读取并绘制完整 DOS，支持总 DOS + 元素 s/p/d 投影、
   堆叠、自旋极化、高斯展宽。

.. note::

   需要真实 ``vasprun.xml`` 文件。支持多元素 DOS。

my_plot_idos
~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.oldplotdos.my_plot_idos
   :no-index:

用途
   从 ``vasprun.xml`` 读取并绘制积分态密度（IDOS）。

my_plot_element_spd_dos
~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.oldplotdos.my_plot_element_spd_dos
   :no-index:

用途
   从 ``vasprun.xml`` 读取并绘制各元素的 s/p/d 投影 DOS。

相关 API
--------

* :mod:`mymetal.universal.plot.oldplotdos`
* :doc:`plot_gallery_plot` （画布创建器）
* :doc:`plot_gallery` （总览）
