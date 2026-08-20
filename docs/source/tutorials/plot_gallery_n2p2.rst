.. _tutorial-plot-gallery-n2p2:

n2p2 模块 —— 训练诊断绘图
==========================

:模块: :mod:`mymetal.universal.plot.n2p2`
:函数数: 8

模块概述
--------

``n2p2`` 模块用于 n2p2 神经网络势训练后处理诊断：学习曲线、DFT-vs-NNP
散点、按 tag 的 RMSE 柱状图、epoch 监控（拉伸量、弹性常数、GSFE、RMSE）、
nnp-LAMMPS 一致性检查。所有函数用
:func:`~mymetal.universal.plot.plot.my_plot` 创建画布，
:func:`~mymetal.universal.plot.general.general_modify_legend` 修饰图例。

示例图
------

.. figure:: /_static/images/generated/plot_gallery_learning.png
   :width: 960px
   :alt: my_plot_learning_curve demonstration

   训练 RMSE 随 epoch 的对数 y 曲线（左：能量，右：力）。

.. figure:: /_static/images/generated/plot_gallery_compare.png
   :width: 960px
   :alt: my_plot_compare demonstration

   DFT-vs-NNP 散点图（上：全范围，下：裁剪百分位），按 tag 着色。

.. figure:: /_static/images/generated/plot_gallery_rmse_tag.png
   :width: 720px
   :alt: my_plot_rmse_by_tag demonstration

   按 tag 的能量 / 力 RMSE 柱状图（对数 y 轴），TOTAL 为参考虚线。

.. figure:: /_static/images/generated/plot_gallery_epoch_rmse.png
   :width: 720px
   :alt: my_plot_epoch_rmse demonstration

   训练集能量 / 力 RMSE 随 epoch 的变化。

公开函数
--------

my_plot_learning_curve
~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.n2p2.my_plot_learning_curve
   :no-index:

用途
   绘制能量 / 力训练 RMSE 随 epoch 的对数 y 曲线。

my_plot_compare
~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.n2p2.my_plot_compare
   :no-index:

用途
   DFT-vs-NNP 散点图，按 tag 着色。上排全范围，下排裁剪百分位避免极端值
   压缩主数据。

my_plot_rmse_by_tag
~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.n2p2.my_plot_rmse_by_tag
   :no-index:

用途
   按 tag 的能量 / 力 RMSE 柱状图（对数 y 轴），柱顶标注数值，TOTAL 行
   为参考虚线。支持可选的不对称误差棒（``df_lo`` / ``df_hi``）。

my_plot_epoch_stretch
~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.n2p2.my_plot_epoch_stretch
   :no-index:

用途
   平衡拉伸量随 epoch 的变化（3 行 × 相列）。

my_plot_epoch_cij
~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.n2p2.my_plot_epoch_cij
   :no-index:

用途
   弹性常数 Cij 随 epoch 的变化（5 行 × 相列）。

my_plot_epoch_gsfe
~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.n2p2.my_plot_epoch_gsfe
   :no-index:

用途
   层错能随 epoch 的变化（2 行 × 滑移系列）。

my_plot_epoch_rmse
~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.n2p2.my_plot_epoch_rmse
   :no-index:

用途
   训练集能量 / 力 RMSE 随 epoch 的变化（2 行 × 1 列）。

my_plot_check_interface
~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.n2p2.my_plot_check_interface
   :no-index:

用途
   nnp 预测 vs LAMMPS hdnnp 一致性判定（2 行 × 1 列）。

相关 API
--------

* :mod:`mymetal.universal.plot.n2p2`
* :doc:`plot_gallery` （总览）
