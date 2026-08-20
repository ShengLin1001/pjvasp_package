.. _tutorial-plot-gallery-workflow:

workflow 模块 —— 工作流后处理出图
==================================

:模块: :mod:`mymetal.universal.plot.workflow`
:函数数: 13

模块概述
--------

``workflow`` 模块是 VASP 工作流后处理出图的集合，涵盖收敛测试、弹性常数
拟合、GSFE、NEB、拉伸分析、HOEC、KPAR/NCORE 基准等。所有函数命名遵循
``my_plot_<workflow>`` 模式，返回 ``(fig, axes)``，用 ``if_save`` /
``savefile`` 控制落盘。

示例图
------

.. figure:: /_static/images/generated/plot_gallery_convergence.png
   :width: 720px
   :alt: my_plot_convergence demonstration

   ``my_plot_convergence`` 绘制能量截止收敛测试（能量差 vs 截止能）。

.. figure:: /_static/images/generated/plot_gallery_relax.png
   :width: 720px
   :alt: my_plot_relax_convergence demonstration

   ``my_plot_relax_convergence`` 绘制离子弛豫的能量和力收敛曲线（对数 y 轴）。

.. figure:: /_static/images/generated/plot_gallery_kpar_ncore.png
   :width: 720px
   :alt: my_plot_kpar_ncore demonstration

   ``my_plot_kpar_ncore`` 绘制 KPAR/NCORE 并行基准的耗时和能量偏差。

.. figure:: /_static/images/generated/plot_gallery_stretch.png
   :width: 720px
   :alt: my_plot_stretch demonstration

   ``my_plot_stretch`` 拟合二次能量 - 拉伸曲线并绘制 c/a 比。

公开函数
--------

my_plot_convergence
~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.workflow.my_plot_convergence
   :no-index:

用途
   绘制能量截止或 k 点网格的收敛测试。``if_difference=True`` 时以最后一个
   值为参考绘制能量差。``y`` 需传 numpy 数组（内部做 ``y - y[-1]``）。

my_plot_cij_energy
~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.workflow.my_plot_cij_energy
   :no-index:

用途
   绘制 Cij 能量 - 密度拟合曲线和应力响应。

my_plot_gsfe
~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.workflow.my_plot_gsfe
   :no-index:

用途
   绘制广义层错能（GSFE）、非弹性法向位移和剪应力。

my_plot_gsfe_displacement
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.workflow.my_plot_gsfe_displacement
   :no-index:

用途
   绘制每个 GSBE 镜像的原子分辨法向位移。

my_plot_relax_convergence
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.workflow.my_plot_relax_convergence
   :no-index:

用途
   绘制离子弛豫的能量和力收敛曲线，支持对数 y 轴和 EDIFFG 参考线。

.. note::

   该函数在 ``if_save=True`` 时用 ``fig.savefig(savefile)`` 保存，并可选
   ``if_close=True`` 关闭图。不需要再在外部调用 ``plt.savefig``。

my_plot_hoec_energy
~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.workflow.my_plot_hoec_energy
   :no-index:

用途
   绘制各变形模式的高阶弹性常数（HOEC）能量 - 应变数据和多项式拟合。

my_plot_hoec_convergence
~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.workflow.my_plot_hoec_convergence
   :no-index:

用途
   绘制每个求解的 HOEC 分量随拟合窗口的变化。

my_plot_kpar_ncore
~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.workflow.my_plot_kpar_ncore
   :no-index:

用途
   绘制 KPAR/NCORE 并行基准的耗时和能量偏差曲线。

.. note::

   ``dict_time`` 和 ``dict_delta_energy`` 的键是 ``(kpar, ncore)`` **元组**，
   不是嵌套 dict。

my_plot_neb
~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.workflow.my_plot_neb
   :no-index:

用途
   绘制 NEB 能量和力曲线。

my_plot_neb_full
~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.workflow.my_plot_neb_full
   :no-index:

用途
   绘制 NEB 全部能量曲线，支持帧切片。

my_plot_neb_xy
~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.workflow.my_plot_neb_xy
   :no-index:

用途
   绘制原子二维轨迹，支持逐原子 / 逐帧模式。

my_plot_stretch
~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.workflow.my_plot_stretch
   :no-index:

用途
   拟合二次能量 - 拉伸曲线，提取平衡值，绘制能量和 c/a 比。

.. note::

   ``rvectors_ref`` 不能为 None（函数内部做 ``rvectors_ref * extr_x``）。

my_plot_E_in_1_2_bulk
~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.workflow.my_plot_E_in_1_2_bulk
   :no-index:

用途
   生成二维等高线图和剖面图，用于双轴拉伸体能量分析。

相关 API
--------

* :mod:`mymetal.universal.plot.workflow`
* :doc:`plot_gallery_plot` （画布创建器）
* :doc:`plot_gallery` （总览）
