.. _tutorial-plot-gallery:

绘图函数总览
============

:Audience: 需要使用 ``mymetal.universal.plot`` 绘制出版级图片的用户
:Time: 20 分钟
:Requires: Python、matplotlib、pjvasp_package
:Runs VASP: No

目标
----

本节系统介绍 :mod:`mymetal.universal.plot` 下全部 11 个模块的公开绘图函数，
覆盖通用样式助手、画布创建器、工作流后处理出图、能量分解、原子结构信息、
n2p2 训练诊断、pymatgen 风格绘图、OVITO 渲染、PowerPoint 导出以及 DOS 绘图。

每个子页面包含：

* 模块概述（用途、设计理念）；
* 每个公开函数的签名、参数说明、返回值、最小示例代码；
* 一张用合成数据生成的示例图（VASP-free）；
* 相关 API 交叉链接。

设计理念
--------

``mymetal.universal.plot`` 的目标是\"出版级、跨图一致、开箱即用\"。所有默认
字号、线宽、画布尺寸、刻度方向、图例样式都在入口函数里定死，用户只需要
调用 :func:`mymetal.universal.plot.plot.my_plot` 或
:func:`mymetal.universal.plot.general.general_set_all_rcParams` 即可获得统一
风格，不要从裸 ``plt.subplots()`` 手调样式。

两个标准入口二选一：

1. :func:`~mymetal.universal.plot.plot.my_plot` —— 直接拿到已配好样式的
   ``fig, axes``，最常用。
2. :func:`~mymetal.universal.plot.general.general_set_all_rcParams` —— 只改
   全局 rcParams，自己 ``plt.subplots()``；返回一个闭包用于修饰图例。

导出统一用 ``fig.savefig(path, bbox_inches='tight')``，不要叠
``tight_layout()``。每个 ``ax.legend()`` 后都要
:func:`~mymetal.universal.plot.general.general_modify_legend`。

运行完整脚本
------------

.. code-block:: console

   $ python docs/examples/plot_gallery_demo.py --output docs/_build/example-plot-gallery

该脚本为每个模块生成一张合成数据示例 PNG，全部 VASP-free、确定性、非空白。

模块导航
--------

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - 模块
     - 函数数
     - 简介
   * - :doc:`plot_gallery_general`
     - 18
     - 通用样式 / 标注 / 采色 / 自省助手
   * - :doc:`plot_gallery_plot`
     - 5
     - 画布创建器（my_plot / my_plot_colorbar / 断轴）
   * - :doc:`plot_gallery_workflow`
     - 13
     - 工作流后处理出图（收敛 / NEB / 拉伸 / HOEC）
   * - :doc:`plot_gallery_energy`
     - 1
     - 能量分量分解绘图
   * - :doc:`plot_gallery_atominfo`
     - 3
     - 层间距 / z 位置 / 径向分布函数
   * - :doc:`plot_gallery_n2p2`
     - 8
     - n2p2 训练诊断（学习曲线 / RMSE / epoch 监控）
   * - :doc:`plot_gallery_plotting`
     - 10
     - pymatgen 风格绘图（pretty_plot / 周期表热图 / van Arkel 三角图）
   * - :doc:`plot_gallery_render`
     - 1
     - OVITO 管线渲染
   * - :doc:`plot_gallery_ppt`
     - 1
     - PowerPoint 幻灯片批量导出图片
   * - :doc:`plot_gallery_oldplotdos`
     - 5
     - DOS / IDOS / 元素投影 DOS 绘图

相关 API
--------

* :mod:`mymetal.universal.plot.general`
* :mod:`mymetal.universal.plot.plot`
* :mod:`mymetal.universal.plot.workflow`
* :mod:`mymetal.universal.plot.energy`
* :mod:`mymetal.universal.plot.atominfo`
* :mod:`mymetal.universal.plot.n2p2`
* :mod:`mymetal.universal.plot.plotting`
* :mod:`mymetal.universal.plot.render`
* :mod:`mymetal.universal.plot.ppt`
* :mod:`mymetal.universal.plot.oldplotdos`
