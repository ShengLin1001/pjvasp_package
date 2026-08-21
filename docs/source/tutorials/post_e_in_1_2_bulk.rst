.. _tutorial-post-e-in-1-2-bulk:

双轴变形 E_in_1/2 体相能量景观
================================

:Audience: 想理解双轴 a1/a2 扫描的能量景观和平衡参数提取的用户
:Time: 8 分钟
:Requires: Python、numpy、pjvasp_package
:Runs VASP: No

目标与边界
----------

演示 :func:`mymetal.post.E_in_1_2_bulk.post_E_in_1_2_bulk` 的分析方法。
使用合成数据展示 2D 能量景观 (E vs a1, a2)、等高线图和 profile 提取。

.. note::

   本教程使用合成数据演示。真实 E_in_1/2 计算需要扫描 a1/a2 网格的
   VASP 计算，函数本身需要 ``y_dir/`` 目录结构和 post data 文件。

背景
--------

E_in_1/2 (in-plane 1/2) 变形是六方晶系的双轴变形模式，独立扫描面内
晶格常数 ``a1`` 和 ``a2``，记录能量：

.. math::

   E(a_1, a_2) = f(a_1, a_2)

平衡参数由 2D 能量极小给出：

.. math::

   (a_1^0, a_2^0) = \\arg\\min_{a_1, a_2} E(a_1, a_2)

合成数据演示
------------

.. figure:: /_static/images/generated/post_E_in_1_2_bulk_demo.png
   :alt: E_in_1/2 bulk 2D energy landscape contour and profile plots

   左：2D 等高线图 E(a1, a2)，星标标记平衡位置。右：沿 a1 和 a2 方向
   的能量 profile，从平衡点截取。

结果含义
--------

* 能量景观呈碗形，平衡点在 ``a1 = a2 = a_eq`` 附近。
* 沿 ``a1 = a2`` 对角线方向曲率较大，反映双轴变形刚度。
* 沿单轴方向曲率较小，反映单轴变形刚度。
* ``c/a`` 比随 ``a1/a2`` 变化而变化，反映 Poisson 耦合。

参数说明
--------

.. list-table:: ``post_E_in_1_2_bulk`` 参数
   :header-rows: 1
   :widths: 25 15 60

   * - 参数
     - 默认
     - 含义
   * - ``jobn``
     - ``None``
     - 作业名列表 (a1-XX-a2-XX 格式)
   * - ``Etot``
     - ``None``
     - 总能量列表 (eV)
   * - ``atoms_ref``
     - ``None``
     - 参考结构
   * - ``latoms``
     - ``None``
     - 原子结构列表
   * - ``save_fig_path``
     - ``'p_post_E_in_1_2_bulk.pdf'``
     - 等高线图路径
   * - ``save_fig_path2``
     - ``'p_post_E_in_1_2_bulk2.pdf'``
     - profile 图路径

相关 API
--------

* :func:`mymetal.post.E_in_1_2_bulk.post_E_in_1_2_bulk`
* :func:`mymetal.post.E_in_1_2_bulk.my_write_E_in_1_2_bulk`
* :func:`mymetal.post.E_in_1_2_bulk.my_read_E_in_1_2_bulk`
