.. _tutorial-cubic-cell-and-stretch:

正交 cell 查找与单轴拉伸
=========================

:Audience: 需要把 primitive 薄膜 cell 转成正交 cell 并施加应变扫描的用户
:Time: 8 分钟
:Requires: Python、ASE、pjvasp_package
:Runs VASP: No

目标与边界
----------

演示两个常用薄膜操作：

1. 用 :func:`mymetal.build.film.findcubic.find_cubic` 把 HCP Mg primitive
   薄膜 cell（六角面内）转成 xy 面正交 cell，面积守恒；
2. 用
   :func:`mymetal.build.film.stretch.stretch_along_direction_to_cell`
   对正交 cell 沿 x 方向施加 ``[0.99, 1.00, 1.01]`` 拉伸因子。

不运行 VASP。最终输出一个 2×2 对比图。

为什么要找正交 cell
-------------------

许多后处理脚本（能量-应变拟合、LAMMPS data 写入）假设面内 cell 是正交的。
但 :func:`mymetal.build.film.stretch.generate_film` 切出的 primitive
HCP/FCC 薄膜面内通常是六角（120° γ）或非正交的。``find_cubic`` 把它转成
``[a', b', c, 90, 90, 90]`` 形式，同时保持 xy 面积不变：

.. math::

   a' = \sqrt{\frac{A}{\sqrt{3}}}, \quad b' = \sqrt{3}\, a',

其中 ``A`` 是 :func:`mymetal.build.film.extrfilm.cal_area` 给出的 xy 投影
面积（Å²）。

运行完整脚本
------------

.. code-block:: console

   $ python docs/examples/cubic_cell_and_stretch.py --output docs/_build/example-cubic-stretch

.. literalinclude:: ../../examples/cubic_cell_and_stretch.py
   :language: python
   :linenos:

预期输出
--------

.. code-block:: text

   label                       formula  atoms  a_A      b_A      c_A      alpha   beta    gamma   area_A2
   primitive (hex)             Mg4      4      3.2100   3.2100   27.8629  90.00   90.00   120.00  8.9236
   orthorhombic (cubic)        Mg4      4      2.2698   3.9314   27.8629  90.00   90.00   90.00   8.9236
   stretch x * 0.99            Mg4      4      2.2471   3.9314   27.8629  90.00   90.00   90.00   8.8344
   stretch x * 1.00            Mg4      4      2.2698   3.9314   27.8629  90.00   90.00   90.00   8.9236
   stretch x * 1.01            Mg4      4      2.2925   3.9314   27.8629  90.00   90.00   90.00   9.0128
   wrote: .../docs/_build/example-cubic-stretch/cubic_cell_and_stretch.png
   OK: all assertions passed.

结果图
------

.. figure:: /_static/images/generated/cubic_cell_and_stretch.png
   :width: 960px
   :alt: Hexagonal primitive cell, orthorhombic cell, and three uniaxially stretched cells

   2×2 对比图。左上：primitive 六角 cell（γ=120°）；右上：正交化后的
   cell（γ=90°，面积守恒）；左下/中下/右下：沿 x 方向拉伸 0.99/1.00/1.01
   后的正交 cell。

结果含义
--------

* **面积守恒**：primitive → orthorhombic 的 xy 面积均为
  ``8.9236 Å²``，验证 ``find_cubic`` 不改变面内面积。
* **正交化**：``a' = 2.2698 Å``、``b' = 3.9314 Å``，满足
  ``b' = sqrt(3) * a'``，γ 从 120° 变为 90°。
* **拉伸**：``stretch x * 0.99`` 把 ``a'`` 缩小 1%（``2.2471 Å``），
  ``b'`` 和 ``c`` 不变；``1.01`` 放大 1%（``2.2925 Å``）。

验证方法
--------

脚本执行以下断言：

* primitive 和 orthorhombic 的 :func:`cal_area` 结果一致（面积守恒）；
* orthorhombic 的 γ = 90°；
* 拉伸后 ``a'`` 按因子缩放，``b'``、``c`` 不变；
* 图片非空白。

常见错误
--------

* ``find_cubic`` 输入 cell 不是 ``[a,b,c,any,any,90]`` 形式（即 c 轴
  不垂直于 xy 面）时，结果会不正确。先用 ``generate_film`` 切出标准
  薄膜再调用。
* ``stretch_along_direction_to_cell`` 的 ``stretch_direction_list``
  只接受 ``'x'``/``'y'``/``'z'``。

下一步
------

* 把拉伸后的 cell 写成 POSCAR 批量提交 VASP 静态计算，扫描能量-应变
  曲线，见 :doc:`../manual/workflows`。
* 双轴拉伸见 :doc:`biaxial_stretch`。

相关 API
--------

* :func:`mymetal.build.film.findcubic.find_cubic`
* :func:`mymetal.build.film.findcubic.find_optimal_cell_shape`
* :func:`mymetal.build.film.stretch.stretch_along_direction_to_cell`
* :func:`mymetal.build.film.stretch.generate_film`
* :func:`mymetal.build.film.extrfilm.cal_area`
