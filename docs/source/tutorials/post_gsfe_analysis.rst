.. _tutorial-post-gsfe:

广义层错能 (GSFE) 分析
======================

:Audience: 想理解 GSFE 后处理如何提取层错能曲线的用户
:Time: 8 分钟
:Requires: Python、numpy、pjvasp_package
:Runs VASP: No

目标与边界
----------

演示 :func:`mymetal.post.gsfe.post_gsfe` 的分析方法。使用 LAMMPS UNNEP 势
计算的 FCC(111) GSFE 真实结果，展示 gamma-surface 曲线、稳定/不稳定层错能
的识别、以及结构约束检查。

.. note::

   本教程的数据来自 LAMMPS 计算（``pj-test-properties-gold/y_gsfe/fcc/FCC_111/``），
   使用 UNNEP 机器学习势。``post_gsfe`` 函数支持 VASP 和 LAMMPS 两种输入格式。

背景公式
--------

广义层错能 (GSFE) 定义为：

.. math::

   \\gamma = \\frac{\\Delta E}{A_{sf}} \\quad [\\text{mJ/m}^2]

其中 ``Delta E`` 是层错构型与参考构型的能量差，``A_sf`` 是滑移面面积。

稳定层错 (SF) 对应 gamma 曲线的局部极小值，不稳定层错 (USF) 对应局部极大值。

真实数据结果
------------

.. list-table:: FCC(111) GSFE 结果 (LAMMPS UNNEP)
   :header-rows: 1
   :widths: 30 25 45

   * - 属性
     - 值
     - 说明
   * - 滑移面面积
     - 7.488 Ang^2
     - Asf
   * - a11
     - 2.941 Ang
     - 晶格常数 a1
   * - a22
     - 2.547 Ang
     - 晶格常数 a2
   * - E0_bulk
     - -3.227 eV
     - 每原子体相能量
   * - USFE
     - 106.49 mJ/m^2
     - 不稳定层错能 (gamma 最大值)

.. figure:: /_static/images/generated/post_gsfe_fcc_111.png
   :alt: FCC(111) GSFE gamma vs displacement curve with SF/USF markers

   FCC(111) GSFE gamma 曲线。21 个数据点覆盖 ``1/3 <11-2>`` 滑移路径。
   USF (不稳定层错) 对应 gamma 最大值 106.49 mJ/m^2。

结果含义
--------

* ``USFE = 106.49 mJ/m^2`` 是 FCC(111) 的关键塑性参数，决定了
  ``1/3 <11-2>`` 分位错的核心宽度和临界分解剪切应力 (CRSS)。
* gamma 曲线在位移 ``0.275 a11`` 处达到最大值，对应
  ``1/3 <11-2>`` 路径的不稳定位置。
* gamma 值随后下降到约 48 mJ/m^2 的局部极小 (SF)，对应稳定层错位置。

结构约束检查
------------

``post_gsfe`` 通过 :func:`check_constraints` 验证：

* 所有构型化学组分一致；
* 面内晶格常数不发生弛豫（仅允许面外弛豫）；
* 原子位移仅沿滑移方向（面内无额外弛豫）；
* 滑移方向沿一条直线。

相关 API
--------

* :func:`mymetal.post.gsfe.post_gsfe`
* :func:`mymetal.post.gsfe.check_constraints`
* :func:`mymetal.post.gsfe.find_sf_usf`
* :func:`mymetal.post.gsfe.write_output`
* :func:`mymetal.post.gsfe.read_output`
