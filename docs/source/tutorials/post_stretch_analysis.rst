.. _tutorial-post-stretch:

单轴拉伸后处理与平衡晶格常数提取
================================

:Audience: 想理解 stretch 后处理如何从能量-应变曲线提取平衡参数的用户
:Time: 8 分钟
:Requires: Python、numpy、pjvasp_package
:Runs VASP: No

目标与边界
----------

演示 :func:`mymetal.post.stretch.post_stretch` 和
:func:`mymetal.post.stretch.post_lammps_stretch` 背后的分析方法。使用
Au HCP (VASP) 和 FCC (LAMMPS UNNEP 势) 的真实计算结果，展示拉伸曲线
的二次多项式拟合和平衡晶格常数的提取。

.. note::

   本教程同时展示 VASP-DFT 和 LAMMPS-NNP 两套真实数据。VASP 数据来自
   zcm6 集群 A11-1 工作流，LAMMPS 数据来自 ``pj-test-properties-gold``。
   两者的拟合方法和后处理流程完全一致。

背景公式
--------

对单轴拉伸，能量-应变关系用二次多项式拟合：

.. math::

   E(f) = a\\,f^2 + b\\,f + c

其中 ``f`` 是拉伸因子。平衡拉伸因子和能量由极值给出：

.. math::

   f_0 = -\\frac{b}{2a}, \\quad E_0 = c - \\frac{b^2}{4a}

平衡晶格常数 ``a_0 = f_0 * a_ref``。

真实数据结果
------------

.. list-table:: 拉伸后处理结果对比
   :header-rows: 1
   :widths: 25 25 25 25

   * - 属性
     - VASP Au HCP
     - LAMMPS FCC (UNNEP)
     - 说明
   * - 拉伸类型
     - xy (双轴)
     - xyz (三轴)
     - 拉伸方向
   * - 原子数
     - 2
     - 4
     - 每个超胞原子数
   * - a_ref (Ang)
     - 2.8597
     - 4.1585
     - 参考晶格常数
   * - a_0 (Ang)
     - 2.8593
     - 4.1585
     - 平衡晶格常数
   * - E_0 (eV/atom)
     - -3.9179
     - -3.2269
     - 平衡能量

.. figure:: /_static/images/generated/post_stretch_comparison.png
   :alt: VASP Au HCP and LAMMPS FCC stretch energy-strain curves with quadratic fits

   上：VASP Au HCP 双轴拉伸能量曲线（17 个数据点 + 二次拟合）。
   下：LAMMPS FCC 三轴拉伸能量曲线（101 个数据点 + 二次拟合）。
   平衡晶格常数由拟合极值给出。

结果含义
--------

* VASP Au HCP 的平衡 ``a_0 = 2.8593 Ang``，与参考值偏差仅 0.01%，
  说明参考结构已接近平衡。
* LAMMPS FCC (UNNEP 势) 的 ``a_0 = 4.1585 Ang``，与 4.1585 Ang 参考值
  完全一致，说明 NNP 势在零应变附近表现良好。
* 两条曲线的开口方向一致（``a > 0``），符合稳定性要求。

参数说明
--------

.. list-table:: ``post_stretch`` / ``post_lammps_stretch`` 关键参数
   :header-rows: 1
   :widths: 25 15 60

   * - 参数
     - 默认
     - 含义
   * - ``dirsurf``
     - ``'y_stretch'``
     - 包含拉伸计算的目录
   * - ``refcontcar``
     - ``'./y_full_relax/CONTCAR'``
     - 参考结构文件
   * - ``save_fig_path``
     - ``'p_post_stretch.pdf'``
     - 图表保存路径
   * - ``save_txt_path``
     - ``'p_post_stretch.txt'``
     - 结果文件路径

相关 API
--------

* :func:`mymetal.post.stretch.post_stretch`
* :func:`mymetal.post.stretch.post_lammps_stretch`
* :func:`mymetal.post.stretch.get_stretch_type`
* :func:`mymetal.post.stretch.my_write_stretch`
* :func:`mymetal.post.stretch.my_read_stretch`
* :doc:`biaxial_stretch` （拉伸结构生成）
* :doc:`eos_curve` （状态方程拟合）
