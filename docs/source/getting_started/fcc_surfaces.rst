.. _getting-started-fcc-surfaces:

构建并对比 FCC Cu(100)/(110)/(111) 表面
==========================================

:Audience: 想用 ``mymetal`` 一次构建多个低指数 FCC 表面的用户
:Time: 10 分钟
:Requires: Python 3.10、pjvasp_package、ASE、Matplotlib
:Runs VASP: No

目标
----

掌握 :func:`mymetal.build.film.stretch.generate_film` 在不同 Miller 指数
下的行为：用一个函数同时构建 Cu(100)、Cu(110) 和 Cu(111) 三个 12-layer
slab，对比面内 cell、c 方向长度和层数诊断，并写出可被 VASP 直接使用的
``POSCAR``。

最终结果
--------

.. figure:: /_static/images/generated/fcc_surfaces.png
   :width: 960px
   :alt: Side and top views of Cu(100), Cu(110) and Cu(111) 12-layer slabs

   上排是三个 slab 沿面内方向 3 倍重复后的侧视图，下排是 2x2 重复后的
   俯视图。三个 slab 共享同一 ``a_fcc = 3.61 Å``，但 (100)/(110)/(111)
   的面内排列和层间距完全不同。图片由
   ``docs/scripts/generate_structure_images.py`` 调用教程脚本确定性生成。

前置条件
--------

完成 :doc:`../user_guide/install` 后，确认：

.. code-block:: console

   $ python -c "import ase, matplotlib, mymetal; print('imports: ok')"
   imports: ok

.. note::

   本教程不运行 VASP，不调用 ``sbatch``，也不需要 ``POTCAR``。

输入和起始目录
--------------

从仓库根目录运行。输入只有一个可下载 Python 脚本：

.. code-block:: text

   pjvasp_package/
   └── docs/
       └── examples/
           └── fcc_surfaces.py

运行完整示例
------------

:download:`下载 fcc_surfaces.py <../../examples/fcc_surfaces.py>`，或在仓库
根目录直接运行：

.. code-block:: console

   $ python docs/examples/fcc_surfaces.py --output docs/_build/example-fcc-surfaces

页面与用户下载复用同一份实现：

.. literalinclude:: ../../examples/fcc_surfaces.py
   :language: python
   :linenos:

预期终端输出
------------

``generate_film`` 会为每个 Miller 指数打印层数诊断，随后脚本给出确定性
摘要表：

.. code-block:: text

   label     formula  atoms  a_in_plane_A  b_in_plane_A  c_total_A  pbc
   Cu(100)   Cu12     12     2.5527        2.5527        49.8550    [True, True, True]
   Cu(110)   Cu12     12     2.5527        3.6100        44.0396    [True, True, True]
   Cu(111)   Cu12     12     2.5527        2.5527        52.9266    [True, True, True]
   wrote: .../docs/_build/example-fcc-surfaces/fcc_surfaces.png

输出目录为：

.. code-block:: text

   docs/_build/example-fcc-surfaces/
   ├── POSCAR_Cu_100
   ├── POSCAR_Cu_110
   ├── POSCAR_Cu_111
   └── fcc_surfaces.png

如何解释结果
------------

* 三个 slab 的 ``formula`` 都为 ``Cu12``，``atoms`` 都为 12，验证了
  ``num_layers=12`` 的语义：它指的是原子层数，与面内排列无关。
* ``a_in_plane_A`` 都等于 ``2.5527`` Å。这是 ``a_fcc / sqrt(2)``，即 FCC
  最短原子间距，对应 (111)/(100) 面内 primitive cell 的边长。
* ``b_in_plane_A`` 在 (110) 处变成 ``3.6100`` Å（= ``a_fcc``），因为 (110)
  面的 rectangular primitive cell 沿 ``b`` 方向需要完整的 cubic 边长。
* ``c_total_A`` 等于 slab 物理厚度加上两侧各 15 Å 真空。
  (111) 最厚（层间距最密），(110) 最薄（层间距最疏）。
* ``pbc = [True, True, True]`` 是 ``generate_film`` 的默认行为：z 方向
  仍周期重复，真空用于隔开相邻周期像。

参数说明
--------

.. list-table:: :func:`generate_film` 在本教程使用的关键参数
   :header-rows: 1
   :widths: 22 12 66

   * - 参数
     - 默认
     - 含义
   * - ``symbols``
     - ``None``
     - 元素符号，例如 ``"Cu"``。仅单元素体系可用内部 bulk 构造路径。
   * - ``structure``
     - ``None``
     - 晶体类型，``"fcc"`` 或 ``"hcp"``。本教程使用 ``"fcc"``。
   * - ``num_layers``
     - ``None``
     - 期望的原子层数。若提供，优先于 ``replic_z``。
   * - ``replic_z``
     - ``None``
     - 沿表面法向的重复数。``num_layers`` 未给时使用。
   * - ``my_vacuum``
     - ``7.5``
     - slab 上下各加的真空厚度，单位 Å。本教程设为 ``15.0``。
   * - ``slice_plane``
     - ``(0,0,1)``
     - Miller 指数 ``(h, k, l)``。本教程依次设为 ``(1,0,0)``、
       ``(1,1,0)`` 和 ``(1,1,1)``。
   * - ``a_fcc``
     - ``2.95*sqrt(2)``
     - FCC 惯用胞边长，单位 Å。本教程统一设为 ``3.61`` （Cu 文献值）。
   * - ``my_periodic``
     - ``True``
     - 是否保持 z 方向周期。Slab 计算通常保持 ``True``，由真空隔开镜像。
   * - ``if_find_prim``
     - ``True``
     - 是否在切面后调用 ``my_find_prim`` 化简为 primitive cell。

验证方法
--------

脚本在写出图片前执行以下断言：

* 三个 slab 的 ``formula`` 都为 ``Cu12`` 且 ``atoms == 12``；
* 三个 slab 的 ``pbc`` 都为 ``[True, True, True]``；
* 三个 ``POSCAR`` 都能被 :func:`mymetal.io.vasp.my_read_vasp` 读回，且
  cell/positions 在 ``1e-10`` 容差内与原始 ``Atoms`` 一致；
* PNG 文件存在且包含足够的非白像素。

常见错误
--------

``ValueError: ... is an invalid structure``
   ``structure`` 只接受 ``"fcc"`` 和 ``"hcp"``。BCC 当前需要外部 ASE
   ``bulk()`` 构造后通过 ``bulk_atoms`` 参数传入。

``num_layers`` 与层数诊断不匹配
   ``generate_film`` 会先用 ``my_find_num_per_slab`` 计算 primitive slab
   的层数，再把 ``num_layers`` 除以它得到 ``num_rep_z``。若除不尽会向下
   取整。例如 Cu(100) 的 primitive slab 含 2 层，``num_layers=12`` 给出
   ``num_rep_z=6``，最终 12 个原子。

看不到真空
   俯视图沿 z 方向观察，不会显示真空；查看侧视图（图片上排），或比较
   ``c_total_A`` 与 slab 物理厚度。

下一步
------

* :doc:`hcp_surfaces` 学习 HCP 体系的多表面构建；
* :doc:`../tutorials/eos_curve` 看如何用 EOS 拟合得到体弹模量；
* :doc:`../tutorials/biaxial_stretch` 看如何对 slab 施加可控应变。

Related API
-----------

* :func:`mymetal.build.film.stretch.generate_film`
* :func:`mymetal.io.vasp.my_read_vasp`
* :func:`mymetal.io.vasp.my_write_vasp`
