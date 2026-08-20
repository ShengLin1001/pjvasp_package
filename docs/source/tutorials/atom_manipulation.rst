.. _tutorial-atom-manipulation:

移动、固定原子并写出 Selective dynamics
=======================================

:Audience: 想用 ``mymetal`` 给 slab 加约束、按子集平移原子并写出
   ``Selective dynamics`` POSCAR 的用户
:Time: 10 分钟
:Requires: Python、ASE、numpy、pjvasp_package
:Runs VASP: No

目标与边界
----------

演示三个最常用的 universal atom helpers：

1. :func:`mymetal.universal.atom.fixatom.fixatoms` -- 用布尔 ``mask``
   或整数 ``indices`` 给 ``Atoms`` 添加 :class:`ase.constraints.FixAtoms`
   约束。
2. :func:`mymetal.universal.atom.moveatom.move_atoms` -- 在分数或笛卡尔
   坐标下，按元素种类和 Cartesian 位置范围筛选原子，然后整体平移。
3. :func:`mymetal.io.vasp.my_write_vasp` -- 把约束写回 ``POSCAR``，使
   下一次 VASP 弛豫只放开未固定的原子。

.. note::

   本教程使用 6-layer Au(111) slab。不运行 VASP，不调用 ``sbatch``，
   也不需要 ``POTCAR``。所有坐标变换都是确定性的。

公式与单位
----------

``move_atoms`` 的 ``translate_matrix`` 在 ``if_scale_position=True``
时按 fractional 坐标加法：

.. math::

   \mathbf{s}_i' = \mathbf{s}_i + \mathbf{t}_{\mathrm{frac}},
   \quad \mathbf{r}_i' = \mathbf{s}_i' \cdot \mathbf{A}.

在 ``if_scale_position=False`` 时按 Cartesian 坐标加法：

.. math::

   \mathbf{r}_i' = \mathbf{r}_i + \mathbf{t}_{\mathrm{cart}}.

筛选条件为：

.. math::

   i \in \mathcal{I} \iff
   (\mathrm{symbol}_i \in \mathrm{types}) \wedge
   (x_i^{\min} \le x_i \le x_i^{\max}) \wedge \ldots

``FixAtoms`` 约束在 ``POSCAR`` 中以 ``Selective dynamics`` 块表示，
每行末尾三个 ``T``/``F`` 标记对应 ``x``/``y``/``z`` 三个方向是否自由
（``T`` = free，``F`` = fixed）。

输入与 provenance
-----------------

无外部文件依赖。Reference slab 由
:func:`mymetal.build.film.stretch.generate_film` 在运行时构建；平移量
固定为 ``[0.25, 0, 0]`` fractional，对应 Au(111) 顶层从 atop 位置移到
bridge 位置。

运行完整脚本
------------

.. code-block:: console

   $ python docs/examples/atom_manipulation.py --output docs/_build/example-manipulation

.. literalinclude:: ../../examples/atom_manipulation.py
   :language: python
   :linenos:

预期输出
--------

.. code-block:: text

   slab formula: Au6
   slab atoms: 6
   frozen atoms: 3
   moving atoms: 3
   top layer shift (fractional): [0.25, 0.0, 0.0]
   POSCAR round-trip: ok
     formula: Au6
     atoms:   6
     constraint: FixAtoms
   wrote: .../docs/_build/example-manipulation/POSCAR
   wrote: .../docs/_build/example-manipulation/atom_manipulation.png

.. figure:: /_static/images/generated/atom_manipulation.png
   :alt: Reference Au(111) slab, slab with bottom half frozen, and slab with top layer shifted to bridge site

   左图：6-layer Au(111) 参考结构，没有任何约束。中图：底层 3 个原子
   被 ``FixAtoms`` 冻结（红色虚线框包围），顶层 3 个原子保持原位。右图：
   在中图基础上把顶层 3 个原子沿 ``x`` 方向平移 ``0.25`` fractional
   （= ``1/4 a``），从 atop 位置移到 bridge 位置；红色箭头标出平移方向。

结果含义
--------

* 6-layer slab 的底层 3 个原子被冻结，顶层 3 个原子保持可动；
* ``POSCAR`` 中第一行注释保留；``Selective dynamics`` 块正确写出；
* ``my_read_vasp`` 读回时识别 ``FixAtoms`` 约束并恢复为
  ``ase.constraints.FixAtoms`` 对象；
* 顶层原子的 fractional 坐标变化恰为 ``[0.25, 0, 0]``，底层原子坐标
  不变（容差 ``1e-10``）。

参数说明
--------

.. list-table:: :func:`fixatoms` 关键参数
   :header-rows: 1
   :widths: 22 14 64

   * - 参数
     - 默认
     - 含义
   * - ``atoms``
     - 必填
     - ``ase.Atoms`` 对象。
   * - ``mask``
     - ``None``
     - 布尔列表，长度等于 ``len(atoms)``。``True`` 表示该原子被冻结。
   * - ``indices``
     - ``None``
     - 整数列表，指定被冻结原子的下标。``mask`` 与 ``indices`` 二选一。

.. list-table:: :func:`move_atoms` 关键参数
   :header-rows: 1
   :widths: 22 14 64

   * - 参数
     - 默认
     - 含义
   * - ``atoms``
     - 必填
     - ``ase.Atoms`` 对象。
   * - ``translate_matrix``
     - ``[0.1, 0.1, 0.0]``
     - 平移向量。``if_scale_position=True`` 时为 fractional，否则为
       Cartesian Å。
   * - ``if_scale_position``
     - ``True``
     - 是否在 fractional 坐标下平移。
   * - ``atom_type``
     - ``None``
     - 元素符号列表，例如 ``["O", "H"]``。``None`` 表示所有元素。
   * - ``position_range``
     - ``((-inf, inf), (-inf, inf), (-inf, inf))``
     - Cartesian 位置范围 ``((xmin, xmax), (ymin, ymax), (zmin, zmax))``，
       单位 Å。只有同时满足三个方向的原子才会被平移。

.. list-table:: :func:`my_write_vasp` Selective dynamics 相关参数
   :header-rows: 1
   :widths: 22 14 64

   * - 参数
     - 默认
     - 含义
   * - ``ignore_constraints``
     - ``False``
     - ``True`` 时忽略 ``atoms.constraints`` 并写出无约束的 ``POSCAR``。
   * - ``wrap``
     - ``True``
     - 写出前是否把原子 wrap 回 cell 内。

验证方法
--------

* slab 原子数为 6（偶数），可对半切分；
* ``mask`` 中 ``True`` 的数量恰好等于 ``len(slab) // 2``；
* ``POSCAR`` round-trip 后 ``atoms.constraints`` 至少包含一个对象；
* 底层原子的 Cartesian 坐标在 ``1e-10`` 容差内不变；
* 顶层原子的 fractional 坐标变化恰为 ``[0.25, 0, 0]``。

常见错误
--------

``RuntimeError: VASP requires that the direction of FixedPlane ...``
   ``POSCAR`` 写出时遇到 ``FixedPlane`` 或 ``FixedLine`` 约束，且方向
   不平行于任何 cell 轴。改用 ``FixAtoms`` 或 ``FixScaled`` 约束。

平移后原子跑出 cell
   ``move_atoms`` 默认会调用 ``atoms.wrap()`` 把原子 wrap 回 cell。
   如果不希望 wrap，需要在调用后手动 ``atoms.set_scaled_positions(...)``。

底层原子也被移动
   检查 ``position_range`` 是否覆盖到底层。本教程用
   ``z_max - 0.1 <= z <= z_max + 0.1`` 严格限制在顶层。

``POSCAR`` 没出现 ``Selective dynamics``
   ``atoms.constraints`` 为空。检查 ``fixatoms`` 是否返回新对象（不是
   原地修改）；旧对象的 ``atoms.constraints`` 不会被更新。

下一步
------

* :doc:`../getting_started/au111_slab` 回到 slab 构建的最小例子；
* :doc:`biaxial_stretch` 看如何对整个 cell 做应变系列；
* :doc:`../manual/vasp` 了解 ``vasp_utils`` 中的 ``pei_vasp_univ_*``
  约束自动生成脚本。

Related API
-----------

* :func:`mymetal.universal.atom.fixatom.fixatoms`
* :func:`mymetal.universal.atom.moveatom.move_atoms`
* :func:`mymetal.io.vasp.my_write_vasp`
* :class:`ase.constraints.FixAtoms`
