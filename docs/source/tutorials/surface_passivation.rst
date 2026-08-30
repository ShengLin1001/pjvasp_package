.. _tut-surface-passivation:

Surface passivation and dangling-bond detection
================================================

目标
----

在 slab 表面添加吸附原子（如 H）钝化悬挂键，并使用悬挂键检测工具判断
哪些表面原子缺少近邻。本教程不需要 VASP、POTCAR 或 SLURM。

前置条件
--------

* 已安装 ``mymetal`` 和 ``ase``；
* 理解 slab、表面方向和 periodic boundary conditions。

物理问题
--------

切割 slab 时，表面原子的配位数低于 bulk。这些 dangling bonds（悬挂键）
会在 DFT 计算中引入人为的表面态。常用对策是在表面添加 H 或其他吸附原子
使表面原子恢复 bulk 配位。

核心代码
--------

``add_hydro_atoms`` —— 手动添加吸附原子
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

最直接的方法是按表面方向（通常是 z）筛选表面原子，然后在每个表面原子上方
或下方放置吸附原子：

.. code-block:: python

   from mymetal.build.film.stretch import generate_film
   from mymetal.build.film.hydroxyl import add_hydro_atoms

   slab = generate_film(
       symbols="Au", structure="fcc", num_layers=6,
       my_vacuum=15.0, slice_plane=(1, 1, 1), a_fcc=4.08,
   )
   z_max = slab.get_positions()[:, 2].max()
   z_min = slab.get_positions()[:, 2].min()

   # Top surface: H 1.8 Å above each Au
   passivated = add_hydro_atoms(
       slab.copy(), add_symbol="H", added_symbol="Au",
       surf_range=[z_max - 2.0, z_max + 5.0],
       shift_distance=[0.0, 0.0, 1.8],
       surf_direction=2,
   )

   # Bottom surface: H 1.8 Å below
   passivated = add_hydro_atoms(
       passivated, add_symbol="H", added_symbol="Au",
       surf_range=[z_min - 5.0, z_min + 2.0],
       shift_distance=[0.0, 0.0, -1.8],
       surf_direction=2,
   )

关键参数
~~~~~~~~

``add_symbol``
   要添加的原子符号（如 ``"H"``）。

``added_symbol``
   被钝化的表面原子符号（如 ``"Au"``、``"O"``、``"Si"``）。

``surf_range``
   沿 ``surf_direction`` 的坐标范围，只有在此范围内的
   ``added_symbol`` 原子才会被添加吸附原子。

``shift_distance``
   吸附原子相对于表面原子的位移向量（Å）。z 分量决定吸附高度。

``surf_direction``
   表面法线轴索引（0=x, 1=y, 2=z）。

``reverse_every``
   交替反转 shift 方向，适用于需要在两种位置间交替放置吸附原子的场景。

Convention
~~~~~~~~~~

* ``add_hydro_atoms`` 使用 Cartesian 坐标，``surf_range`` 和
  ``shift_distance`` 单位均为 Å；
* 吸附原子通过 ASE ``add_adsorbate`` 放置，PBC 设为 ``True``；
* 函数返回深拷贝，不修改输入 ``atoms``；
* 结果按 atomic number 排序并 ``wrap()``。

悬挂键检测
----------

``passivate_surface_custom`` 通过比较 bulk 和 slab 的近邻来自动检测
悬挂键。它需要 bulk 和 slab 的原子位置在相同坐标系下对齐：

.. code-block:: python

   from ase.build import bulk
   from mymetal.build.film.hydroxyl import passivate_surface_custom

   si_bulk = bulk("Si", "diamond", a=5.43, cubic=True)
   # slab 必须与 bulk 位置对齐 — 不要 center/vacuum-shift 后再传入

   passivated = passivate_surface_custom(
       bulk=si_bulk, slab=si_slab,
       adsorbates={"Si": "H"},
       cutoff=2.5,
       weights={"Si": 1.0},
       top_or_bottom=["top", "bottom"],
   )

.. warning::

   ``passivate_surface_custom`` 要求 slab 原子的 Cartesian 位置与 bulk
   中对应原子匹配（容差 1e-5 Å）。如果 slab 经过 ``center(vacuum=...)``
   或坐标平移，匹配会失败并打印 "Cannot find matching bulk atom"。
   正确做法：通过扩展 cell 添加真空，而不是平移原子。

底层工具
~~~~~~~~

* ``get_neighbors`` —— 给定原子索引、位置数组、cell 和 PBC，返回 cutoff
  内的近邻索引、位移向量和距离。
* ``find_matching_atom_in_bulk`` —— 在 bulk 位置数组中查找与 slab 原子
  位置匹配的原子，返回索引或 -1。
* ``find_unsaturated_neighbors`` —— 比较 bulk 和 slab 的近邻位移向量，
  返回 slab 中缺失的近邻方向（即悬挂键方向）。

完整可运行脚本
--------------

.. literalinclude:: ../../examples/surface_passivation.py
   :language: python

运行：

.. code-block:: shell

   python docs/examples/surface_passivation.py --output docs/_build/example-passivation

输出检查
--------

* Au(111) slab 公式为 ``Au6``；
* 顶面钝化后为 ``HAu6`` — 7 atoms；
* 双面钝化后为 ``H2Au6`` — 8 atoms；
* Si 原子 0 有 4 个近邻（diamond 配位），距离均为 2.351 Å；
* 悬挂键检测正确识别出 1 个缺失近邻。

常见错误
--------

* **Cannot find matching bulk atom**：slab 坐标与 bulk 不对齐。检查是否
  在传入前做了 ``center(vacuum=...)``。
* **TypeError: atom_index must be non-negative**：``find_matching_atom_in_bulk``
  返回 -1 但代码未检查。已在 ``passivate_surface_custom`` 中修复
  （``if idx_bulk >= 0``）。
* **吸附原子数量为 0**：``surf_range`` 没有覆盖到表面原子。检查
  ``z_max`` 和 ``z_min`` 的实际值。

相关 API
--------

* :func:`mymetal.build.film.hydroxyl.add_hydro_atoms`
* :func:`mymetal.build.film.hydroxyl.passivate_surface_custom`
* :func:`mymetal.build.film.hydroxyl.get_neighbors`
* :func:`mymetal.build.film.hydroxyl.find_matching_atom_in_bulk`
* :func:`mymetal.build.film.hydroxyl.find_unsaturated_neighbors`

下一步
------

* :doc:`atom_manipulation` —— 更通用的原子操作（FixAtoms、筛选、移动）；
* :doc:`../getting_started/au111_slab` —— Au(111) slab 构建；
* :doc:`surface_energy` —— 钝化后如何计算表面能。
