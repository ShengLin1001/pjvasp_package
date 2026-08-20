.. _tutorial-kpoints-sampling:

生成与对比 VASP k 点网格
========================

:Audience: 已经会用 ``ase.build.bulk`` 构建 bulk/slab，但不知道如何选
   ``KPOINTS`` 的用户
:Time: 10 分钟
:Requires: Python、ASE、numpy、pjvasp_package
:Runs VASP: No

目标与边界
----------

演示 :mod:`mymetal.calculate.calqm.kpoints` 提供的三类常用 k 点工具：

1. :func:`get_kpoints_by_size` 同时生成 Monkhorst-Pack 与 Gamma-centered
   k 点，并对比它们在 (kx, ky) 平面上的位置。
2. :func:`get_size_by_distance` 用 ``RK`` 产品（旧 VASP ``round`` 法则）
   或 ``KSPACING``（新 VASP ``ceil`` 法则）为 slab 自动选 k 点。
3. :func:`cal_reciprocal_matrix` 与 :func:`cal_reciprocal_matrix2`
   两种等价的倒易点阵计算（叉乘 vs 矩阵求逆），用于验证实现一致性。

.. note::

   本教程使用一个 4-layer Cu(111) slab（含 20 Å 真空）作为参考结构。
   所有 k 点都由函数在运行时生成，不读取外部 ``KPOINTS`` 文件，也不
   调用 VASP。

公式与单位
----------

Monkhorst-Pack 网格在分数坐标下的位置为：

.. math::

   \mathbf{k}_{i,j,l} = \left(
     \frac{2i - n_x - 1}{2 n_x},\,
     \frac{2j - n_y - 1}{2 n_y},\,
     \frac{2l - n_z - 1}{2 n_z}
   \right), \quad
   i \in [1, n_x].

当某个方向的网格数 ``n`` 为偶数时，Gamma-centered 网格会沿该方向整体
平移 ``1/(2n)``，使 ``Γ`` 点本身进入网格。

VASP 的 ``KSPACING`` 与 ``RK`` 产品关系为：

.. math::

   R_k = \frac{2\pi}{\mathrm{KSPACING}}, \quad
   n_i = \mathrm{round}\!\left(R_k \cdot |\mathbf{b}_i|\right)
   \quad (\text{old VASP}),

.. math::

   n_i = \left\lceil \frac{2\pi \cdot |\mathbf{b}_i|}{\mathrm{KSPACING}}
   \right\rceil \quad (\text{new VASP}).

其中 ``|b_i|`` 是倒易点阵矢量的长度，单位 ``2π/Å``。

输入与 provenance
-----------------

无外部文件依赖。Reference slab 由
:func:`mymetal.build.film.stretch.generate_film` 在运行时构建，
``RK`` 列表固定为 ``[20, 40, 60, 80, 100, 120]``。

运行完整脚本
------------

.. code-block:: console

   $ python docs/examples/kpoints_sampling.py --output docs/_build/example-kpoints

.. literalinclude:: ../../examples/kpoints_sampling.py
   :language: python
   :linenos:

预期输出
--------

.. code-block:: text

   slab formula: Cu4
   slab atoms: 4
   slab cell lengths (A): [2.5527, 2.5527, 46.2527]
   MP k-points (first 5): [[-0.4167, -0.4167, 0.0], ...]
   Gamma k-points (first 5): [[-0.3333, -0.3333, 0.0], ...]
   MP count: 36, Gamma count: 36
   RK list: [20, 40, 60, 80, 100, 120]
   old VASP mesh per RK: [[9, 9, 1], [18, 18, 1], [27, 27, 1], [36, 36, 2], [45, 45, 2], [54, 54, 3]]
   new VASP mesh per RK: [[10, 10, 1], [19, 19, 1], [28, 28, 2], [37, 37, 2], [46, 46, 3], [55, 55, 3]]
   wrote: .../docs/_build/example-kpoints/kpoints_sampling.png

.. figure:: /_static/images/generated/kpoints_sampling.png
   :alt: Monkhorst-Pack vs Gamma-centered k-points, RK scan, and in-plane/out-of-plane ratio

   左图：6x6x1 Monkhorst-Pack（蓝圆）与 Gamma-centered（红叉）k 点在
   (kx, ky) 平面的位置。Gamma-centered 网格包含 Γ 点，Monkhorst-Pack
   不包含。中图：随 ``RK`` 增大，``n_x`` 与 ``n_z`` 的增长曲线；``n_z``
   受 20 Å 真空压制增长很慢。右图：``n_x / n_z`` 比值随 ``RK`` 的变化，
   反映 slab 计算中面内/法向 k 点密度的不对称。

结果含义
--------

* 6x6x1 网格共 36 个 k 点。Monkhorst-Pack 与 Gamma-centered 数量相同，
  但位置不同：前者最靠近 Γ 的点是 ``1/12``，后者直接包含 Γ。
* ``RK = 80`` 时旧 VASP 给出 ``[36, 36, 2]``，新 VASP 给出
  ``[37, 37, 2]``。差异来自 ``round`` 与 ``ceil`` 的边界处理。
* ``n_z`` 在 ``RK = 20-60`` 区间保持为 1，因为 slab 的 c 方向有 20 Å
  真空，倒易 ``|b_z|`` 很小。直到 ``RK = 80`` 才升级为 2。
* 倒易点阵的两种实现（叉乘 vs 矩阵求逆）在 ``1e-10`` 容差内一致。

参数说明
--------

.. list-table:: :func:`get_kpoints_by_size` 关键参数
   :header-rows: 1
   :widths: 22 14 64

   * - 参数
     - 默认
     - 含义
   * - ``size``
     - ``(1, 1, 1)``
     - ``(nx, ny, nz)`` 网格大小。本教程设为 ``(6, 6, 1)``。
   * - ``offset``
     - ``(0.5, 0.5, 0.5)``
     - 偏移量。当某方向 ``n`` 为偶数时，实际偏移为 ``offset[i] / n``。

.. list-table:: :func:`get_size_by_distance` 关键参数
   :header-rows: 1
   :widths: 22 14 64

   * - 参数
     - 默认
     - 含义
   * - ``atoms``
     - 必填
     - ``ase.Atoms``，函数会从中读取 ``cell.reciprocal()``。
   * - ``rk``
     - ``100``
     - ``R_k`` 产品，单位 ``2π/Å``。与 ``KSPACING = 2π/R_k`` 等价。
   * - ``kspacing``
     - ``None``
     - 直接指定 ``KSPACING`` (Å⁻¹)。若提供则覆盖 ``rk``。

验证方法
--------

* Monkhorst-Pack 与 Gamma-centered 网格的 k 点总数必须相等；
* 两种网格在偶数 ``size`` 下必须有差异（不能完全重合）；
* 倒易点阵的两种实现必须一致到 ``1e-10``；
* ``n_x`` 随 ``RK`` 单调不减；
* slab 的 ``n_z`` 受真空压制，在 ``RK <= 120`` 时 ``n_z <= 3``。

常见错误
--------

``ValueError: atoms must be ASE Atoms``
   ``get_size_by_distance`` 必须传 ``Atoms`` 对象，不能传文件名。若需要
   从 ``CONTCAR`` 读取，先 ``atoms = ase.io.read('CONTCAR')``。

``n_z`` 远大于预期
   检查 slab 是否有足够真空。若 c 方向没有真空，``|b_z|`` 会按 bulk
   比例缩放，``n_z`` 会与 ``n_x`` 接近。

Gamma-centered 与 Monkhorst-Pack 完全重合
   这只发生在所有方向网格数都为奇数时。偶数 ``size`` 才能看到差异。

下一步
------

* :doc:`../getting_started/bulk_structures` 回到 bulk 构建的最小例子；
* :doc:`../manual/vasp` 了解 ``vasp_utils`` 中 ``KPOINTS`` 自动生成脚本；
* :doc:`surface_energy` 看 k 点密度如何影响表面能收敛。

Related API
-----------

* :func:`mymetal.calculate.calqm.kpoints.get_kpoints_by_size`
* :func:`mymetal.calculate.calqm.kpoints.get_size_by_distance`
* :func:`mymetal.calculate.calqm.kpoints.cal_reciprocal_matrix`
* :func:`mymetal.calculate.calqm.kpoints.cal_reciprocal_matrix2`
* :class:`ase.dft.kpoints.monkhorst_pack`
