.. _tutorial-neighbor-distances:

计算近邻距离与径向分布函数
==========================

:Audience: 想用 ``mymetal`` 计算 bulk 结构的近邻距离、RDF 与累计配位数
   的用户
:Time: 10 分钟
:Requires: Python、ASE、numpy、pjvasp_package
:Runs VASP: No

目标与边界
----------

演示 :func:`mymetal.universal.atom.neighbor.get_neighbor_distances`：

1. 构建 FCC Cu、BCC Fe、HCP Mg 三个 2x2x2 supercell。
2. 用 ``cutoff = 8 Å`` 的截断半径收集每个结构的全部 pair distance
   （含周期镜像）。
3. 把距离分箱为 RDF 直方图，归一化后绘制为标准 ``g(r)`` 曲线。
4. 同时绘制累计配位数曲线（integrated RDF），直观显示前几个壳层的
   配位数。

.. note::

   本教程只使用 ``ase.build.bulk`` 与 ``mymetal`` 的 neighbor helper，
   不读取任何 ``CONTCAR`` 或 ``OUTCAR``。所有数值都是确定性的。

公式与单位
----------

径向分布函数 ``g(r)`` 定义为：

.. math::

   g(r) = \frac{1}{N \rho}\,
          \frac{1}{4\pi r^2 \Delta r}\,
          \left\langle N(r, r+\Delta r) \right\rangle,

其中 ``N`` 是 supercell 中原子数，``ρ = N/V`` 是数密度，``Δr`` 是 bin
宽度。本教程把 ``g(r)`` 进一步除以 ``max(g)`` 让三种结构在同一 y 轴
尺度上对比。

累计配位数（integrated RDF）：

.. math::

   N(r) = \int_0^r 4\pi r'^2 \rho\, g(r')\, \mathrm{d}r'
        \approx \sum_{r_i \le r} \mathrm{count}(r_i).

输入与 provenance
-----------------

无外部文件依赖。三个 bulk 结构在
:func:`docs.examples.neighbor_distances.build_structures` 中用
``ase.build.bulk`` 直接构建；``cutoff = 8 Å``、``bin_width = 0.05 Å``
均为教学值。

运行完整脚本
------------

.. code-block:: console

   $ python docs/examples/neighbor_distances.py --output docs/_build/example-neighbor

.. literalinclude:: ../../examples/neighbor_distances.py
   :language: python
   :linenos:

预期输出
--------

.. code-block:: text

   FCC Cu: atoms=32, NN distance=2.5750 A, max bin count=5632
   BCC Fe: atoms=16, NN distance=2.4750 A, max bin count=2688
   HCP Mg: atoms=16, NN distance=3.1750 A, max bin count=1376
   wrote: .../docs/_build/example-neighbor/neighbor_distances.png

.. figure:: /_static/images/generated/neighbor_distances.png
   :width: 960px
   :alt: Normalized RDF and cumulative coordination number for FCC Cu, BCC Fe and HCP Mg

   左图：三种结构的归一化 ``g(r)`` 曲线。FCC Cu 第一峰在 ``2.55 Å``
   （= ``a/sqrt(2)``），BCC Fe 第一峰在 ``2.48 Å``
   （= ``a*sqrt(3)/2``），HCP Mg 第一峰在 ``3.21 Å``（= ``a``）。每个
   峰对应一个 coordination shell。右图：累计配位数曲线，归一化到
   cutoff 球内的总数。FCC 的第一壳层配位数为 12，BCC 为 8，HCP 为 12。

结果含义
--------

* FCC Cu 第一近邻距离 ``2.55 Å`` 与 ``a/sqrt(2) = 3.61/sqrt(2)`` 一致，
  对应 12 配位的面心立方密堆。
* BCC Fe 第一近邻距离 ``2.48 Å`` 与 ``a*sqrt(3)/2 = 2.87*sqrt(3)/2``
  一致，对应 8 配位的体心立方。
* HCP Mg 第一近邻距离 ``3.21 Å`` 等于 ``a``，对应 12 配位的密堆六方
  （理想 HCP 与 FCC 第一配位数相同，但次近邻排列不同）。
* 第一峰越尖锐，说明该壳层配位数越大且距离分布越集中。
* ``cutoff = 8 Å`` 对 FCC Cu 大约能覆盖 5--6 个壳层。

参数说明
--------

.. list-table:: :func:`get_neighbor_distances` 关键参数
   :header-rows: 1
   :widths: 22 14 64

   * - 参数
     - 默认
     - 含义
   * - ``atoms``
     - 必填
     - ``ase.Atoms`` 对象。本教程传入 2x2x2 supercell 以避免 cutoff
       超出 cell。
   * - ``cutoff``
     - ``10``
     - 近邻搜索半径，单位 Å。本教程设为 ``8``。
   * - ``self_interaction``
     - ``False``
     - 是否包含 ``i == j`` 的自相互作用。
   * - ``bothways``
     - ``True``
     - 是否双向搜索邻居。``True`` 会同时记录 ``(i, j)`` 和 ``(j, i)``。
   * - ``skin``
     - ``0.0``
     - 邻居列表的 skin 厚度，单位 Å。设为 0 以保证距离精确。

.. list-table:: :func:`build_rdf` 关键参数
   :header-rows: 1
   :widths: 22 14 64

   * - 参数
     - 默认
     - 含义
   * - ``distances``
     - 必填
     - ``get_neighbor_distances`` 返回的 pair distance 数组。
   * - ``bin_width``
     - ``0.05``
     - RDF bin 宽度，单位 Å。
   * - ``r_max``
     - ``CUTOFF``
     - RDF 最大半径，单位 Å。

验证方法
--------

* FCC Cu 第一峰位置必须接近 ``a/sqrt(2) = 2.5527 Å``（容差 ``0.10 Å``）；
* BCC Fe 第一峰位置必须接近 ``a*sqrt(3)/2 = 2.4826 Å``（容差 ``0.10 Å``）；
* HCP Mg 第一峰位置必须接近 ``a = 3.21 Å``（容差 ``0.10 Å``）；
* 三种结构的 RDF 必须至少有一个非零 bin。

常见错误
--------

``ValueError: cell has not been set``
   ``atoms`` 没有 cell。检查 ``bulk(...)`` 是否正确返回；如果用
   ``Atoms(...)`` 手动构造，必须传 ``cell=`` 参数。

RDF 第一峰位置远大于预期
   ``cutoff`` 太小，或 supercell 不够大。规则：``cutoff < L/2``，
   其中 ``L`` 是 supercell 最短边长。

RDF 全为零
   ``bothways=False`` 且 ``self_interaction=False`` 可能漏掉所有 pair。
   保持默认 ``bothways=True``。

下一步
------

* :doc:`atom_manipulation` 看如何基于 ``z`` 坐标筛选子集原子；
* :doc:`surface_energy` 看如何把 neighbor 距离用于表面能收敛性检查；
* :doc:`../manual/vasp` 了解 ``vasp_utils`` 中分析 ``OUTCAR`` 力的脚本。

Related API
-----------

* :func:`mymetal.universal.atom.neighbor.get_neighbor_distances`
* :class:`ase.neighborlist.NeighborList`
* :func:`numpy.histogram`
