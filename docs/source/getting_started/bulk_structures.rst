.. _getting-started-bulk-structures:

构建并对比 FCC / BCC / HCP / Diamond 四种常见 bulk
===================================================

:Audience: 第一次使用 ``mymetal`` 与 ``ase.build.bulk`` 的用户
:Time: 5--10 分钟
:Requires: Python 3.10、pjvasp_package、ASE、Matplotlib
:Runs VASP: No

目标
----

掌握 ASE :func:`ase.build.bulk` 的最常用调用方式，一次性生成 FCC Cu、
BCC Fe、HCP Mg 与 Diamond Si 四个典型 bulk cell，对比它们的 cell 参数、
原子数与角度，并写出可被 VASP 直接使用的 ``POSCAR``。这个案例是后续
surface、slab、k-point 教程的最小入口。

最终结果
--------

.. figure:: /_static/images/generated/bulk_structures.png
   :width: 960px
   :alt: 3D and side views of FCC Cu, BCC Fe, HCP Mg and Diamond Si conventional cells

   上排是四个 bulk cell 的 3D 倾斜视图，下排是沿 ``y`` 方向观察的侧视图。
   ``FCC Cu`` 和 ``BCC Fe`` 都用 ``cubic=True`` 给出正交 conventional
   cell；``HCP Mg`` 用 primitive cell 直接显示 120° ``gamma`` 角；
   ``Diamond Si`` 用 ``cubic=True`` 给出 8 原子金刚石 conventional cell。

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
           └── bulk_structures.py

运行完整示例
------------

:download:`下载 bulk_structures.py
<../../examples/bulk_structures.py>`，或在仓库根目录直接运行：

.. code-block:: console

   $ python docs/examples/bulk_structures.py --output docs/_build/example-bulk

页面与用户下载复用同一份实现：

.. literalinclude:: ../../examples/bulk_structures.py
   :language: python
   :linenos:

预期终端输出
------------

.. code-block:: text

   label                      formula  atoms  a_A     b_A     c_A     alpha   beta    gamma   volume_A3  pbc
   FCC Cu (cubic)              Cu4      4      3.6100  3.6100  3.6100  90.00   90.00   90.00   47.0459    [True, True, True]
   BCC Fe (cubic)              Fe2      2      2.8700  2.8700  2.8700  90.00   90.00   90.00   23.6399    [True, True, True]
   HCP Mg (primitive)          Mg2      2      3.2100  3.2100  5.2100  90.00   90.00  120.00   46.4920    [True, True, True]
   Diamond Si (cubic)          Si8      8      5.4300  5.4300  5.4300  90.00   90.00   90.00  160.1030    [True, True, True]
   wrote: .../docs/_build/example-bulk/bulk_structures.png

输出目录为：

.. code-block:: text

   docs/_build/example-bulk/
   ├── POSCAR_FCC
   ├── POSCAR_BCC
   ├── POSCAR_HCP
   ├── POSCAR_Diamond
   └── bulk_structures.png

如何解释结果
------------

* ``FCC Cu`` 用 ``cubic=True`` 得到 4 原子 conventional cell（顶点 + 面心），
  体积 ``a^3 = 47.05 Å³``；这是 VASP ``POSCAR`` 中最常见的 FCC 输入。
* ``BCC Fe`` 用 ``cubic=True`` 得到 2 原子 conventional cell（顶点 + 体心），
  体积 ``a^3 = 23.64 Å³``；``a = 2.87 Å`` 是室温 α-Fe 文献值。
* ``HCP Mg`` 用 primitive cell，2 个原子，``gamma = 120°``；``c = 5.21 Å``
  与 ``a = 3.21 Å`` 给出 ``c/a = 1.623``，接近理想 HCP ``sqrt(8/3) = 1.633``。
* ``Diamond Si`` 用 ``cubic=True`` 得到 8 原子金刚石 cell，每个原子 4 个
  最近邻；``a = 5.43 Å`` 是 Si 实验值。

参数说明
--------

.. list-table:: :func:`ase.build.bulk` 关键参数
   :header-rows: 1
   :widths: 22 14 64

   * - 参数
     - 默认
     - 含义
   * - ``symbol``
     - 必填
     - 元素符号，例如 ``"Cu"``。ASE 会自动查找对应的实验晶格常数。
   * - ``structure``
     - 必填
     - 晶体类型。本教程依次使用 ``"fcc"``、``"bcc"``、``"hcp"``、
       ``"diamond"``。ASE 还支持 ``"sc"``（简单立方）、``"rocksalt"`` 等。
   * - ``a``
     - ``None``
     - 立方晶系的晶格常数，单位 Å。提供时覆盖实验值。
   * - ``c``
     - ``None``
     - HCP 的 ``c`` 轴长度，单位 Å。HCP 必须额外提供 ``covera`` 或 ``c``。
   * - ``covera``
     - ``None``
     - HCP 的 ``c/a`` 比例。本教程设为 ``5.21/3.21``。
   * - ``cubic``
     - ``False``
     - 是否返回 conventional cell（而非 primitive cell）。FCC 返回 4 原子
       立方 cell，BCC 返回 2 原子立方 cell，Diamond 返回 8 原子立方 cell。

.. list-table:: :func:`mymetal.io.vasp.my_write_vasp` 关键参数
   :header-rows: 1
   :widths: 22 14 64

   * - 参数
     - 默认
     - 含义
   * - ``filename``
     - 必填
     - 输出 ``POSCAR`` 路径。
   * - ``atoms``
     - 必填
     - 待写出的 ``ase.Atoms`` 对象。
   * - ``label``
     - ``None``
     - ``POSCAR`` 第一行注释。
   * - ``lattice_scale_factor``
     - ``1``
     - 全局 scale factor。``1`` 表示 cell 矢量直接以 Å 写出。
   * - ``direct``
     - ``True``
     - ``True`` 写分数坐标（``Direct``），``False`` 写笛卡尔坐标。
   * - ``vasp5``
     - ``True``
     - 是否输出 VASP 5+ 格式（带元素行）。

验证方法
--------

脚本在写出图片前执行以下断言：

* 四个 cell 的 ``atoms > 0`` 且 ``pbc == [True, True, True]``；
* FCC / BCC / Diamond 三个 cubic cell 的 ``alpha``/``beta``/``gamma`` 都
  为 90°；
* HCP primitive cell 的 ``gamma`` 为 120°（在 ``1e-3`` 容差内）；
* 每个 ``POSCAR`` 都能被 :func:`mymetal.io.vasp.my_read_vasp` 读回，且
  cell / positions 在 ``1e-10`` 容差内与原始 ``Atoms`` 一致。

常见错误
--------

``ValueError: covera must be supplied for hcp``
   HCP 晶系必须显式提供 ``covera`` 或同时提供 ``a`` 和 ``c``。ASE 不会
   使用元素数据库中的 ``c/a`` 值。

POSCAR 中元素行缺失
   ``my_write_vasp`` 默认 ``vasp5=True``。如果 ``POSCAR`` 第一行没有元素
   符号，检查是否手动传了 ``vasp5=False``。

看到 ``Diamond Si`` 有 8 个原子
   ``cubic=True`` 给出 conventional cell。改用 ``cubic=False`` 会返回 2
   原子 primitive cell，但坐标在菱方 cell 中，可能不符合 VASP 习惯。

下一步
------

* :doc:`au111_slab` 学习如何把 bulk 切成 slab；
* :doc:`fcc_surfaces` 同时构建 Cu(100)/(110)/(111)；
* :doc:`../tutorials/kpoints_sampling` 看如何为这些 bulk 选择 k 点。

Related API
-----------

* :func:`ase.build.bulk`
* :func:`mymetal.io.vasp.my_write_vasp`
* :func:`mymetal.io.vasp.my_read_vasp`
