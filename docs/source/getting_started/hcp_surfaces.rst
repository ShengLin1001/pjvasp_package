.. _getting-started-hcp-surfaces:

构建并对比 HCP Mg 低指数表面
=============================

:Audience: 想用 ``mymetal`` 构建 HCP 体系 basal / prism I / prism II 表面的用户
:Time: 10 分钟
:Requires: Python 3.10、pjvasp_package、ASE、Matplotlib
:Runs VASP: No

目标
----

掌握 :mod:`mymetal.build.bulk.create` 中针对 HCP 体系的四个构造函数：
:func:`create_hcp_basal`、:func:`create_hcp_prism1`（含 wide/narrow 两种
层间距）和 :func:`create_hcp_prism2`，构建 Mg(0001)、Mg(10-10) wide、
Mg(10-10) narrow 与 Mg(10-11) 四个 slab，对比面内 cell、原子数和层序列。

最终结果
--------

.. figure:: /_static/images/generated/hcp_surfaces.png
   :alt: Side and top views of Mg(0001), Mg(10-10) prism I (wide/narrow), Mg(10-11) prism II slabs

   上排是四个 slab 沿面内 2 倍重复后的侧视图，下排是 2x2 重复后的俯视图。
   四个 slab 共享同一 ``a = 3.21 Å``、``c = 5.21 Å`` （Mg 文献值）。图片由
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
           └── hcp_surfaces.py

运行完整示例
------------

:download:`下载 hcp_surfaces.py <../../examples/hcp_surfaces.py>`，或在仓库
根目录直接运行：

.. code-block:: console

   $ python docs/examples/hcp_surfaces.py --output docs/_build/example-hcp-surfaces

页面与用户下载复用同一份实现：

.. literalinclude:: ../../examples/hcp_surfaces.py
   :language: python
   :linenos:

预期终端输出
------------

.. code-block:: text

   label                          formula  atoms  a_A      b_A      c_total_A  pbc
   Mg(0001) basal                  Mg16     16     3.2100   3.2100   69.0750    [True, True, True]
   Mg(10-10) prism I (wide)        Mg20     20     3.2100   5.2100   55.9461    [True, True, True]
   Mg(10-10) prism I (narrow)      Mg20     20     3.2100   5.2100   56.8728    [True, True, True]
   Mg(10-11) prism II              Mg24     24     5.5599   5.2100   47.6550    [True, True, True]
   wrote: .../docs/_build/example-hcp-surfaces/hcp_surfaces.png

输出目录为：

.. code-block:: text

   docs/_build/example-hcp-surfaces/
   ├── POSCAR_Mg_0001
   ├── POSCAR_Mg_10_10
   ├── POSCAR_Mg_10_10
   ├── POSCAR_Mg_10_11
   └── hcp_surfaces.png

如何解释结果
------------

* 四个 slab 的 ``formula`` 都以 ``Mg`` 开头，原子数随 ``size`` 重复数线性
  增加：basal 是 2 atoms/cell × 8 重复 = 16，prism I/II 是 4 atoms/cell ×
  5 或 6 重复。
* basal 的 ``a_A`` 和 ``b_A`` 都等于 ``3.21`` Å，对应 (0001) 面内六角
  primitive cell。prism I 把 ``b`` 方向变成 ``5.21`` Å（= ``c``），prism II
  把 ``a`` 方向变成 ``5.56`` Å（= ``sqrt(3)*a``）。
* ``c_total_A`` 包括 slab 厚度加上两侧各 15 Å 真空；脚本通过
  ``cell.center(vacuum=15.0, axis=2)`` 在切面后加真空。
* prism I 的 ``wide`` 和 ``narrow`` 模式对应同一 (10-10) 面但不同的层
  序列：wide 的层间距为 ``a/sqrt(3) ≈ 1.85 Å``，narrow 为其一半。
  HCP 的 prism I slip system 在两种 termination 上的临界剪应力不同，因此
  需要分别建模。

参数说明
--------

.. list-table:: :func:`create_hcp_basal` / :func:`create_hcp_prism1` /
   :func:`create_hcp_prism2` 在本教程使用的关键参数
   :header-rows: 1
   :widths: 18 12 70

   * - 参数
     - 默认
     - 含义
   * - ``a``
     - ``None``
     - HCP ``a`` 晶格常数，单位 Å。本教程设为 ``3.21`` （Mg 文献值）。
   * - ``c``
     - ``None``
     - HCP ``c`` 晶格常数，单位 Å。本教程设为 ``5.21``。
   * - ``size``
     - ``(1,1,1)``
     - 沿 (a1, a2, a3) 的超胞重复数。本教程调整 ``size[2]`` 控制 slab
       厚度。
   * - ``symbol``
     - ``"Au"``
     - 元素符号。本教程统一设为 ``"Mg"``。
   * - ``pbc``
     - ``(1,1,1)``
     - 三个方向的周期边界条件。
   * - ``mode`` (仅 prism1)
     - ``"wide"``
     - prism I 的 termination：``"wide"`` 或 ``"narrow"``，对应不同的
       层间距。

验证方法
--------

脚本在写出图片前执行以下断言：

* 四个 slab 的 ``formula`` 都以 ``Mg`` 开头，``atoms > 0``；
* 四个 slab 的 ``pbc`` 都为 ``[True, True, True]``；
* 四个 ``POSCAR`` 都能被 :func:`mymetal.io.vasp.my_read_vasp` 读回，且
  cell/positions 在 ``1e-10`` 容差内与原始 ``Atoms`` 一致；
* PNG 文件存在且包含足够的非白像素。

常见错误
--------

``vasp_create_hcp_prism1`` 报错或结构异常
   ``create_hcp_prism1`` 内部调用 ``make_SFP_xy`` 和 ``make_a3_ortho``
   完成 cell 正交化。这两个 helper 来自 companion 包 ``myvasp``，使用前
   确认已安装 ``myalloy_package``。

prism I 的 wide 和 narrow 看起来一样
   侧视图的差异很小（层间距差两倍）。查看 ``c_total_A`` 列：相同
   ``size=(1,1,5)`` 下 narrow 的 c 长度比 wide 略大，因为 narrow 模式
   在 slab 顶部留出了额外的 half-layer 间距。

看不到真空
   ``create_hcp_*`` 函数本身不加真空。本教程脚本在调用它们之后用
   ``cell.center(vacuum=15.0, axis=2)`` 补上真空。

下一步
------

* :doc:`au111_slab` 看 ``generate_film`` 的另一种用法（带 primitive
  cell 化简）；
* :doc:`../tutorials/biaxial_stretch` 看如何对 slab 施加单轴/双轴应变；
* :doc:`../manual/vasp` 了解 VASP workflow 的完整 lifecycle。

Related API
-----------

* :func:`mymetal.build.bulk.create.create_hcp_basal`
* :func:`mymetal.build.bulk.create.create_hcp_prism1`
* :func:`mymetal.build.bulk.create.create_hcp_prism2`
* :func:`mymetal.io.vasp.my_read_vasp`
* :func:`mymetal.io.vasp.my_write_vasp`
