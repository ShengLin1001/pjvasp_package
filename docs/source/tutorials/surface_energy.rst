.. _tutorial-surface-energy:

用已有结果计算 surface energy
================================

:Audience: 已有 bulk/slab VASP 结果的用户
:Time: 10 分钟
:Requires: Python、ASE、pjvasp_package、仓库内 tracked fixture
:Runs VASP: No

目标与边界
----------

从仓库已有 Au FCC(111) bulk/slab 结构和能量复算 surface energy，并把手算结果
与 :func:`mymetal.calculate.calenergy.surfenergy.cal_surface_energy` 对照。
这是一条 parser/API 工程示例，不是独立复核过的材料学 benchmark。

.. warning::

   仓库当前没有顶层 LICENSE 文件。本教程不复制或提供 fixture 下载，只读取同一
   clone 中已经跟踪的文件并展示派生结果。把结构或 VASP 输出重新打包到公开下载
   资产前，维护者仍需确认数据再分发条款。

公式与单位
----------

当前 API 实现：

.. math::

   \gamma =
   \frac{E_\mathrm{slab}
   - (E_\mathrm{bulk}/N_\mathrm{bulk})N_\mathrm{slab}}
   {f A}.

输入能量均为 eV，``A`` 为一个表面的面积（Å²），``f`` 是 excess energy
对应的等价表面数量。``energy_unit="eV"`` 返回 eV/Å²；
``energy_unit="J"`` 把同一结果转换为 J/m²。

.. important::

   本 fixture 沿用仓库 notebook 的 ``factor=2``，即把 excess energy 分给两个
   等价表面。非对称 slab、不同 termination 或 adsorption model 不能机械复用 2。

输入与 provenance
------------------

从仓库根目录使用以下已跟踪文件：

.. code-block:: text

   mymetal/example/test-surface-energy/fcc/
   ├── energy.txt
   └── 1.000-2.8485/
       ├── bulk/CONTCAR
       └── full_relaxed_surface_111/CONTCAR

本案例在 ``energy.txt`` 中选择这一对：

.. list-table:: Fixture inputs
   :header-rows: 1

   * - Quantity
     - Value
     - Unit
   * - ``E_bulk``
     - -47.04775262
     - eV
   * - ``N_bulk``
     - 12
     - atoms
   * - ``E_slab``
     - -46.24580704
     - eV
   * - ``N_slab``
     - 12
     - atoms
   * - ``factor``
     - 2
     - equivalent surfaces

SHA-256：

.. code-block:: text

   27D8C2A4999C5B11BA2626AD108A27360DD75AA2896866DEC154816CCFF5A50B  fcc/energy.txt
   C2C35710DF4FB43F2BD13EEEF972177E109B31BE8E43DDDA0079D4957B94A3A1  bulk/CONTCAR
   61455310E168E32A7450225677366F3505F5632CA66D846E75D5D93852CD3B9A  full_relaxed_surface_111/CONTCAR

运行完整脚本
------------

.. code-block:: console

   $ python docs/examples/surface_energy.py

.. literalinclude:: ../../examples/surface_energy.py
   :language: python
   :linenos:

预期输出
--------

.. code-block:: text

   bulk atoms: 12
   slab atoms: 12
   area_A2: 7.0268887736
   factor: 2
   manual_eV_A2: 0.0570626351
   api_eV_A2: 0.0570626351
   api_J_m2: 0.9142441870
   check: ok

手算/API 对照
-------------

``cal_area`` 实际返回 ``abs((a × b)[2])``。本 fixture 的 a、b 都在 xy 平面，
因此 projected area 就是完整表面面积 ``7.0268887736 Å²``。

.. list-table:: Calculated result
   :header-rows: 1

   * - Method
     - eV/Å²
     - J/m²
   * - Direct formula
     - 0.0570626351
     - —
   * - ``cal_surface_energy``
     - 0.0570626351
     - 0.9142441870

脚本用 ``1e-12`` tolerance 断言手算和 API 一致，并对仓库已保存的 eV/Å²
结果做回归检查。J/m² 数值使用当前代码中的转换常数 ``16.021766``，因此可能与
旧 notebook 采用更多有效位的常数在最后几位略有差异。

结果含义
--------

这里确认的是：给定这组结构、能量、面积与 ``factor``，当前 API 能稳定复算仓库
记录。它不证明 DFT 设置、bulk reference、slab thickness 或 surface termination
足以支持新的科研结论；这些仍需领域维护者复核。

常见错误
--------

``area must be a positive number``
   检查 cell 的前两个矢量；``cal_area`` 只计算 xy projection，不适用于任意倾斜
   表面而不做额外几何判断。

结果大约相差两倍
   首先核对 ``factor``，再核对 slab 是否有两个等价表面；不要只为匹配目标数值而
   修改 ``factor``。

结果单位异常
   输入能量必须是 eV、面积必须是 Å²。``energy_unit`` 只选择输出单位，不表示输入
   能量已经是 joule。

找不到 fixture
   在完整仓库 clone 的根目录运行，或用 ``-path_case`` 指向包含 ``bulk/`` 和
   ``full_relaxed_surface_111/`` 的 case。

下一步
------

继续阅读 :doc:`outcar_batch`，或在 :doc:`../manual/vasp` 中了解如何把已有
VASP 目录接入后处理。真实 VASP 计算仍需用户自行准备合法 POTCAR。

Related API
-----------

* :func:`mymetal.build.film.extrfilm.cal_area`
* :func:`mymetal.calculate.calenergy.surfenergy.cal_surface_energy`
