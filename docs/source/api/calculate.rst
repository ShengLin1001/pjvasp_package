Calculations
============

.. automodule:: mymetal.calculate

Surface energy
--------------

.. automodule:: mymetal.calculate.calenergy.surfenergy

.. autofunction:: mymetal.calculate.calenergy.surfenergy.cal_surface_energy

Contract
~~~~~~~~

* ``bulk_energy`` 与 ``relaxed_surface_energy``：eV；
* ``bulk_atoms_number`` 与 ``surface_atoms_number``：正整数；
* ``area``：一个表面的面积，Å²，必须大于 0；
* ``factor``：excess energy 对应的等价表面数，正整数；
* ``energy_unit="eV"`` 返回 eV/Å²；``"J"`` 返回 J/m²。输入能量始终按
  eV 解释。

函数实现的公式为：

.. math::

   \frac{E_\mathrm{slab}
   - E_\mathrm{bulk}N_\mathrm{slab}/N_\mathrm{bulk}}
   {fA}.

返回 ``float``/NumPy scalar。非法 atom count、area、factor 或
``energy_unit`` 抛 ``ValueError``；无法相除的非数值输入会抛 Python 原生
``TypeError``。

最小示例
~~~~~~~~

.. code-block:: python

   from mymetal.calculate.calenergy.surfenergy import cal_surface_energy

   gamma = cal_surface_energy(
       bulk_energy=-47.04775262,
       bulk_atoms_number=12,
       relaxed_surface_energy=-46.24580704,
       surface_atoms_number=12,
       area=7.026888773593904,
       factor=2,
   )
   assert abs(gamma - 0.05706263510343279) < 1e-12

Related tutorials
~~~~~~~~~~~~~~~~~

See :doc:`../tutorials/surface_energy`.

Mechanics and mismatch
----------------------

.. automodule:: mymetal.calculate.calmechanics.deformation
   :members:
   :show-inheritance:

.. automodule:: mymetal.calculate.calmechanics.strain
   :members:
   :show-inheritance:

.. automodule:: mymetal.calculate.calmechanics.stress
   :members:
   :show-inheritance:

``cal_stretch`` / ``cal_relative_stretch``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.calculate.calmechanics.stretch.cal_stretch

.. autofunction:: mymetal.calculate.calmechanics.stretch.cal_relative_stretch

用途
   基于 ``hetbuilder.Interface`` 的超胞匹配结果，计算 bottom/top 层相对于
   heterostructure stack 的主应变方向拉伸因子。``cal_relative_stretch``
   返回 ``[factor - 1, direction]``，``cal_stretch`` 返回绝对值。

.. note::

   这两个函数依赖可选的 ``hetbuilder`` 包。

Linear algebra (calmath)
------------------------

``hermite_normal_form``
~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.calculate.calmath.matrix.hermite_normal_form

用途
   计算整数矩阵的 Hermite Normal Form (HNF)，用于分析超胞变换矩阵、
   寻找 commensurate cell。使用纯整数算术保证结果为整数矩阵。

参数
   ``matrix`` 为 ``np.ndarray``（会被转为 ``int64``）。

返回
   HNF 矩阵（``np.ndarray``，``int64``）。

最小示例
   .. code-block:: python

      import numpy as np
      from mymetal.calculate.calmath.matrix import hermite_normal_form

      m = np.array([[2, 1, 0], [0, 2, 1], [0, 0, 3]])
      hnf = hermite_normal_form(m)

Related tutorials
   :doc:`../tutorials/deformation_and_hnf`

Mismatch analysis (calmismatch)
-------------------------------

.. automodule:: mymetal.calculate.calmismatch.calhetero

Principal entry point: ``mymetal.calculate.calmismatch.calhetero.cal_mismatch``。

``cal_mismatch``
~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.calculate.calmismatch.calhetero.cal_mismatch

用途
   计算 bottom、top 与 heterostructure stack 三者之间的最大 mismatch，
   返回 ``[a_mismatch, b_mismatch, gamma_mismatch]``。每个分量为
   ``max(|bottom-hetero|/bottom, |top-hetero|/top)``。

.. note::

   ``calhetero`` 模块依赖可选的 ``hetbuilder`` 包（``Interface`` 类）。
   ``cal_mismatch`` 本身只依赖 ASE，可独立使用。

``compare_atoms``
~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.calculate.calmismatch.calhetero.compare_atoms

用途
   逐项比较两个 ``ase.Atoms`` 的 cell 长度、角度、笛卡尔坐标和 fractional
   坐标，返回 4 个 bool 的列表。

``cal_stretch_lattice``
~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.calculate.calmismatch.calhetero.cal_stretch_lattice

用途
   根据 bottom/top 参考层和 heterostructure stack 的 cell 参数，计算
   两层各自的 ``[a, b, gamma]`` 拉伸量 ``(hetero - ref) / ref``。

Electronic and materials analysis
---------------------------------

``cal_reciprocal_matrix`` / ``cal_reciprocal_matrix2``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.calculate.calqm.kpoints.cal_reciprocal_matrix
   :no-index:

.. autofunction:: mymetal.calculate.calqm.kpoints.cal_reciprocal_matrix2
   :no-index:

用途
   两种方法计算倒格子矢量：``cal_reciprocal_matrix`` 用叉积公式，
   ``cal_reciprocal_matrix2`` 用矩阵求逆。两者结果一致（在数值精度内）。

Related tutorials
   :doc:`../tutorials/reciprocal_lattice`、 :doc:`../tutorials/kpoints_sampling`

``get_size_by_distance``
~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.calculate.calqm.kpoints.get_size_by_distance
   :no-index:

用途
   根据 RK 乘积或 KSPACING 计算 VASP 自动 k 点网格。RK 法用
   ``round``，KSPACING 法用 ``ceil``。返回 ``(old_auto_kpoints, new_kspacing_kpoints)``。

Related tutorials
   :doc:`../tutorials/reciprocal_lattice`、 :doc:`../tutorials/kpoints_sampling`

.. automodule:: mymetal.calculate.calqm.kpoints
   :members:
   :show-inheritance:

.. automodule:: mymetal.calculate.material_science.schmid

Higher-order elastic constants (HOEC)
-------------------------------------

.. automodule:: mymetal.calculate.calmechanics.hoec
   :members:
   :show-inheritance:

用途
   基于 Wang & Li (PRB 79, 224102, 2009) 的 energy-strain 方法，计算任意点群
   对称性的二阶/三阶/四阶弹性常数。核心思想：应变能密度按
   engineering-Voigt Green-Lagrange 应变 ``eta`` 展开，点群对称性把完全对称
   的 Brugger 常数缩减为独立集合（立方：3/6/11，六方：5/10/19）。

关键概念
   * ``eta`` 是 engineering-Voigt 应变（``eta4=2*eps_yz`` 等）；
   * 单参数变形模式 ``eta = xi * d`` 下，``[U(xi)-U(0)]/V0 = P2*xi² + P3*xi³ + P4*xi⁴``；
   * 每个 ``P_n`` 是独立 ``n`` 阶常数的固定线性组合（``B_n_k`` 系数）；
   * ``HOECModel`` 缓存每个对称性的缩减模型，``get_model`` 返回缓存实例；
   * ``get_deformation_gradient`` 给出 ``xi * d`` 对应的对称变形梯度 ``F``；
   * ``get_mode_strain_lists`` 按各模式 severity 预算分配 ``xi`` 列表。

.. note::

   此模块仅依赖 ``numpy`` 和 Python 标准库，不需要 VASP 或外部可执行文件。
   ``selftest_hoec`` 可验证立方对称性下的系数是否与 Wang-Li Table I 一致。

最小示例
   .. code-block:: python

      from mymetal.calculate.calmechanics.hoec import get_model, get_hoec_modes, get_strain_list

      model = get_model('cubic')
      modes = get_hoec_modes('cubic')
      xi_list = get_strain_list(window=0.06)

Related tutorials
   :doc:`../tutorials/strain_deformation`、
   :doc:`../tutorials/post_hoec_energy`

Electronic structure
--------------------

.. automodule:: mymetal.calculate.electronic_structure.universal
   :members:
   :show-inheritance:

.. note::

   ``electronic_structure`` 模块依赖可选的 ``pymatgen`` 包。
   ``plotter.py`` 提供布里渊区可视化（``plot_brillouin_zone_from_kpath``
   等），``universal.py`` 提供能带数据提取（``get_n_band``、
   ``summarize_band_structure_info``、``get_n_kpoints_band``）。
