Structure Building
==================

.. automodule:: mymetal.build

Film construction
-----------------

``generate_film``
~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.build.film.stretch.generate_film

用途
   从内部 FCC/HCP bulk 或传入的 ``ase.Atoms`` 切取指定 Miller plane，
   生成 slab/periodic cell，并可加真空和 primitive-cell 处理。

关键参数
   ``symbols``/``structure`` 定义内部 bulk；``num_layers`` 优先于
   ``replic_z``；``my_vacuum`` 单位为 Å；``slice_plane`` 是
   ``(h, k, l)``；``a_fcc``/``a_hcp`` 单位为 Å；``move_atom`` 与
   ``before_move_atom`` 是 fractional translation。

返回
   一个 ``ase.Atoms``。函数会执行 ``wrap()``；PBC 和真空行为由
   ``my_periodic``/``my_vacuum`` 共同决定。

异常
   非 ``fcc``/``hcp`` 且没有 ``bulk_atoms`` 时抛 ``ValueError``；
   ``num_layers`` 与 ``replic_z`` 都未给出时也抛 ``ValueError``。
   ASE/spglib 的几何错误会继续向上传播。

最小示例
   .. code-block:: python

      from mymetal.build.film.stretch import generate_film

      slab = generate_film(
          symbols="Au", structure="fcc", num_layers=12,
          my_vacuum=20.0, slice_plane=(1, 1, 1),
      )
      assert slab.get_chemical_formula() == "Au12"

Related tutorials
   :doc:`../getting_started/au111_slab`

``cal_area``
~~~~~~~~~~~~

.. autofunction:: mymetal.build.film.extrfilm.cal_area

用途与单位
   返回第一、第二 cell vector 叉积的 z 分量绝对值，即 xy projected area，
   单位 Å²。只有这两个矢量位于 xy 平面时，它才等于完整 parallelogram area。

参数/返回
   ``atoms`` 是带非退化 cell 的 ``ase.Atoms``；返回 ``float``/NumPy scalar。
   ``atoms=None`` 或无有效 cell 时，底层属性/数组错误会向上传播。

最小示例
   .. code-block:: python

      from mymetal.build.film.extrfilm import cal_area

      area_a2 = cal_area(slab)
      assert area_a2 > 0

Related tutorials
   :doc:`../tutorials/surface_energy`

其他结构入口
------------

``mymetal.build.bulk.create`` 提供 ``(hkl)`` 取向的 FCC/HCP 晶胞构造器，
``mymetal.build.bulk.gsfe`` 在此基础上生成 GSFE（广义层错能）模型，
``mymetal.build.film.findcubic`` 把 primitive 薄膜 cell 转成正交 cell，
``mymetal.build.film.stretch`` 还提供沿指定方向的拉伸函数。

``create_fcc_111``
~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.build.bulk.create.create_fcc_111

用途
   构造 FCC(111) 取向晶胞，使用 ASE ``FaceCenteredCubic`` + direction/miller
   方式，并做微小 in-plane 平移避免原子落在倾斜 cell 的 xlo/ylo 面上。

关键参数
   ``a`` 为晶格常数（Å，``None`` 时用 ASE 默认值）；``size`` 为
   ``(nx, ny, nz)`` 超胞重复；``symbol`` 为元素符号；``pbc`` 为三方向
   周期性。层间距约为 ``a/sqrt(3)``。

返回
   一个 ``ase.Atoms``，已 ``wrap()``。

最小示例
   .. code-block:: python

      from mymetal.build.bulk.create import create_fcc_111

      atoms = create_fcc_111(a=4.08, size=(1, 1, 7), symbol='Au')
      assert len(atoms) > 0

Related tutorials
   :doc:`../tutorials/gsfe_models`

``create_hcp_basal``
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.build.bulk.create.create_hcp_basal

用途
   构造 HCP(0001) 基面晶胞，使用 ASE ``HexagonalClosedPacked``。
   层间距为 ``c/2``。同样做 in-plane 微平移避免角原子在 LAMMPS
   minimize 时被推出 box。

关键参数
   ``a``、``c`` 为 HCP 晶格常数（Å）；``size``、``symbol``、``pbc``
   同上。

Related tutorials
   :doc:`../tutorials/gsfe_models`

``create_hcp_prism1``
~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.build.bulk.create.create_hcp_prism1

用途
   构造 HCP(10-10) prism I 晶胞，支持 ``mode='wide'``（宽层，间距
   ``a/sqrt(3)``）和 ``mode='narrow'``（窄层，间距 ``a/(2*sqrt(3))``）。

``create_gsfe_model``
~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.build.bulk.gsfe.create_gsfe_model

用途
   根据 ``gsfe_type`` 选择晶体取向，生成用于 GSFE 计算的超胞。
   支持的类型：``"FCC_111"``、``"FCC_100"``、``"HCP_basal"``、
   ``"HCP_prism1w"``、``"HCP_pyr1w"``、``"HCP_pyr2"``。

关键参数
   ``a``、``c``（Å）为晶格常数；``size`` 为超胞尺寸。未给 ``size``
   时使用各类型的内置默认值。

.. note::

   ``"HCP_pyr1w"``、``"HCP_pyr2"``、``"FCC_100"`` 路径调用
   ``vasp_create_*`` 系列函数，依赖可选的 ``myvasp`` 包。
   ``"FCC_111"``、``"HCP_basal"``、``"HCP_prism1w"`` 路径仅依赖 ASE。

最小示例
   .. code-block:: python

      from mymetal.build.bulk.gsfe import create_gsfe_model

      atoms = create_gsfe_model(gsfe_type="FCC_111", a=4.08, c=None, size=[1, 1, 7])
      assert len(atoms) > 0

Related tutorials
   :doc:`../tutorials/gsfe_models`

``find_cubic``
~~~~~~~~~~~~~~

.. autofunction:: mymetal.build.film.findcubic.find_cubic

用途
   把 primitive hcp/fcc 薄膜 cell（输入 ``[a,b,c,any,any,90]``）
   转成 xy 面正交 cell：``cell[0,0] = sqrt(area/sqrt(3))``，
   ``cell[1,1] = sqrt(3)*cell[0,0]``，并清零 off-diagonal xy 分量。

参数
   ``prim`` 为输入 ``ase.Atoms``；``type`` 目前仅支持 ``'hcp'``/``'fcc'``。

返回
   正交化后的 ``ase.Atoms`` （深拷贝，不修改输入）。

最小示例
   .. code-block:: python

      from mymetal.build.film.findcubic import find_cubic
      from mymetal.build.film.stretch import generate_film

      film = generate_film(symbols='Mg', structure='hcp', num_layers=4,
                           my_vacuum=10.0, slice_plane=(0,0,1), a_hcp=3.21)
      ortho = find_cubic(film, type='hcp')

Related tutorials
   :doc:`../tutorials/cubic_cell_and_stretch`

``find_optimal_cell_shape``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.build.film.findcubic.find_optimal_cell_shape

用途
   搜索整数变换矩阵 ``P``，使 ``P @ cell`` 最接近目标形状（``'sc'``
   简单立方或 ``'fcc'`` 面心立方）。基于 ASE 实现，加入 ``doped`` 包
   的旋转不变改进。

``stretch_along_direction_to_cell``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.build.film.stretch.stretch_along_direction_to_cell

用途
   沿 ``x`` / ``y`` / ``z`` 方向按给定拉伸因子列表拉伸 cell，原子 fractional
   坐标不变。返回拉伸后的 ``ase.Atoms`` （单个因子）或列表。

关键参数
   ``stretch_direction_list`` 如 ``['x']``；``stretch_factor_list``
   如 ``[0.99, 1.00, 1.01]``。

最小示例
   .. code-block:: python

      from mymetal.build.film.stretch import stretch_along_direction_to_cell

      stretched = stretch_along_direction_to_cell(
          film, stretch_direction_list=['x'], stretch_factor_list=[1.01])

Related tutorials
   :doc:`../tutorials/cubic_cell_and_stretch`、
   :doc:`../tutorials/biaxial_stretch`

.. automodule:: mymetal.build.workflow.hoec
   :members:

.. automodule:: mymetal.build.workflow.kpar_ncore
   :members:
