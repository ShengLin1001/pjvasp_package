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

``mymetal.build.bulk.create`` 还包含 ``create_fcc_111`` 与
``create_hcp_basal``。heterostructure 入口依赖可选 ``hetbuilder``，本轮未把
无法在当前环境实测的路径包装成正式教程。

.. automodule:: mymetal.build.workflow.hoec
   :members:

.. automodule:: mymetal.build.workflow.kpar_ncore
   :members:
