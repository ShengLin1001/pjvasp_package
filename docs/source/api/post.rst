Post-Processing
===============

.. automodule:: mymetal.post

.. automodule:: mymetal.post.newmain

Stable parser classes used by the tutorial
------------------------------------------

.. autoclass:: mymetal.post.newmain.PostTime
   :members: read_OUTCAR

.. autoclass:: mymetal.post.newmain.PostData
   :members: read_OUTCAR

.. autoclass:: mymetal.post.newmain.PostData2
   :members: read_OUTCAR

Shared contract
~~~~~~~~~~~~~~~

构造参数 ``my_path`` 是 case 根目录，``post_path`` 是将被创建/覆盖的摘要文件，
``config`` 可替换 parser marker。``read_OUTCAR`` 的 ``job_list`` 虽然默认
``None``，当前实现会直接迭代它，因此调用方必须传入明确 list。

三个 ``read_OUTCAR`` 都返回 ``None`` 并写文本：

.. list-table::
   :header-rows: 1

   * - Class
     - Extracted fields
     - Units
   * - ``PostTime``
     - convergence, iteration, elapsed, CPU groups, memory, LOOP real time
     - s, kB for reported memory field
   * - ``PostData``
     - ``energy(sigma->0)``, EENTRO, six stress components
     - eV, eV, kB
   * - ``PostData2``
     - volume, external pressure, maximum force
     - Å³, kB, eV/Å

缺少 case file 时类会打印 warning 并继续；``job_list=None`` 会抛
``TypeError``；输出路径不可写时传播 ``OSError``；marker 不匹配时可能产生空字段。

最小示例
~~~~~~~~

.. code-block:: python

   from mymetal.post.newmain import PostData

   parser = PostData(
       my_path="./y_dir/",
       post_path="./p_post_data.txt",
   )
   parser.read_OUTCAR(job_list=["0.997", "1.000"])

.. warning::

   ``post_general()`` 当前不是推荐入口：它的默认调用路径不能为这些 class 提供
   必需的 ``job_list``。先使用上面三个 parser class 或
   ``docs/examples/outcar_summary.py``。

Related tutorials
~~~~~~~~~~~~~~~~~

See :doc:`../tutorials/outcar_batch`.

Runtime-dependent modules
-------------------------

The modules below are part of the post-processing API, but they import the
external ``myvasp`` helper package and are not imported during the online
documentation build:

* ``mymetal.post.general``
* ``mymetal.post.neb``
* ``mymetal.post.gsfe``
* ``mymetal.post.stretch``
* ``mymetal.post.Cij_energy``
* ``mymetal.post.convergence``

Elastic constants (Cij_energy)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: mymetal.post.Cij_energy
   :members:
   :show-inheritance:

Principal entry point: ``mymetal.post.Cij_energy.fit_cij_energy``。

用途
   用 energy-strain 方法拟合二阶弹性常数 ``Cij``。从 ``y_dir/<strain>``
   目录读取变形能量，对 ``U/V0`` vs ``eta²`` 做二次多项式拟合，提取
   ``Cij`` 分量。配套 ``read_cij_energy``/``write_cij_energy`` 读写结果，
   ``read_deform_data``/``calc_strain`` 处理变形数据。

.. note::

   这些函数依赖可选的 ``myvasp`` 包和 VASP 输出目录结构。文档构建时
   ``myvasp`` 被 mock。

Related tutorials
   :doc:`../tutorials/cij_energy_fitting`

Convergence post-processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: mymetal.post.convergence
   :members:
   :show-inheritance:

用途
   后处理 ENCUT/k-point 收敛测试。``post_convergence`` 扫描
   ``y_convergence_encuts``/``y_convergence_kpoints`` 目录，提取能量，
   排序后绘制收敛曲线并写摘要文件。``my_write_convergence``/
   ``my_read_convergence`` 读写结果。

Stretch post-processing
~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: mymetal.post.stretch
   :members:
   :show-inheritance:

用途
   后处理 VASP/LAMMPS 拉伸计算。``post_stretch``/``post_lammps_stretch``
   从 ``y_dir/<strain>`` 提取能量-应变数据，``my_write_stretch``/
   ``my_read_stretch`` 读写结果。

GSFE / NEB post-processing
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::

   ``mymetal.post.gsfe`` 和 ``mymetal.post.neb`` 提供 GSFE 和 NEB 后处理
   （曲线绘制、spline 拟合、摘要文件）。它们的 docstring 当前包含 RST
   格式问题，暂不自动渲染完整 API；函数签名和用法请直接查看源码或
   :doc:`../manual/workflows` 中的 Advanced Workflow 概览。

.. automodule:: mymetal.post.kpar_ncore
   :members:
