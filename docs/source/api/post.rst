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
   :doc:`../tutorials/post_cij_comparison`

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

Related tutorials
   :doc:`../tutorials/post_convergence`

Stretch post-processing
~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: mymetal.post.stretch
   :members:
   :show-inheritance:

用途
   后处理 VASP/LAMMPS 拉伸计算。``post_stretch``/``post_lammps_stretch``
   从 ``y_dir/<strain>`` 提取能量-应变数据，``my_write_stretch``/
   ``my_read_stretch`` 读写结果。

Related tutorials
   :doc:`../tutorials/post_stretch_analysis`

GSFE / NEB post-processing
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::

   ``mymetal.post.gsfe`` 和 ``mymetal.post.neb`` 提供 GSFE 和 NEB 后处理
   （曲线绘制、spline 拟合、摘要文件）。它们的 docstring 当前包含 RST
   格式问题，暂不自动渲染完整 API；函数签名和用法请直接查看源码或
   :doc:`../manual/workflows` 中的 Advanced Workflow 概览。

.. automodule:: mymetal.post.gsfe
   :members: post_gsfe, check_constraints, find_sf_usf, write_output
   :show-inheritance:

.. automodule:: mymetal.post.neb
   :members: post_neb, my_copy_neb_files, my_write_neb, my_add_head, my_read_neb, my_spline_neb
   :show-inheritance:

Related tutorials
   :doc:`../tutorials/post_gsfe_analysis`

.. automodule:: mymetal.post.kpar_ncore
   :members:

HOEC energy post-processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: mymetal.post.hoec_energy
   :members: post_hoec_energy, run_univ_post, load_hoec_manifest, collect_mode_data, fit_P, solve_constants, write_hoec_energy, read_hoec_energy
   :show-inheritance:

用途
   用 energy-strain 方法拟合 2/3/4 阶弹性常数 (SOEC/TOEC/FOEC)。
   从 ``y_hoec_energy/y_dir/<mode>`` 目录读取变形能量，对每个模式
   做 4 阶多项式拟合，提取 P2/P3/P4 系数，然后通过线性系统求解
   各阶独立弹性常数。支持固定 SOEC 拟合高阶项、窗口扫描和
   交叉验证。

Principal entry point: ``mymetal.post.hoec_energy.post_hoec_energy``。

.. note::

   这些函数依赖可选的 ``myvasp`` 包和 ``mymetal.calculate.calmechanics.hoec``
   模型。文档构建时 ``myvasp`` 被 mock。

Related tutorials
   :doc:`../tutorials/post_hoec_energy`

Relax convergence post-processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: mymetal.post.relax_convergence
   :members: my_univ_post_convergence, read_univ_post_convergence
   :show-inheritance:

用途
   后处理 VASP ionic relaxation 的每步能量和力收敛轨迹。读取
   ``pei_vasp_univ_extract_convergence`` 生成的数据文件，绘制
   对数纵轴的能量/力收敛曲线。

Related tutorials
   :doc:`../tutorials/post_relax_convergence`

General post-processing helpers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: mymetal.post.general
   :members:
   :show-inheritance:

用途
   ``my_ployfit`` 对数据做多项式拟合并返回系数；
   ``my_read_y_dir_contcar`` 从 ``y_dir`` 子目录读取 CONTCAR；
   ``my_sort`` 按自然排序排列字符串列表（``"10"`` 排在 ``"2"`` 之后）；
   ``get_structure_info`` 提取结构的 cell 和原子数信息。

.. note::

   这些函数不依赖 VASP 可执行文件，但 ``my_read_y_dir_contcar`` 需要磁盘
   上存在 CONTCAR 文件。

E_in_1_2_bulk post-processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: mymetal.post.E_in_1_2_bulk
   :members: post_E_in_1_2_bulk, my_write_E_in_1_2_bulk
   :show-inheritance:

用途
   后处理 E_in_1/2 双轴变形体相计算。从 ``y_dir/<a1-a2>`` 提取
   能量，生成 2D 等高线图和 profile 图，提取平衡参数。

Related tutorials
   :doc:`../tutorials/post_e_in_1_2_bulk`
