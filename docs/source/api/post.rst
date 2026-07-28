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

.. automodule:: mymetal.post.kpar_ncore
   :members:
