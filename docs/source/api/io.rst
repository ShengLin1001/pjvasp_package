Input and Output
================

.. automodule:: mymetal.io

.. automodule:: mymetal.io.vasp

``my_read_vasp``
----------------

.. autofunction:: mymetal.io.vasp.my_read_vasp

参数
   ``filename`` 是 POSCAR/CONTCAR 路径或由 ASE ``@reader`` 接受的 file-like
   对象；默认 ``"CONTCAR"``。

返回
   ``(atoms, lattice_scale_factor)``。``atoms`` 是 ``ase.Atoms``，
   cell/position 单位为 Å；第二项是 VASP 文件第二行的原始 dimensionless
   scale factor。源码上的旧 return annotation 只写 ``Atoms``，实际合同以此
   tuple 为准。

异常
   路径不存在时抛 ``FileNotFoundError``；无法从 VASP 4 格式 header 或相邻
   POTCAR/OUTCAR 推断元素时抛 ``ase.io.ParseError``；格式错误产生相应
   ``ValueError``/``IndexError``。

``my_write_vasp``
-----------------

.. autofunction:: mymetal.io.vasp.my_write_vasp

参数
   ``filename`` 为输出路径；``atoms`` 是一个 ``ase.Atoms`` 或只含一帧的
   list/tuple；``lattice_scale_factor`` 为 dimensionless VASP scale；
   ``direct`` 选择 Direct/Cartesian；``sort``、``symbol_count``、
   ``ignore_constraints`` 与 ``wrap`` 控制写出细节。

返回与副作用
   返回 ``None``，创建或覆盖明确的 VASP structure file。cell/positions
   以 Å 表示，除以 ``lattice_scale_factor`` 后按 VASP 约定写出。

异常
   多于一帧、退化 cell，或 VASP 不能表示的 ``FixedPlane``/``FixedLine``
   constraint 抛 ``RuntimeError``。

最小 round trip
----------------

.. code-block:: python

   from mymetal.io.vasp import my_read_vasp, my_write_vasp

   my_write_vasp("POSCAR", slab, label="Au(111) slab")
   slab_read, scale = my_read_vasp("POSCAR")
   assert scale == 1.0
   assert slab_read.get_chemical_formula() == slab.get_chemical_formula()

Related tutorials
-----------------

See :doc:`../getting_started/au111_slab`.

其他 I/O 模块
-------------

``extxyz_to_atomlist``
~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.io.extxyz.extxyz_to_atomlist

用途
   把 Extended XYZ 轨迹文件读成 ``ase.Atoms`` 列表，每帧执行 ``wrap()``。
   底层用 ``ase.io.read(file, index=':', format='extxyz')``。

参数/返回
   ``file`` 为文件路径；返回 ``list[Atoms]``，每帧一个 ``Atoms``。

最小示例
   .. code-block:: python

      from mymetal.io.extxyz import extxyz_to_atomlist

      atomlist = extxyz_to_atomlist('trajectory.xyz')
      print(len(atomlist), atomlist[0].get_chemical_formula())

Related tutorials
   :doc:`../tutorials/io_extxyz_and_general`

``general_read`` / ``general_write``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.io.general.general_read

.. autofunction:: mymetal.io.general.general_write

用途
   ``general_read`` 把分隔符文件读成 ``pandas.DataFrame``，支持自定义
   header、分隔符、注释符、index 列。``general_write`` 把 DataFrame 按
   列类型（int/float/bool/str）格式化写出为纯文本。

参数
   ``general_read``: ``filepath``、``has_header``、``header_names``、
   ``sep``（默认 ``r"\s+"``）、``comment_char``（默认 ``#``）、
   ``index_col``、``header_row``。
   ``general_write``: ``filename``、``dfc``、``int_format``/``float_format``/
   ``str_format``/``bool_format``、``if_write_col_num``、``if_write_row_num``。

最小示例
   .. code-block:: python

      import pandas as pd
      from mymetal.io.general import general_read, general_write

      df = pd.DataFrame({'encut': [300, 400, 500], 'energy': [-8.1, -8.2, -8.21]})
      general_write('data.txt', dfc=df, if_write_col_num=True)
      df_read = general_read('data.txt', has_header=True)

Related tutorials
   :doc:`../tutorials/io_extxyz_and_general`

.. automodule:: mymetal.io.post.construct
   :members:

.. automodule:: mymetal.io.post.write
   :members:
