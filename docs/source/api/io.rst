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

.. automodule:: mymetal.io.general

.. automodule:: mymetal.io.extxyz

.. automodule:: mymetal.io.post.construct
   :members:

.. automodule:: mymetal.io.post.write
   :members:
