.. _tutorial-io-extxyz-and-general:

Extended XYZ 轨迹与通用文本 I/O
================================

:Audience: 需要读写轨迹文件或格式化文本数据的用户
:Time: 6 分钟
:Requires: Python、ASE、pandas、pjvasp_package
:Runs VASP: No

目标
----

演示两个 I/O 工具：

1. :func:`mymetal.io.extxyz.extxyz_to_atomlist` — 把 Extended XYZ 轨迹文件
   读成 ``ase.Atoms`` 列表；
2. :func:`mymetal.io.general.general_read` /
   :func:`mymetal.io.general.general_write` — 格式化文本 DataFrame 读写。

运行完整脚本
------------

.. code-block:: console

   $ python docs/examples/io_extxyz_and_general.py --output docs/_build/example-io

.. literalinclude:: ../../examples/io_extxyz_and_general.py
   :language: python
   :linenos:

预期输出
--------

.. code-block:: text

   ================ Part A: extxyz trajectory I/O ================
   extxyz file written : .../trajectory.xyz
   frames read back : 5
   frame 0 formula : Cu4
   first-atom trajectory (Å):
     frame 0: [0.000 0.000 0.000]
     frame 1: [0.010 0.000 0.000]
     ...

   ================ Part B: general_read / general_write ================
   original DataFrame:
      encut   energy  pressure
   0  300.0 -15.6320      0.12
   ...
   round-trip DataFrame:
      encut   energy  pressure
   0  300.0 -15.6320      0.12
   ...
   numeric round-trip match per column:
     encut     : OK
     energy    : OK
     pressure  : OK

结果图
------

.. figure:: /_static/images/generated/io_extxyz_and_general.png
   :alt: Extended XYZ trajectory displacement and convergence data round-trip

   左图：5 帧 Cu 轨迹中第一个原子的位移随帧变化。右图：general_write 写出
   再 general_read 读回的 ENCUT 收敛数据（能量 vs ENCUT）。

验证方法
--------

* ``len(atomlist) == 5``，所有帧 formula 为 ``Cu4``；
* general_write → general_read 的数值列 ``allclose``；
* 图片非空白。

相关 API
--------

* :func:`mymetal.io.extxyz.extxyz_to_atomlist`
* :func:`mymetal.io.general.general_read`
* :func:`mymetal.io.general.general_write`
* :func:`mymetal.io.vasp.my_read_vasp`、 :func:`mymetal.io.vasp.my_write_vasp`
  （POSCAR/CONTCAR 读写，见 :doc:`../getting_started/au111_slab`）
