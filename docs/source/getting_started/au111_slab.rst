.. _getting-started-au111:

生成并验证 Au(111) slab
========================

:Audience: 第一次使用 ``mymetal`` 的用户
:Time: 10--15 分钟
:Requires: Python 3.10、pjvasp_package、ASE、Matplotlib
:Runs VASP: No

目标
----

使用真实入口 :func:`mymetal.build.film.stretch.generate_film` 构建
12-layer FCC Au(111) slab，打印可检查信息，写出并重新读取 ``POSCAR``，
最后生成非空结构图。

最终结果
--------

.. figure:: /_static/images/generated/au111_slab.png
   :width: 960px
   :alt: Side and top views of a twelve-layer Au(111) slab with vacuum along z

   左图把同一个 primitive slab 沿面内重复六次以便观察，并保留上下各
   20 Å 真空；右图显示 Au(111) 面内排列。源结构来自本页引用的唯一脚本，
   图片由 ``docs/scripts/generate_structure_images.py`` 复用该脚本生成。

前置条件
--------

先完成 :doc:`../user_guide/install`，并确认：

.. code-block:: console

   $ python -c "import ase, matplotlib, mymetal; print('imports: ok')"
   imports: ok

.. note::

   本教程不运行 VASP，不调用 ``sbatch``，也不需要或生成 POTCAR。

输入和起始目录
--------------

从仓库根目录运行。输入只有一个可下载 Python 脚本：

.. code-block:: text

   pjvasp_package/
   └── docs/
       └── examples/
           └── getting_started_au111.py

运行完整示例
------------

:download:`下载 getting_started_au111.py
<../../examples/getting_started_au111.py>`，或在仓库根目录直接运行：

.. code-block:: console

   $ python docs/examples/getting_started_au111.py --output docs/_build/example-au111

页面与用户下载复用同一份实现：

.. literalinclude:: ../../examples/getting_started_au111.py
   :language: python
   :linenos:

预期终端输出
------------

``generate_film`` 会先打印层数诊断，随后脚本给出确定性摘要：

.. code-block:: text

   number_per_layer: 1.0
   atom_number_per_slab: 1
   layer_number_per_slab: 1.0
   num_rep_z: 12
   len(my_slab): 12
   formula: Au12
   atoms: 12
   cell (angstrom):
     [2.950000, 0.000000, 0.000000]
     [-1.475000, 2.554775, 0.000000]
     [0.000000, 0.000000, 66.495314]
   pbc: [True, True, True]
   round-trip: ok
   image-check: nonblank
   wrote: .../POSCAR
   wrote: .../au111_slab.png

输出目录为：

.. code-block:: text

   docs/_build/example-au111/
   ├── POSCAR
   └── au111_slab.png

如何解释结果
------------

* ``Au12`` 与 ``atoms: 12`` 验证一个 primitive 面内单元包含 12 层、每层
  一个原子。
* c 方向 cell 长度约 66.495 Å；原子 slab 厚度约 26.495 Å，ASE 在两侧各
  加入约 20 Å 真空。
* PBC 在三个方向均为 ``True``。z 方向仍周期重复，真空用于隔开相邻周期像。
* ``round-trip: ok`` 表示项目的 VASP writer/reader 保留 formula、cell、
  positions 与 PBC。
* ``image-check: nonblank`` 会检查 PNG 存在、大小合理且不是近似纯白图。

验证方法
--------

脚本在写出结果前执行以下断言：

* ``len(slab) == 12`` 且 formula 为 ``Au12``；
* PBC 为 ``[True, True, True]``；
* 重新读取的 cell 和 positions 在 ``1e-10`` 容差内一致；
* PNG 文件存在且包含足够的非白像素。

任一检查失败时进程会以非零状态退出，因此同一脚本可以直接作为 CI smoke test。

常见错误
--------

``ModuleNotFoundError: mymetal``
   回到仓库根目录执行 ``python -m pip install -e .``，并确认当前 shell 使用
   同一个 Python。

``ModuleNotFoundError: adjustText`` 或 ``brokenaxes``
   图片使用仓库统一 Matplotlib 样式入口；安装 ``docs/requirements.txt``，
   或安装这两个直接依赖。

看不到真空
   顶视图沿 z 方向观察，不会显示真空；查看图片左侧 side view，或比较 atom
   z-span 与 cell c-length。

输出目录已经存在
   脚本只覆盖明确的 ``POSCAR`` 和 ``au111_slab.png``，不会递归清理目录。
   如需保留旧结果，请换一个 ``--output``。

下一步
------

* :doc:`../tutorials/surface_energy`
* :doc:`../tutorials/outcar_batch`
* :doc:`../manual/vasp`

Related API
-----------

* :func:`mymetal.build.film.stretch.generate_film`
* :func:`mymetal.io.vasp.my_read_vasp`
* :func:`mymetal.io.vasp.my_write_vasp`
