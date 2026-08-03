neb_utils NEB 后处理工具集
============================

``vasp_utils/neb_utils`` 是 NEB（Nudged Elastic Band）计算的后处理工具集。
这些脚本从 VTST tools 衍生，并扩展了 Python 重写版本：原版 perl 脚本保留
VTST 兼容性，Python 版本用 ASE 读 ``CONTCAR`` 并输出更完整的数据。

.. note::

   这些脚本不注册 ``console_scripts``。请将 ``vasp_utils/neb_utils/`` 加入
   ``PATH``，或使用完整路径调用。首次使用应先 ``-h`` 或阅读源码头部注释。
   perl 脚本需要与 VTST 的 ``Vasp.pm`` 同目录（或通过 ``FindBin`` 定位）。

.. contents:: 本页内容
   :local:
   :depth: 2

数据流
------

NEB 后处理的数据流如下：

.. code-block:: text

   00/OUTCAR  01/OUTCAR  ...  NN/OUTCAR
        │         │              │
        └─────────┴──────────────┘
                  ▼
        nebbarrier / nebbarrier.py      ← 提取距离、力、能量
                  ▼
        neb.dat / p_neb_py.dat / p_neb_full.dat
                  ▼
        neb_plot.py                     ← 绘图
                  ▼
        p_post_neb_full_selected.pdf

- ``nebef`` 提取每 image 的 max force + energy（``nebef.dat``）
- ``nebbarrier`` 提取累积距离、力、能量（``neb.dat``）
- ``nebbarrier.py`` / ``nebbarrier_full.py`` 是 Python 重写版
- ``neb_plot.py`` 读 ``p_neb_full.dat`` 绘完整数据图
- ``neb_select.py`` 从 full 数据提取指定 frame
- ``nebmovie`` 生成 xyz movie
- ``nebresults`` 是一站式入口：解压 OUTCAR + 调 nebef/nebbarrier/nebmovie

vtst 原版 vs Python 重写
------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 20 25 30

   * - 脚本
     - 语言
     - 输出
     - 说明
   * - ``pei_vasp_univ_neb_nebbarrier``
     - perl (vtst)
     - ``neb.dat``
     - 原版，读 ``00..NN/OUTCAR``
   * - ``pei_vasp_univ_neb_nebbarrier.py``
     - python
     - ``p_neb_py.dat``
     - 用 ASE 读 ``CONTCAR``，累积距离/力/能量
   * - ``pei_vasp_univ_neb_nebbarrier_full.py``
     - python
     - ``p_neb_full.dat``
     - 完整版：每帧每 image 的力/能量/距离
   * - ``pei_vasp_univ_neb_nebef``
     - perl (vtst)
     - ``nebef.dat``
     - 每 image max force + energy
   * - ``pei_vasp_univ_neb_nebmovie``
     - bash
     - ``movie_POSCAR.xyz`` / ``movie_CONTCAR.xyz``
     - 用 ASE 合并 xyz
   * - ``pei_vasp_univ_neb_nebresults``
     - perl (vtst)
     - 调用上述三者
     - 解压 OUTCAR + 一站式后处理
   * - ``pei_vasp_univ_neb_plot.py``
     - python
     - ``p_post_neb_full_selected.pdf``
     - 绘完整 NEB 数据图，支持切片
   * - ``pei_vasp_univ_neb_select.py``
     - python
     - ``p_neb_select_*.dat``
     - 提取指定 frame 的 nebef 数据

pei_vasp_univ_neb_nebbarrier
----------------------------

原版 VTST perl 脚本，读 ``00..NN/OUTCAR``，输出 ``neb.dat``。关闭了原版的
zip 部分（本仓库不需要）。脚本需与 VTST 的 ``Vasp.pm`` 同目录。

.. code-block:: bash

   pei_vasp_univ_neb_nebbarrier            # 自动发现 00..NN 子目录
   pei_vasp_univ_neb_nebbarrier 01 02 03   # 指定 image 目录

自动检测 ``LNEBCELL`` （SSNEB）标志以区分固定晶格 NEB 与可变晶格 SSNEB。

pei_vasp_univ_neb_nebbarrier.py
-------------------------------

Python 重写版，用 ASE 读 ``CONTCAR``（而非 perl 版的 OUTCAR-only），输出
``p_neb_py.dat``。提取每个 image 的累积距离、力与能量。

.. code-block:: bash

   pei_vasp_univ_neb_nebbarrier.py

自动发现当前目录下 ``[0-9][0-9]`` 子目录并排序。输出列：image id、累积距离
（Å）、能量（eV）、力（eV/Å）。

pei_vasp_univ_neb_nebbarrier_full.py
------------------------------------

完整版 Python 脚本，输出 ``p_neb_full.dat``：**每帧每 image** 的力/能量/距离。
与 ``nebbarrier.py``（只取最后一帧）不同，``_full`` 保留所有 ionic 步，供
``neb_plot.py`` 做时序切片与动画分析。

.. code-block:: bash

   pei_vasp_univ_neb_nebbarrier_full.py

输出包含：image id、frame、累积距离（Å）、相对能量（eV）、绝对能量（eV）、
力（eV/Å）。

pei_vasp_univ_neb_nebef
-----------------------

原版 VTST perl 脚本，提取每 image 的 max force + total energy，输出
``nebef.dat``。扫描当前目录下 ``[0-9][0-9]/OUTCAR``。

.. code-block:: bash

   pei_vasp_univ_neb_nebef

输出列：image id、max force（eV/Å）、total energy（eV）、relative energy（eV）。

pei_vasp_univ_neb_nebmovie
--------------------------

bash 脚本，生成 NEB movie 的 xyz 文件（POSCAR 或 CONTCAR），用 ASE 合并。

.. code-block:: bash

   pei_vasp_univ_neb_nebmovie              # 默认用 POSCAR
   pei_vasp_univ_neb_nebmovie 1            # 用 CONTCAR

输出 ``movie_POSCAR.xyz`` 或 ``movie_CONTCAR.xyz``。收集所有数字子目录
（``00``、``01``、…）中的目标文件，用 ASE 的 ``read_vasp`` + ``write``
合并为多帧 xyz。

pei_vasp_univ_neb_nebresults
----------------------------

原版 VTST perl 脚本，一站式 NEB 后处理入口：

1. 解压所有 image 目录的 ``OUTCAR.gz`` / ``OUTCAR.bz2``
2. 调用 ``nebef``
3. 调用 ``nebbarrier``
4. 调用 ``nebmovie`` （用 CONTCAR）

.. code-block:: bash

   pei_vasp_univ_neb_nebresults

.. note::

   假设 ``vfin.pl`` 已运行过（OUTCAR 已压缩）。若 OUTCAR 未压缩也能工作。

pei_vasp_univ_neb_plot.py
-------------------------

Python 脚本，绘制完整 NEB 数据图，支持 ``-start``/``-end`` 切片选择特定
frames。读 ``p_neb_full.dat``（由 ``nebbarrier_full.py`` 生成）。

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - 参数
     - 说明
   * - ``-start``
     - 切片起始位置（包含），支持负数，默认 ``None`` （从头）
   * - ``-end``
     - 切片结束位置（不包含），支持负数，默认 ``None`` （到尾）

.. code-block:: bash

   pei_vasp_univ_neb_plot.py                       # 全部 frames
   pei_vasp_univ_neb_plot.py -start 0 -end 10      # 前 10 帧
   pei_vasp_univ_neb_plot.py -start -5             # 最后 5 帧

输出 ``p_post_neb_full_selected.pdf``。

pei_vasp_univ_neb_select.py
---------------------------

Python 脚本，从 ``p_neb_full.dat`` 提取指定 frame 的 nebef 数据，输出
``p_neb_select_*.dat``。

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - 参数
     - 说明
   * - ``-frame``
     - **必需**。要提取的 frame 编号

.. code-block:: bash

   pei_vasp_univ_neb_select.py -frame 5

生成 ``p_neb_select_neb_*.dat`` 与 ``p_neb_select_nebef_*.dat``，可用于单独
分析某个 ionic 步的 NEB 能垒形状。

典型流程
--------

.. code-block:: bash

   # 1. NEB 计算完成后，在 y_neb 同级（含 00..NN 目录）运行
   pei_vasp_univ_neb_nebresults              # 一站式：解压 + nebef + nebbarrier + nebmovie

   # 2. 或分步运行
   pei_vasp_univ_neb_nebef                   # → nebef.dat
   pei_vasp_univ_neb_nebbarrier              # → neb.dat（perl 原版）
   pei_vasp_univ_neb_nebbarrier.py           # → p_neb_py.dat（python 重写）
   pei_vasp_univ_neb_nebbarrier_full.py      # → p_neb_full.dat（完整版）

   # 3. 绘图
   pei_vasp_univ_neb_plot.py                 # → p_post_neb_full_selected.pdf

   # 4. 提取特定 frame
   pei_vasp_univ_neb_select.py -frame 5      # → p_neb_select_*.dat

安全检查清单
------------

1. perl 脚本需与 VTST ``Vasp.pm`` 同目录——确认 ``FindBin`` 能定位
2. ``nebresults`` 会解压 ``OUTCAR.gz``——确认磁盘空间
3. ``nebbarrier_full.py`` 输出可能很大（每帧每 image）——先小规模测试
4. ``neb_plot.py`` 的 ``-start/-end`` 支持负数——注意越界

.. seealso::

   - :doc:`vasp_workflow_bulk` — NEB 目录准备脚本 ``pei_vasp_run_neb``
   - :doc:`vasp_universal` — 单目录 runner/post 工具
