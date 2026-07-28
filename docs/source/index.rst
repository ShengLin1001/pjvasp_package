pjvasp_package documentation
=============================

``pjvasp_package`` 是面向计算材料学的研究工具仓库；``mymetal`` 是其中可导入的
核心 Python package。仓库还提供 VASP、SLURM、LAMMPS 与 n2p2 的脚本和目录型
workflow。叙事内容使用中文，API、模块、命令与文件名保留英文。

它适合需要用 ASE/Python 构建结构、在 HPC 上组织 VASP 任务、从 ``OUTCAR``
提取结果，或准备机器学习势数据的研究用户。新用户可以先完成一个不需要 VASP、
POTCAR、SLURM 或集群账号的 Au(111) 案例。

可以完成什么
------------

* 构建 bulk、film、surface 与 heterostructure，并读写 POSCAR/CONTCAR；
* 组织 VASP/LAMMPS 输入、批量目录与 SLURM dry-run；
* 汇总能量、压力、体积、力和收敛状态；
* 准备 n2p2 数据，并查阅稳定 Python/脚本接口。

从这里开始
----------

.. container:: task-grid

   .. container:: task-card

      **Getting Started**

      :doc:`生成、写出并验证 12-layer Au(111) slab
      <getting_started/au111_slab>`。

   .. container:: task-card

      **Tutorials**

      :doc:`按真实 fixture 学习 surface energy 与 OUTCAR 汇总
      <tutorials/index>`。

   .. container:: task-card

      **Workflow Guides**

      :doc:`理解 VASP、SLURM、LAMMPS 与 n2p2 的协作边界
      <manual/workflows>`。

   .. container:: task-card

      **API / Script Reference**

      查 :doc:`Python API <api>` 或 :doc:`当前 pei_* 脚本
      <reference/scripts>`。

结果预览
--------

.. figure:: /_static/images/generated/au111_slab.png
   :width: 960px
   :alt: Side and top views of a twelve-layer Au(111) slab with vacuum along z

   Getting Started 生成的真实 Au(111) 结构。左图以六个面内重复单元显示
   12 个原子层和 z 方向真空，右图显示 (111) 面内排列。图片由
   ``docs/scripts/generate_structure_images.py`` 从教程脚本确定性生成。

.. toctree::
   :maxdepth: 2
   :caption: Start

   user_guide/install
   getting_started/index

.. toctree::
   :maxdepth: 2
   :caption: Learn by task

   tutorials/index
   manual/workflows
   manual/vasp
   manual/slurm
   manual/lammps
   manual/n2p2
   user_guide/examples
   user_guide/troubleshooting

.. toctree::
   :maxdepth: 2
   :caption: Reference and development

   api
   reference/scripts
   reference/dependencies
   reference/development
