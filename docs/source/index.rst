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
      <getting_started/au111_slab>`，或一次对比 :doc:`Cu(100)/(110)/(111)
      <getting_started/fcc_surfaces>` 与 :doc:`HCP Mg 多表面
      <getting_started/hcp_surfaces>`。

   .. container:: task-card

      **Tutorials**

      :doc:`按真实 fixture 学习 surface energy 与 OUTCAR 汇总
      <tutorials/index>`，或用合成数据走通 :doc:`EOS 拟合
      <tutorials/eos_curve>` 与 :doc:`单轴/双轴应变 <tutorials/biaxial_stretch>`。

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

.. figure:: /_static/images/generated/fcc_surfaces.png
   :width: 960px
   :alt: Side and top views of Cu(100), Cu(110) and Cu(111) 12-layer slabs

   FCC 多表面教程同时构建 Cu(100)/(110)/(111) 三个 12-layer slab。上排为
   侧视图，下排为俯视图。三个 slab 共享同一 ``a_fcc = 3.61 Å``，但面内
   排列和层间距完全不同。

.. figure:: /_static/images/generated/eos_curve.png
   :width: 960px
   :alt: Murnaghan and Birch-Murnaghan EOS fits to synthetic Cu-like data

   EOS 教程用 Cu-like 合成数据演示 Murnaghan 与 Birch-Murnaghan 拟合。
   左图为拟合曲线，右图为残差。

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
