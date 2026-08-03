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

.. figure:: /_static/images/generated/bulk_structures.png
   :width: 960px
   :alt: 3D and side views of FCC Cu, BCC Fe, HCP Mg and Diamond Si conventional cells

   :doc:`getting_started/bulk_structures` 一次性构建四种常见 bulk cell
   (FCC/BCC/HCP/Diamond)，对比 cell 参数与原子排列。

.. figure:: /_static/images/generated/eos_curve.png
   :width: 960px
   :alt: Murnaghan and Birch-Murnaghan EOS fits to synthetic Cu-like data

   EOS 教程用 Cu-like 合成数据演示 Murnaghan 与 Birch-Murnaghan 拟合。
   左图为拟合曲线，右图为残差。

.. figure:: /_static/images/generated/kpoints_sampling.png
   :width: 960px
   :alt: Monkhorst-Pack vs Gamma-centered k-points and RK scan for a Cu(111) slab

   :doc:`tutorials/kpoints_sampling` 同时展示 MP/Gamma k 点位置、
   ``RK`` 扫描与 ``n_x/n_z`` 比值，覆盖 VASP ``KSPACING`` 与 ``R_k``
   两种自动网格。

.. figure:: /_static/images/generated/schmid_factor.png
   :width: 960px
   :alt: FCC Schmid factor polar map and per-slip-system bar chart for [1, 1, 6] loading

   :doc:`tutorials/schmid_factor` 在 (001) 立方极图上扫描最大 Schmid
   因子分布，并列出 ``[1, 1, 6]`` 加载方向下 12 个 FCC 滑移系的 ``m``。

.. figure:: /_static/images/generated/neighbor_distances.png
   :width: 960px
   :alt: Normalized RDF and cumulative coordination number for FCC Cu, BCC Fe and HCP Mg

   :doc:`tutorials/neighbor_distances` 同时给出 FCC/BCC/HCP 三种结构的
   归一化 RDF 与累计配位数曲线。

.. figure:: /_static/images/generated/atom_manipulation.png
   :width: 960px
   :alt: Reference Au(111) slab, slab with bottom half frozen, and slab with top layer shifted to bridge site

   :doc:`tutorials/atom_manipulation` 演示 ``FixAtoms`` 约束、按位置
   范围筛选原子、写出 ``Selective dynamics`` POSCAR 的完整流程。

.. figure:: /_static/images/generated/strain_deformation.png
   :width: 960px
   :alt: Lagrangian strain tensors for uniaxial x, biaxial xy and simple shear xy

   :doc:`tutorials/strain_deformation` 把三种参考变形转换为 deformation
   gradient ``F`` 与 Lagrangian 应变 ``E``，并以热图标注每个分量。

.. figure:: /_static/images/generated/gsfe_models.png
   :width: 960px
   :alt: Side views of FCC(111), HCP(0001) basal and HCP(10-10) prism I GSFE supercells

   :doc:`tutorials/gsfe_models` 一次性构建 FCC(111)、HCP(0001) 基面和
   HCP(10-10) prism I 三种 GSFE 超胞。

.. figure:: /_static/images/generated/cubic_cell_and_stretch.png
   :width: 960px
   :alt: Hexagonal primitive cell, orthorhombic cell, and three uniaxially stretched cells

   :doc:`tutorials/cubic_cell_and_stretch` 把 primitive 六角 cell 转成正交
   cell（面积守恒），再沿 x 方向施加拉伸扫描。

.. figure:: /_static/images/generated/deformation_and_hnf.png
   :width: 960px
   :alt: Deformation gradient schematic and Hermite Normal Form heatmaps

   :doc:`tutorials/deformation_and_hnf` 计算变形梯度 ``F`` 和整数矩阵的
   Hermite Normal Form。

.. figure:: /_static/images/generated/reciprocal_lattice.png
   :width: 960px
   :alt: Reciprocal lattice vector lengths and RK-based k-point meshes for FCC/BCC/HCP

   :doc:`tutorials/reciprocal_lattice` 对比叉积法和矩阵求逆法计算倒格子
   矢量，并展示 RK → k 点网格映射。

.. figure:: /_static/images/generated/cij_energy_fitting.png
   :width: 960px
   :alt: Strain-energy curves and input vs fitted Cij bar chart for synthetic Cu

   :doc:`tutorials/cij_energy_fitting` 用合成 Cu-like 数据演示 energy-strain
   法拟合二阶弹性常数 C11/C12/C44。

.. figure:: /_static/images/generated/periodic_table_heatmap.png
   :width: 960px
   :alt: Periodic table heatmap of bulk moduli for selected pure elements

   :doc:`tutorials/periodic_table_and_arkel` 在周期表上按元素值着色，
   并把二元化合物绘制在 van Arkel 三角图上。

.. figure:: /_static/images/generated/miller_index_and_density.png
   :width: 960px
   :alt: HCP Miller index conversion table and density bar chart

   :doc:`tutorials/miller_index_and_density` 演示 HCP 方向 3↔4 指数转换
   和 FCC/BCC/HCP/diamond 密度计算。

.. figure:: /_static/images/generated/slurm_script_generation.png
   :width: 960px
   :alt: Slurm script generation pipeline flowchart and rendered script text

   :doc:`tutorials/slurm_script_generation` 以 dry-run 模式生成 Slurm
   作业脚本，不调用 ``sbatch``。

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
   manual/vasp_universal
   manual/vasp_workflow_bulk
   manual/neb_utils
   manual/slurm
   manual/slurm_utils
   manual/lammps
   manual/lmp_utils
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
