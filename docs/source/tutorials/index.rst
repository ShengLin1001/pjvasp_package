Tutorials
=========

每个教程都从仓库 fixture 或合成数据复算，且不需要运行 VASP。研究 notebook
不直接作为文档源；页面使用整理后的脚本、确定性输出和明确验证。

.. toctree::
   :maxdepth: 1

   surface_energy
   outcar_batch
   eos_curve
   biaxial_stretch
   kpoints_sampling
   schmid_factor
   neighbor_distances
   atom_manipulation
   strain_deformation
   gsfe_models
   cubic_cell_and_stretch
   deformation_and_hnf
   reciprocal_lattice
   io_extxyz_and_general
   cij_energy_fitting
   miller_index_and_density
   periodic_table_and_arkel
   slurm_script_generation
   plot_gallery
   plot_gallery_general
   plot_gallery_plot
   plot_gallery_workflow
   plot_gallery_energy
   plot_gallery_atominfo
   plot_gallery_n2p2
   plot_gallery_plotting
   plot_gallery_render
   plot_gallery_ppt
   plot_gallery_oldplotdos

后续案例
--------

strained films 的真实 VASP 复算、SLURM dry-run、n2p2 dataset 与
heterostructure 仍在规划中。GSFE、decohesion 和 HOEC 当前只在
:doc:`../manual/workflows` 提供 Advanced Workflow 概览，尚未伪造缺少
审核 fixture 的完整教程。
