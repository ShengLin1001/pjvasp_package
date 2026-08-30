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
   surface_passivation
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

Post-processing 后处理分析
--------------------------

以下教程基于真实集群计算数据 (zcm6 / LAMMPS) 或合成数据，
演示 ``mymetal.post`` 各模块的后处理分析流程。

.. toctree::
   :maxdepth: 1

   post_hoec_energy
   post_cij_comparison
   post_stretch_analysis
   post_gsfe_analysis
   post_convergence
   post_relax_convergence
   post_kpar_ncore
   post_e_in_1_2_bulk

后续案例
--------

strained films 的真实 VASP 复算、SLURM dry-run、n2p2 dataset 与
heterostructure 仍在规划中。
