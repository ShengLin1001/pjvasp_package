# 阶段状态追踪

> 此文件是多阶段文档升级的唯一状态源。cron job 每次运行时先读此文件，
> 判断当前阶段，执行完毕后更新状态。

## 当前状态

- **工作分支**: `docs/mymetal-website-upgrade`
- **当前阶段**: Stage 4（待开始）
- **上一阶段完成时间**: Stage 3 — 2026-08-02
- **最后更新**: 2026-08-02 stage-3-complete

## 阶段进度

| 阶段 | 状态 | 开始时间 | 完成时间 | commit | 说明 |
|------|------|----------|----------|--------|------|
| Stage 0 | completed | 2026-08-02 | 2026-08-02 | 见下 | 主题/CSS 升级 + 基线构建 |
| Stage 1 | completed | 2026-08-02 | 2026-08-02 | 见下 | build 子包 |
| Stage 2 | completed | 2026-08-02 | 2026-08-02 | 见下 | calculate 子包 |
| Stage 3 | completed | 2026-08-02 | 2026-08-02 | 见下 | io + post 子包 |
| Stage 4 | pending | - | - | - | universal 子包 |
| Stage 5 | pending | - | - | - | ml + slurm + cr 子包 |
| Stage 6 | pending | - | - | - | 审查收尾 + 部署 |

## 阶段定义

### Stage 0 — 主题/CSS 升级 + 基线构建测试
- 更新 `docs/source/_static/css/custom.css` 实现 ASE 风格色板（见 ASE_AESTHETIC_STUDY.md 路径 A）
- 更新 `docs/source/conf.py`（如需加 extension）
- 运行基线 Sphinx 构建确认无回归
- 产出：`ai-guide/multistage-docs/handoff/stage_0_handoff.md`

### Stage 1 — build 子包
- 覆盖模块：bulk/create.py, bulk/gsfe.py, film/extrfilm.py, film/findcubic.py,
  film/findhetero.py, film/findprim.py, film/hydroxyl.py, film/stretch.py,
  workflow/hoec.py, workflow/kpar_ncore.py, workflow/general.py
- 重点函数：generate_film（已有）, create_fcc_111, create_hcp_basal,
  create_gsfe_model, cal_area（已有）, find_cubic, stretch_along_direction_to_cell,
  add_hydro_atoms, prepare_hoec_reference, generate_hoec_dirs
- 新增示例：bulk_structures（已有）, gsfe_model, hydroxyl_passivation,
  cubic_cell_finding, hoec_deformation_setup
- 产出：`ai-guide/multistage-docs/handoff/stage_1_handoff.md`

### Stage 2 — calculate 子包
- 覆盖模块：calenergy/surfenergy.py（已有）, calmath/matrix.py,
  calmechanics/deformation.py, calmechanics/strain.py（已有 tutorial）,
  calmechanics/stretch.py, calmechanics/hoec.py, calmismatch/calhetero.py,
  calqm/kpoints.py（已有 tutorial）, electronic_structure/plotter.py,
  electronic_structure/universal.py, material_science/schmid.py（已有 tutorial）
- 重点函数：cal_surface_energy（已有）, hermite_normal_form, cal_deform_matrix,
  cal_strain_matrix, cal_von_mises_strain, cal_mismatch, get_kpoints_by_size,
  cal_reciprocal_matrix, cal_fcc_schmid_factors（已有）, plot_brillouin_zone_from_kpath
- 新增示例：hermite_normal_form, mismatch_analysis, brillouin_zone_plotting,
  band_structure_summary
- 产出：`ai-guide/multistage-docs/handoff/stage_2_handoff.md`

### Stage 3 — io + post 子包
- io: vasp.py, extxyz.py, general.py, post/construct.py, post/write.py
- post: Cij_energy.py, convergence.py, gsfe.py, hoec_energy.py, neb.py,
  newmain.py, stretch.py, relax_convergence.py, E_in_1_2_bulk.py, kpar_ncore.py
- 重点：my_read_vasp, my_write_vasp（已有）, extxyz_to_atomlist, general_read,
  PostData/PostTime（已有 API）, post_gsfe, post_neb, post_stretch,
  post_convergence, fit_cij_energy
- 新增示例：vasp_io_roundtrip, extxyz_trajectory, outcar_batch（已有）,
  convergence_post, cij_energy_post
- 产出：`ai-guide/multistage-docs/handoff/stage_3_handoff.md`

### Stage 4 — universal 子包
- atom: delatom, density, fixatom, miller, moveatom, neighbor（已有 tutorial）, rotate
- check: atoms, checkinput, convergence, type
- data: dataadjust, datachange, patterntrans
- math: operations
- matrix: adjust
- plot: atominfo, energy, general, n2p2, plot, plotting, render, workflow
- 重点：fixatoms（已有 tutorial）, move_atoms, cal_density, get_neighbor_distances,
  three_index_to_four_index, check_phase_transition, my_plot_convergence,
  my_plot_cij_energy, my_plot_gsfe, periodic_table_heatmap, van_arkel_triangle
- 新增示例：miller_index_conversion, cna_analysis, density_calculation,
  universal_plotting_showcase, periodic_table_heatmap
- 产出：`ai-guide/multistage-docs/handoff/stage_4_handoff.md`

### Stage 5 — ml + slurm + cr 子包
- ml: dataset, model, confusionmatrix, plot, n2p2/dataset, n2p2/workflow,
  n2p2/calculate/cur, n2p2/calculate/sf, n2p2/calculate/post
- slurm: submit.py（pei_slurm_univ_submit）
- cr: crsum3, crplotkcontactgraph
- 重点：PeiN2p2 workflow, cur_select, generate_radial_blocks, train_model,
  generate_slurm_script_base, K_handle/Freq_handle
- 新增示例：n2p2_sf_generation, cur_selection_demo, slurm_script_generation,
  ml_model_training（合成数据）, cr_stiffness_fit（合成数据）
- 产出：`ai-guide/multistage-docs/handoff/stage_5_handoff.md`

### Stage 6 — 审查收尾 + 部署
- 全站交叉链接检查
- 首页 index.rst 更新（加新图廊、新 tutorial 入口）
- 移动端 390px 检查
- linkcheck
- 全量 Sphinx 构建 -W
- compileall mymetal
- git commit（分功能边界）
- git push origin docs/mymetal-website-upgrade
- 检查 GitHub Actions 构建
- 产出：`ai-guide/multistage-docs/handoff/stage_6_handoff.md`

## 更新规则

cron job 每次运行时：
1. 读取本文件确定当前 pending 阶段
2. 执行该阶段
3. 完成后更新本文件：标记阶段 completed，填入时间和 commit hash
4. 写 handoff 文档
5. commit 本文件和 handoff
