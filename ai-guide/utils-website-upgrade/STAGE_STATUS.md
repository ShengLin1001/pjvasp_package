# 阶段状态追踪 — utils-website-upgrade

> 此文件是多阶段文档升级的唯一状态源。cron job 每次运行时先读此文件，
> 判断当前阶段，执行完毕后更新状态。

## 当前状态

- **工作分支**: `docs/utils-website-upgrade`
- **当前阶段**: Stage 2 pending（等待 cron 执行）
- **上一阶段完成时间**: Stage 1 — 2026-08-03
- **最后更新**: 2026-08-03 stage-1-complete

## 阶段进度

| 阶段 | 状态 | 开始时间 | 完成时间 | commit | 说明 |
|------|------|----------|----------|--------|------|
| Stage 0 | completed | 2026-08-03 | 2026-08-03 | — | 分支创建 + 计划 + 基线构建 + cron |
| Stage 1 | completed | 2026-08-03 | 2026-08-03 | — | vasp_utils/vasp_universal (22个脚本) |
| Stage 2 | pending | — | — | — | vasp_utils/vasp_workflow_bulk + neb_utils |
| Stage 3 | pending | — | — | — | slurm_utils |
| Stage 4 | pending | — | — | — | lmp_utils |
| Stage 5 | pending | — | — | — | n2p2_utils |
| Stage 6 | pending | — | — | — | 审查收尾 + 部署 |

## 阶段定义

### Stage 0 — 分支创建 + 计划文档 + 基线构建 + cron 设置
- [x] 从 main 创建分支 `docs/utils-website-upgrade`
- [x] 创建 `ai-guide/utils-website-upgrade/MASTER_PLAN.md`
- [x] 基线 Sphinx 构建成功（0 warnings）
- [x] 创建 STAGE_STATUS.md
- [x] 设置 cron job

### Stage 1 — vasp_utils/vasp_universal 子包
- 覆盖脚本：pei_vasp_univ_sbatch, pei_vasp_univ_post, pei_vasp_univ_clean_*,
  pei_vasp_univ_cp_contcar_cartesian_poscar, pei_vasp_univ_extract_*,
  pei_vasp_univ_find_and_change, pei_vasp_univ_increase_nbands,
  pei_vasp_univ_load_env, pei_vasp_univ_monitor_*, pei_vasp_univ_resubmit_*,
  pei_vasp_univ_transfer_*, pei_vasp_univ_get_size_by_distance.py,
  pei_vasp_univ_check_phase_transition, pei_vasp_univ_get_struct_infos,
  pei_vasp_plot_convergence.py
- 产出：`ai-guide/utils-website-upgrade/handoff/stage_1_handoff.md`

### Stage 2 — vasp_utils/vasp_workflow_bulk + neb_utils
- 覆盖脚本：pei_vasp_run_eos.py, pei_vasp_run_stretch.py, pei_vasp_run_neb,
  pei_vasp_run_convergence, pei_vasp_run_cij_energy, pei_vasp_run_cohesive,
  pei_vasp_run_dos_band, pei_vasp_run_surface_energy, pei_vasp_run_hoec_energy,
  pei_vasp_run_kpar_ncore, pei_vasp_plot_*, neb_utils/*
- 产出：`ai-guide/utils-website-upgrade/handoff/stage_2_handoff.md`

### Stage 3 — slurm_utils 子包
- 覆盖脚本：slurm_universal/*, slurm_vasp/*, slurm_n2p2/*
- 产出：`ai-guide/utils-website-upgrade/handoff/stage_3_handoff.md`

### Stage 4 — lmp_utils 子包
- 覆盖文件：template/*, post/*, lmp_universal/*
- 产出：`ai-guide/utils-website-upgrade/handoff/stage_4_handoff.md`

### Stage 5 — n2p2_utils 子包
- 覆盖脚本：data_input.py, data_read.py, sfs_gen_basic_SF.py,
  active_sf_0_*, n2p2_universal/*
- 产出：`ai-guide/utils-website-upgrade/handoff/stage_5_handoff.md`

### Stage 6 — 审查收尾 + 部署
- 全站交叉链接检查
- 首页 index.rst 更新
- 全量 Sphinx 构建
- git commit + push
- 部署建议
- 产出：`ai-guide/utils-website-upgrade/handoff/stage_6_handoff.md`

## 更新规则

cron job 每次运行时：
1. 读取本文件确定当前 pending 阶段
2. 执行该阶段
3. 完成后更新本文件：标记阶段 completed，填入时间和 commit hash
4. 写 handoff 文档
5. commit 本文件和 handoff
