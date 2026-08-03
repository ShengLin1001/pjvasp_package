# 阶段状态追踪 — utils-website-upgrade

> 此文件是多阶段文档升级的唯一状态源。cron job 每次运行时先读此文件，
> 判断当前阶段，执行完毕后更新状态。

## 当前状态

- **工作分支**: `docs/utils-website-upgrade`
- **当前阶段**: 全部完成（Stage 0-6 completed）
- **最后更新**: 2026-08-03 stage-4-figure-supplement + 状态文件同步修正
- **远程同步**: HEAD = origin/docs/utils-website-upgrade = f6fed87（已 push）

> **状态同步说明**：此前版本的状态文件停留在 Stage 2，但 git log 显示
> Stage 3/4/5/6 均已由前序 cron 完成 commit。本次 cron 触发时发现该不一致，
> 经核实 git 历史（commit ffcbb12/bad9640/7657857）确认 Stage 4/5/6 产出
> 完整且已在 origin，故不重复执行（遵守"不覆盖已有产出"约束），仅同步状态
> 文件并对 Stage 4 lmp_utils.rst 做小幅 figure 补强（原页面缺 figure 引用，
> 复用已生成的 lmp_utils_overview.png）。

## 阶段进度

| 阶段 | 状态 | 开始时间 | 完成时间 | commit | 说明 |
|------|------|----------|----------|--------|------|
| Stage 0 | completed | 2026-08-03 | 2026-08-03 | 2809f95 | 分支创建 + 计划 + 基线构建 + cron |
| Stage 1 | completed | 2026-08-03 | 2026-08-03 | 7be199c | vasp_utils/vasp_universal (22个脚本) |
| Stage 2 | completed | 2026-08-03 | 2026-08-03 | 111c8d1 | vasp_utils/vasp_workflow_bulk + neb_utils |
| Stage 3 | completed | 2026-08-03 | 2026-08-03 | 99cdfd2 | slurm_utils (8个脚本) |
| Stage 4 | completed | 2026-08-03 | 2026-08-03 | ffcbb12 | lmp_utils (模板+后处理+runner，20个脚本) |
| Stage 5 | completed | 2026-08-03 | 2026-08-03 | bad9640 | n2p2_utils (数据/SF/训练/清理，15个脚本) |
| Stage 6 | completed | 2026-08-03 | 2026-08-03 | 7657857 | 首页图廊 + 全量构建验证 + 部署准备 |

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
- 产出：`ai-guide/utils-website-upgrade/handoff/stage_3_handoff.md`（注：前序 cron 未写此 handoff 文件，但 commit 99cdfd2 产出完整：slurm_utils_overview.py + slurm_utils.rst + slurm_script_generation.rst + PNG + toctree/CI 注册）

### Stage 4 — lmp_utils 子包
- 覆盖文件：template/*, post/*, lmp_universal/*
- 产出：`ai-guide/utils-website-upgrade/handoff/stage_4_handoff.md`（本 cron 补写）
- 补强：本 cron 在 lmp_utils.rst 补充 figure 引用（复用 lmp_utils_overview.png），原页面仅有 seealso 交叉链接缺结果图

### Stage 5 — n2p2_utils 子包
- 覆盖脚本：data_input.py, data_read.py, sfs_gen_basic_SF.py,
  active_sf_0_*, n2p2_universal/*
- 产出：`ai-guide/utils-website-upgrade/handoff/stage_5_handoff.md`（注：前序 cron 未写此 handoff 文件，但 commit bad9640 产出完整）

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
