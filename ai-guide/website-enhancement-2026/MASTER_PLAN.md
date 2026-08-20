# 网站增强编排计划 (2026-08-20)

## 目标

三大任务，按优先级排列：

1. **优化现有图表** — 网页部署中图片字体过大、超出边界等问题
2. **新增 plot 函数文档页** — mymetal/universal/plot 下的函数图文并茂介绍
3. **扩充 workflow 文档** — vasp/lammps/n2p2 workflow 页面从简短目的扩充为图文并茂的详细介绍

## 执行约束

- 当前 session 只负责编排（委托 + 审阅 + 修复），不做具体工作
- 每个工作 subagent 配独立的审阅 subagent
- 重大问题→重新调起工作 subagent；小问题→in-session 修复
- 删除文件慎重，不批量删除
- 服务器文件只读；如需写，复制到 /public3/home/scg6928/mywork/temp/ 子目录
- **服务器任务一律在 zcm6 登录节点直接运行，绝不提交 Slurm 作业（sbatch/srun）**
- 工作文档放 ai-guide/website-enhancement-2026/ 子目录

## 环境

- Repo: F:\BaiduSyncdisk\version20240608\main_code_space\pjvasp_package
- Python (codex venv): C:/Users/louis/mysoft/env/pyenv/codex/Scripts/python.exe
- Sphinx build: <py> -m sphinx -E -b html -W --keep-going docs/source docs/_build/html
- SSH: ssh zcm6 (scg6928@ZC-M6)
- 当前分支: main (clean)
- utils-website-upgrade 分支已有 overview 图和 manual 页面（未合并到 main）

## 任务分解

### Task A: 优化现有图表 (worker + auditor)

**问题**: 网页中图片字体大小不合适、超出边界。首页 17 张 figure 都用 `:width: 960px`，
在 RTD 主题 (max-width: 1120px) 下可能过大。各 tutorial 页面的 figure 也有类似问题。

**工作内容**:
- 检查 docs/source/_static/css/custom.css，增加 figure 响应式样式
- 调整 index.rst 中 figure 的 width 参数（从固定 960px 改为响应式或更合理的值）
- 检查各 tutorial RST 中的 figure directive，确保 width 合理
- 可能需要重新生成部分 overview 图（调整 fig 字体大小/dpi）
- utils-website-upgrade 分支上的 overview 图需要合并到工作分支

**审阅重点**: 确认图片不溢出、字体可读、移动端友好

### Task B: plot 函数文档页 (worker + auditor)

**源代码**: mymetal/universal/plot/ 下 11 个模块
- general.py (general_set_all_rcParams, my_plot 入口, general_modify_legend 等)
- plot.py (my_plot, my_plot_colorbar, my_plot_brokenaxed 等)
- workflow.py (my_plot_convergence, my_plot_cij_energy, my_plot_neb 等)
- energy.py (my_plot_energy_components)
- atominfo.py (my_plot_interlayer_distance, my_plot_rdf 等)
- n2p2.py (my_plot_learning_curve, my_plot_rmse_by_tag 等)
- plotting.py (periodic_table_heatmap, van_arkel_triangle 等)
- render.py (my_render)
- ppt.py (PPT 导出工具)
- oldplotdos.py (DOS 绘制)

**工作内容**:
- 新建 docs/source/tutorials/plot_gallery.rst (或 plot_reference.rst)
- 为每个 plot 模块写一个图文并茂的介绍页
- 每个函数配最小示例代码 + 生成 PNG 示例图
- 注册到 toctree 和 generate_structure_images.py
- 遵循 p-plot-figure skill 的风格约束

**审阅重点**: 函数签名准确、示例可运行、图非空白、Sphinx build 0 warning

### Task C: 扩充 workflow 文档 (worker + auditor)

**当前状态**: vasp.rst (65行)、lammps.rst (13行)、n2p2.rst (20行) 都非常简短

**服务器数据源**:
- VASP: /public3/home/scg6928/mywork/20250521_au_n2p2/construct_dataset/calculate/
- n2p2: /public3/home/scg6928/mywork/20250521_au_n2p2/train/
- LAMMPS: /public3/home/scg6928/mysoft/lammps/lammps/lammps-nc/pj-test-properties-gold/

**工作内容**:
- SSH 到 zcm6 收集 workflow 实际运行结果（目录结构、输入输出文件、结果数据）
- 在登录节点直接运行 python 脚本生成结果图表（不提交 Slurm）
- 复制必要文件到 /public3/home/scg6928/mywork/temp/ 如需写操作
- 为每个 workflow 写图文并茂的详细介绍：
  - 应用场景
  - 输入文件结构
  - 执行流程（配流程图）
  - 实际结果展示（配真实数据图表）
  - 参数说明
  - 与 mymetal 库函数的关联

**审阅重点**: 数据真实（来自服务器）、图非空白、不泄露敏感信息、Sphinx build 0 warning

## 编排流程

1. 创建工作分支 docs/website-enhancement
2. 先合并 utils-website-upgrade 分支的 overview 图和 manual 页面到工作分支
3. 并行委派 Task A (优化图表) + Task B (plot 文档) + Task C (workflow 文档)
   - 每个 task 一个 worker subagent
   - 每个 worker 完成后配一个 auditor subagent
   - 重大问题→重新调起 worker；小问题→in-session 修复
4. 各 task 完成后 commit
5. 最终统一 build 验证 + 链接检查

## 状态跟踪

| Task | Worker | Auditor | Status | Notes |
|------|--------|---------|--------|-------|
| A: 优化图表 | pending | pending | not_started | |
| B: plot 文档 | pending | pending | not_started | |
| C: workflow 文档 | pending | pending | not_started | |
