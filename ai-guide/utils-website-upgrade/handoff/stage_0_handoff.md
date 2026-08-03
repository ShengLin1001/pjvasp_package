# Stage 0 Handoff — 分支创建 + 计划文档 + 基线构建 + cron 设置

## 完成时间
2026-08-03

## 完成内容

1. **分支创建**：从 main 创建 `docs/utils-website-upgrade` 分支
   - main 分支保持不动
   - 当前所有工作在此分支进行

2. **计划文档**：创建 `ai-guide/utils-website-upgrade/MASTER_PLAN.md`
   - 定义 7 个阶段（Stage 0-6）
   - 覆盖 vasp_utils/vasp_universal, vasp_utils/vasp_workflow_bulk+neb_utils,
     slurm_utils, lmp_utils, n2p2_utils
   - 每阶段通用产出规范
   - 子agent分工模式
   - ASE 审美规范（沿用上一轮）

3. **基线构建**：Sphinx 构建成功，0 warnings
   - 命令：`python -m sphinx -E -b html --keep-going docs/source docs/_build/html`
   - 构建产出在 `docs/_build/html/`

4. **状态追踪**：创建 `STAGE_STATUS.md`
   - Stage 0 标记 completed
   - Stage 1-6 标记 pending

5. **Python 环境**：
   - Sphinx 构建 + 示例运行：`C:/Users/louis/mysoft/env/pyenv/codex/Scripts/python.exe` (3.10.7)
   - 已安装 sphinx 8.1.3, ase 3.28, mymetal, myvasp

## 仓库当前状态

- 分支：`docs/utils-website-upgrade`
- 基线构建：成功
- 现有文档：mymetal package 文档已完成（上一轮 Stage 0-6）
- 待做：vasp_utils, slurm_utils, n2p2_utils, lmp_utils 四个脚本包的文档

## 安全约束（每个 cron session 必须遵守）

1. 禁止递归删除命令
2. 不得修改 setup.py
3. 不得运行真实 VASP/LAMMPS/n2p2/sbatch
4. 不得发布 POTCAR/secrets
5. 示例脚本用 dry-run/合成数据
6. 不覆盖已有产出
7. 每阶段结束运行 Sphinx 构建
8. 中文叙事，英文 API/命令名

## 给 Stage 1 的指示

Stage 1 覆盖 `vasp_utils/vasp_universal/` 子包。这是最大的子包，包含 20+ 个脚本。
建议用 subagent 并行处理：
- Agent A：读脚本源码，写示例脚本，生成图片
- Agent B：写 RST 文档页，更新 manual/vasp.rst
- Agent C：审查

图片类型建议：流程图、目录树图、脚本内容展示、dry-run 输出展示。
