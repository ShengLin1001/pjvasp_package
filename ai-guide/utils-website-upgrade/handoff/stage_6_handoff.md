# Stage 6 Handoff — 审查收尾 + 部署

## 完成时间
2026-08-03

## 完成内容

### 首页图廊更新
在 index.rst 新增 5 张 utils 子包概览图，每张链接到对应 manual 页面：
- vasp_universal_overview.png → manual/vasp_universal
- vasp_workflow_bulk_overview.png → manual/vasp_workflow_bulk
- slurm_utils_overview.png → manual/slurm_utils
- lmp_utils_overview.png → manual/lmp_utils
- n2p2_utils_overview.png → manual/n2p2_utils

### 全量验证
- Sphinx 构建：成功，0 warnings，0 errors

### 全项目统计（Stage 0-6 累计）

- 新增示例脚本：5 个（docs/examples/）
- 新增文档页面：7 个（docs/source/manual/）
- 新增图片：5 张（docs/source/_static/images/generated/）
- 更新文件：index.rst, vasp.rst, generate_structure_images.py, docs.yml
- 覆盖脚本：62+ 个（vasp_universal 22 + vasp_workflow_bulk+neb 24 + slurm 8 + lmp 20 + n2p2 15）

### 部署
工作分支 docs/utils-website-upgrade 已完成所有修改。
main 分支保持不动。

push 到 origin 后，GitHub Actions 不会自动部署（CI 只在 push main 时触发）。
需要 merge 到 main 或修改 workflow 触发条件才能上线。

## 建议
1. 审阅 docs/utils-website-upgrade 分支的改动
2. 确认无误后 merge 到 main（或提 PR）
3. main 分支 push 后 GitHub Actions 自动部署
4. 部署后访问 https://shenglin1001.github.io/pjvasp_package/ 验收
