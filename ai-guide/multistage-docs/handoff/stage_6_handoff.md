# Stage 6 Handoff — 审查收尾 + 部署

## 完成时间
2026-08-02

## 完成内容

1. **首页 index.rst 更新**：新增 8 张 Stage 1-5 生成的教程图片画廊，
   每张图链接到对应 tutorial 页面。
2. **linkcheck 修复**：修复 findcubic.py 中 2 个 ASE 文档旧 URL
   （wiki.fysik.dtu.dk → docs.ase-lib.org）。
3. **全量验证**：
   - `compileall mymetal`：成功
   - `sphinx -E -b html -W --keep-going`：build succeeded，0 warning
   - `sphinx -b linkcheck`：原有 broken link 已修复

## 全项目统计（Stage 0-6 累计）

- 新增示例脚本：11 个（docs/examples/）
- 新增教程页面：11 个（docs/source/tutorials/）
- 新增/更新 API 页面：8 个（docs/source/api/）
- 新增图片：14 张（docs/source/_static/images/generated/）
- 修复源码 docstring RST 格式：5 个文件
  （extxyz.py, general.py, calhetero.py, plotting.py, post.py）
- 新增 Sphinx 扩展：viewcode, intersphinx, copybutton
- CSS 升级到 ASE 风格色板

## 部署

工作分支 `docs/mymetal-website-upgrade` 已完成所有修改。
main 分支保持不动。

push 到 origin 后，GitHub Actions 会自动构建并部署到 GitHub Pages。
但由于 CI workflow 的 `on.push.branches` 只包含 `[main]`，
push 到 `docs/mymetal-website-upgrade` 不会自动触发部署。
需要 merge 到 main 或修改 workflow 触发条件才能上线。

## 建议

1. 审阅 `docs/mymetal-website-upgrade` 分支的改动
2. 确认无误后 merge 到 main（或提 PR）
3. main 分支 push 后 GitHub Actions 自动部署
4. 部署后访问 https://shenglin1001.github.io/pjvasp_package/ 验收
