---
name: p-page-creator
description: 为 Python package 创建、维护和部署 Sphinx manual 风格的 GitHub Pages 官方文档站点。用于制作或更新 package 网站、重构为类似 LAMMPS Manual 的技术手册、发布 Pages、排查 Sphinx Actions 构建失败、验证远程文档页面，以及处理 companion package 的 autodoc 依赖。
---

# P Page Creator

## 概述

针对当前 Python package 仓库，创建或更新 Sphinx 官方文档站点，并通过 GitHub
Pages 部署为可远程访问的 manual 网站。目标不是简单 landing page，而是类似
LAMMPS Manual 的分章节技术手册，并完成本地构建、CI 部署和线上页面验收闭环。

## 一、适用场景
- 用户要求为当前 Python package 创建官方文档。
- 用户要求用 Sphinx 生成可浏览的 manual。
- 用户要求部署到 GitHub Pages 并给出远程访问 URL。
- 用户要求维护、补充内容或修复已经部署的 package 文档站点。
- 现有静态首页过于简单，需要重构成技术手册结构。
- 包含外部 companion package，autodoc 构建时必须正确处理其来源与安装方式。
- 本地 Sphinx 构建成功，但 GitHub Actions 或 Pages 发布失败。

## 二、工作原则
- 先检查当前仓库真实结构、已有 `docs/`、`setup.py` / `pyproject.toml`、模块目录和 GitHub remote。
- 只在用户指定分支上操作；若用户要求保持 `main`，不得切换到其他自动生成分支。
- 不要默认删除整个旧 `docs/`；优先在现有文档结构中增量重构。
- 不要把实际运行时依赖简单 mock 掉来掩盖问题。若依赖有明确来源，应在 CI 中安装该依赖。
- optional dependency 才适合放入 `autodoc_mock_imports`。
- GitHub Pages 是否可用、workflow 是否成功、远程 URL 是否真实返回目标页面，都必须实时验证。
- 若远端分支被其他工具修改，默认不要强推；只有用户明确授权覆盖时才使用 `--force-with-lease`，并先确认远端最新 hash。
- 若用户只要求更新现有站点，先阅读已有 `docs/` 与 workflow，只修改受影响页面或依赖，不重复搭建整套站点。
- 本 skill 的主流程对 Python package 通用；以下出现 `pjvasp_package`、`mymetal` 或 `myalloy_package` 的内容是已验证案例，应先核对目标仓库后再复用。

## 三、文档站点推荐结构
将 Sphinx 首页重构为 manual 目录入口：

`docs/source/index.rst`
- 标题：`<package> Manual`
- 简介：说明这是科学计算/技术软件手册
- 顶层分为三组 toctree：
  1. `User Guide`
  2. `Workflow Guide`
  3. `Reference`

推荐新增页面：

`docs/source/user_guide/`
- `overview.rst`：项目组成、核心模块、外部依赖
- `install.rst`：包安装、运行时依赖安装、文档构建安装
- `quickstart.rst`：最小 Python 示例
- `examples.rst`：仓库已有 example 的索引
- `troubleshooting.rst`：导入失败、HPC 命令缺失、Sphinx warning 处理

`docs/source/manual/`
- `workflows.rst`：整体工作流和目录约定
- `vasp.rst`：VASP 相关工具与提交/后处理入口
- `lammps.rst`：LAMMPS 模板与适用场景
- `slurm.rst`：调度器脚本检查项
- `n2p2.rst`：机器学习势数据流程

`docs/source/reference/`
- `scripts.rst`：脚本目录说明
- `dependencies.rst`：核心依赖、companion package、optional dependencies
- `development.rst`：本地文档构建和验证命令

## 四、Sphinx 配置要点
推荐配置：
- 使用 `sphinx_rtd_theme`，提供左侧导航形式的 manual 阅读体验。
- `html_theme_options` 推荐：
  ```python
  html_theme_options = {
      'navigation_depth': 4,
      'collapse_navigation': False,
      'sticky_navigation': True,
  }
  ```

- 开启常用 extension：

  ```python
  extensions = [
      'sphinx.ext.autodoc',
      'sphinx.ext.napoleon',
      'sphinx.ext.mathjax',
      'sphinx.ext.todo',
  ]
  ```

- 旧生成页或被替换页面可加入 exclude_patterns，避免重复导航和构建 warning。
- autodoc 页应区分：
    - 可直接导入的核心 API：使用 automodule。
    - 强依赖运行环境、在在线构建中不适宜 import 的模块：写明确的文本参考说明，不要为了展示完整 API 导致 Pages 构建失败。
    - 明确的外部运行时依赖：安装依赖，而不是伪装为本地源码或随意 mock。

## 五、外部依赖处理原则与案例

对任意 package，先从 import 失败日志和源码导入链确认 companion package 的
真实来源。若它是文档构建必须导入的运行时包，在 workflow 中安装它；若只需要
其 importable helper 而不需要完整计算栈，可审慎使用 `--no-deps`，并将文档
实际需要的直接依赖显式写入 `docs/requirements.txt`。

### `pjvasp_package` 已验证案例

`pjvasp_package` 中的 `mymetal` 部分模块依赖 `myalloy` 和 `myvasp`，它们不是
当前仓库内需要寻找的 Python 源包，而是来自 companion repository：

https://github.com/ShengLin1001/myalloy_package

该仓库的 master 分支同时提供 importable package：

- myalloy
- myvasp

该案例的文档 workflow 中应安装：

- name: Install myalloy runtime helpers
  run: python -m pip install --no-deps "git+https://github.com/ShengLin1001/myalloy_package.git@master"

采用 --no-deps 的原因：

- 文档构建只需要 importable helper 包本身。
- companion package 的完整依赖中可能包含较重或不适合 GitHub Pages runner 的运行时组件。
- 文档需要的直接科学计算依赖应明确写入 docs/requirements.txt。

该案例中 `mymetal` 的 autodoc 导入链还会加载 plotting helpers，因此
`docs/requirements.txt` 至少需要包含：

```text
brokenaxes>=0.6
adjustText>=1.3
```

`adjustText` 不能仅因本地环境已安装而遗漏。GitHub runner 在导入
`mymetal.ml.n2p2.calculate.sf` 时会加载
`mymetal.universal.plot.general` 的顶层 `adjustText` import；缺少该依赖时，
Sphinx 报出 `No module named 'adjustText' [autodoc.import_object]`，并在
`-W` 下使构建失败。

## 六、GitHub Pages Workflow 模板

以下模板提供 workflow 结构与已验证基线，不代表 action major version 永久
适用。创建或更新站点时：

1. 查看仓库现有 workflow 与最近 run 的 runner warnings。
2. 通过 GitHub 官方 action 仓库或 Pages 文档确认当前受支持 major version。
3. 如果日志提示 Node.js runtime 即将强制切换或 action 已弃用，先升级相关
   `actions/*` 引用并重新验证，再发布站点。

`pjvasp_package` 在 2026-05-27 成功部署时使用了下列 action 版本；同一次
runner 日志提示 `actions/checkout@v4`、`actions/setup-python@v5` 和
`actions/configure-pages@v5` 的 Node.js 20 执行环境将于 2026-06-02 起默认
切换到 Node.js 24，因此后续更新不得机械复制版本号。

```yaml
name: Publish Sphinx documentation

on:
  push:
    branches: [main]
    paths:
      - ".github/workflows/docs.yml"
      - "docs/**"
      - "<package>/**"
      - "setup.py"
  workflow_dispatch:

permissions:
  contents: read
  pages: write
  id-token: write

concurrency:
  group: pages
  cancel-in-progress: false

jobs:
  build:
    runs-on: ubuntu-latest
    steps:
      - name: Checkout
        uses: actions/checkout@v4

      - name: Set up Python
        uses: actions/setup-python@v5
        with:
          python-version: "3.10"
          cache: pip
          cache-dependency-path: docs/requirements.txt

      - name: Install documentation dependencies
        run: python -m pip install -r docs/requirements.txt

      - name: Install runtime helper package
        run: python -m pip install --no-deps "git+<companion-repo-url>@<branch>"

      - name: Configure Pages
        uses: actions/configure-pages@v5

      - name: Build Sphinx documentation
        run: python -m sphinx -b html -W --keep-going docs/source docs/_build/html

      - name: Disable Jekyll processing
        run: touch docs/_build/html/.nojekyll

      - name: Upload Pages artifact
        uses: actions/upload-pages-artifact@v3
        with:
          path: docs/_build/html

  deploy:
    environment:
      name: github-pages
      url: ${{ steps.deployment.outputs.page_url }}
    runs-on: ubuntu-latest
    needs: build
    steps:
      - name: Deploy Pages
        id: deployment
        uses: actions/deploy-pages@v4
```

## 七、本地验证门槛
提交前至少运行：

```bash
python -m sphinx -E -b html -W --keep-going docs/source /tmp/<project>-manual-build
python -m compileall -q <package>
git diff --check
```

验证渲染首页包含目标导航文本，例如：

- <package> Manual
- Manual contents
- User Guide
- Workflow Guide
- Reference

如果本地构建使用了已包含大量依赖的开发环境，而 CI 使用
`docs/requirements.txt` 创建干净环境，还必须将 autodoc 实际 import 到的
直接依赖与文档 requirements 对照检查。不要把“本地成功”当作 CI 依赖完备
的证据。

## 八、发布与部署验证门槛
本地成功不等于任务完成。必须继续验证：

1. 推送目标分支。
2. 检查最新 GitHub Actions Publish Sphinx documentation workflow。
3. 若失败，定位具体失败 step 和日志，继续修复。
4. workflow 成功后，实际访问 GitHub Pages URL。
5. 只有 URL 返回真实 manual 页面，且包含首页标题/导航文本，才可报告部署完成。

## 九、GitHub Actions 失败诊断闭环

当 Pages workflow 失败时，不要根据失败步骤名称猜测原因：

1. 确认失败 run、job 和具体 step，例如 `Build Sphinx documentation`。
2. 优先使用已认证的 GitHub connector 或 GitHub Actions 工具读取 job steps
   和完整日志；若环境中已有并认证 `gh`，可用 `gh run view <run-id> --log`。
3. 未认证 `curl` 请求 Actions job logs 可能返回 `403`。不能因此跳过日志
   读取，也不能仅根据本地构建成功宣布完成。
4. 如果可用 connector 只能按 PR 查 run 而部署由 `push` 触发，可从公开
   workflow 页面定位新 run ID，再交由已认证工具核验 jobs 与日志。
5. 从日志识别最小修复面：缺少直接依赖时补 `docs/requirements.txt`；只有
   明确 optional 且在线构建不应运行的依赖才加入 `autodoc_mock_imports`；
   RST warning 则修正文档本身。
6. 修复后重新运行第七节的严格验证，提交并推送目标分支，等待新 workflow
   的 `build` 和 `deploy` jobs 都成功。
7. 直接访问 Pages URL 并搜索 `<package> Manual`、`Manual contents` 等目标
   文本，记录部署验收结果。

## 十、提交和分支行为

- 提交前先检查实际 diff，并调用 `$p-git-commit` 生成符合仓库约定的中文
  Conventional Commit + gitmoji shortcode message；不要自行绕过该步骤。
- 根据变更意图，`$p-git-commit` 可生成类似以下 message：
  - 建站或内容重构：`docs: :memo: 重构 Sphinx manual 文档`
  - 修复文档 CI 依赖：`ci: :green_heart: 补充文档构建所需依赖`
- 提交前只暂存该次网站更新相关文件，保留无关工作区改动。
- 只在用户指定分支操作；用户限定 `main` 时，禁止切换到 Copilot 或其他临时分支。
- 对常规更新使用正常 push，不因旧失败 run 直接强推。

如果远端 `main` 被自动生成分支的 merge 覆盖，而用户明确授权覆盖，则：

1. 先 fetch 并确认远端 `main` 当前 hash。
2. 使用 `git push --force-with-lease=refs/heads/main:<confirmed-hash> origin main:main`。
3. 不使用无保护的 `--force`。

## 十一、完成报告

完成时报告：

- 实际修改的文档、Sphinx 配置或 CI 依赖文件。
- 本地验证命令及结果。
- 提交 hash 和推送分支（如果用户要求发布）。
- 成功的 Actions run URL 或 run ID，以及 `build` / `deploy` 结论。
- 已实际访问并验收的 GitHub Pages URL 和匹配到的首页文本。
