# 下一步任务：实施 pjvasp_package 文档网站第一轮升级

## 使用说明

在仓库根目录启动新的 Codex session，让新 session 完整读取并执行本文。

默认授权范围是：**修改、构建和本地验收，不 commit、不 push、不部署**。完成本地审阅后，如果确定要直接上线，再由用户明确追加本文末尾的“上线授权提示词”。

---

## 给下一 session 的任务提示词

你是一名熟悉科学计算软件、计算材料学工作流、Python API 文档和 Sphinx/GitHub Pages 的资深文档工程师。

当前仓库：

```text
F:\BaiduSyncdisk\version20240608\main_code_space\pjvasp_package
```

任务目标：

依据已经完成的审阅报告，实际实施 `pjvasp_package` 文档网站的第一轮可读性改版。这是“修改并验证”的任务，不要重新写一份审阅报告，也不要只给建议。

### 必须首先完整阅读

1. 仓库根目录适用的 `AGENTS.md` 指令；
2. `ai-guide/docs_readability_improvement_plan.md`；
3. `$p-page-creator` skill 的 `SKILL.md`；
4. 当前 `docs/source/`、`docs/requirements.txt`、`.github/workflows/docs.yml`、`README.md`；
5. 本轮涉及的真实 Python API 和 shell scripts。

开始前执行并记录：

```powershell
git status --short --untracked-files=all
git branch --show-current
git rev-parse HEAD
```

当前已知工作树状态：

- `setup.py` 存在用户原有修改，必须保留，不得修改、还原或暂存；
- `ai-guide/docs_readability_improvement_plan.md` 和本文是本次改版依据，必须保留；
- 只修改文档改版直接需要的文件；
- 禁止递归或批量删除；
- 本轮不要删除 `docs/build/`，只记录后续由维护者手工处理。

## 一、实施范围

以报告第 14 节“第一轮具体修改清单”为实施依据，按以下顺序完成。

### 1. 首页和信息架构

重写 `docs/source/index.rst`，使首页在首个视口内回答：

- `pjvasp_package` 是什么；
- `mymetal` 与 `pjvasp_package` 是什么关系；
- 适合哪些用户；
- 能完成哪些计算材料学任务；
- 新用户从哪里开始。

首页必须直接链接：

- Getting Started；
- Tutorials；
- Workflow Guides；
- API/Script Reference。

保留 `sphinx_rtd_theme`，不切换 Furo、PyData Theme，不引入复杂前端框架。

### 2. 完整 Getting Started

新增或调整：

```text
docs/source/getting_started/index.rst
docs/source/getting_started/au111_slab.rst
docs/examples/getting_started_au111.py
```

建立一个无需真实 VASP 即可完成的连续案例：

1. 使用真实入口 `mymetal.build.film.stretch.generate_film`；
2. 构建 12-layer FCC Au(111) slab；
3. 打印 formula、atom count、cell 和 PBC；
4. 使用 `mymetal.io.vasp.my_write_vasp` 写出 `POSCAR`；
5. 重新读取 `POSCAR` 并验证；
6. 使用 ASE/Matplotlib 生成结构图；
7. 提供确定性断言、预期终端输出和常见错误；
8. 页面使用 `literalinclude` 引用同一个脚本，禁止在页面和脚本中维护两份实现。

不得运行 VASP，不得要求 POTCAR。

### 3. 图片和轻量视觉改进

新增：

```text
docs/scripts/generate_structure_images.py
docs/source/_static/images/generated/
docs/source/_static/css/custom.css
```

要求：

- 从真实 Au(111) 示例生成非空结构图；
- 图片有 caption、alt text 和生成来源；
- 检查图片确实包含有效内容，不得提交纯白图；
- CSS 只处理首页内容宽度、任务入口、figure、宽表和移动端；
- 不做装饰性 banner；
- 不创建暂时不用的空目录或模板。

### 4. Installation、README 和 Examples

更新：

```text
README.md
docs/source/user_guide/install.rst
docs/source/user_guide/examples.rst
docs/source/reference/dependencies.rst
docs/source/reference/development.rst
```

要求：

- README 作为中文项目入口，不再重复完整 Tutorial 和 Workflow；
- Sphinx 是详细内容的 canonical source；
- 说明 Python package、optional dependency、external executable 和 shell script 的区别；
- `setup.py` 当前没有 `console_scripts`，不得声称 pip 会自动安装 `pei_*` 命令；
- 说明在 CentOS HPC 上缺少命令时先检查 `module avail`；
- Examples 页面必须提供可点击、经过验证的案例入口；
- 不得直接发布已知存在旧 import、保存错误或纯白图片的 notebook；
- 写明“中文叙事内容、英文 API/模块/命令名称”的语言策略。

### 5. 修正文档事实和核心 API

从当前源码验证并为以下入口提供可读 reference：

- `mymetal.build.film.stretch.generate_film`
- `mymetal.io.vasp.my_read_vasp`
- `mymetal.io.vasp.my_write_vasp`
- `mymetal.build.film.extrfilm.cal_area`
- `mymetal.calculate.calenergy.surfenergy.cal_surface_energy`
- `mymetal.post.newmain.PostTime`
- `mymetal.post.newmain.PostData`
- `mymetal.post.newmain.PostData2`
- `mymetal.slurm.submit.pei_slurm_univ_submit`

要求：

- 函数签名、参数、返回值、单位、异常和最小示例均以真实源码为准；
- 修复 `surfenergy.py` docstring 中面积和表面能单位乱码；
- 单位和 `factor` 含义必须根据代码与现有数据谨慎说明，不能虚构科研结论；
- 不把本轮实际失败的 `post_general()` 宣传为推荐入口；
- API 页面增加 Related tutorials；
- Tutorial 增加 Related API。

### 6. Workflow 和真实命令

重点更新：

```text
docs/source/manual/workflows.rst
docs/source/manual/vasp.rst
docs/source/manual/slurm.rst
docs/source/reference/scripts.rst
```

复用并核对：

```text
slurm_utils/README.md
slurm_utils/docs/submission-architecture-and-modes.md
vasp_utils/
slurm_utils/
lmp_utils/
```

要求：

- 所有命令必须对应当前真实存在的 `pei_*` 文件；
- 不继续把不存在的 `yin_*` 命令作为当前推荐命令；
- 说明 `parallel`、`each_subdir`、`single_alloc`；
- SLURM 示例第一遍必须是 dry-run；
- 不传 `-if_sbatch` 时不得提交作业；
- 不在本轮调用 `sbatch`；
- 不运行真实 VASP、LAMMPS 或 n2p2 training；
- EOS、NEB、GSFE、decohesion 暂时作为 Advanced Workflow 概览，不在缺少审核后 fixture 时伪造完整教程。

### 7. 两个 fixture-based Tutorial

在数据来源、许可和结果能够核对的前提下新增：

```text
docs/source/tutorials/surface_energy.rst
docs/source/tutorials/outcar_batch.rst
docs/examples/outcar_summary.py
```

要求：

- 优先复用仓库已有结构、OUTCAR 和结果文件；
- 不复制或发布 POTCAR；
- OUTCAR 只选择满足教学需要的最小 fixture；
- surface-energy 教程包含公式、输入、单位、结果表和手算/API 对照；
- OUTCAR 教程输出 convergence、energy、pressure、volume、force 等确定性表格；
- 每个教程必须包含前置条件、完整步骤、预期输出、验证方法、常见错误、下一步和相关 API；
- 如果数据授权或物理含义无法确认，停止发布该数据，记录具体 blocker，不得编造替代结果。

### 8. Sphinx 配置与 CI

按最小必要范围更新：

```text
docs/source/conf.py
docs/requirements.txt
.github/workflows/docs.yml
```

要求：

- 保留 `sphinx_rtd_theme`；
- 所有新增 extension 必须进入 `docs/requirements.txt`；
- 优先使用 Sphinx 内置能力；
- 可以增加 `intersphinx`、`viewcode`；只有确有收益时增加 `sphinx-copybutton`；
- 暂不引入 `sphinx-gallery`、`nbsphinx`、`myst-nb`、Mermaid 或新主题；
- workflow path filter 应覆盖：
  - `docs/**`
  - `mymetal/**`
  - `vasp_utils/**`
  - `slurm_utils/**`
  - `lmp_utils/**`
  - `myvasp/**`
  - `n2p2_utils/**`
- CI 增加无需 VASP 的 Getting Started smoke；
- 保持 GitHub Pages 自动部署；
- 不使用 mocked dependency 掩盖实际必须安装的 runtime dependency；
- 如需升级 GitHub Actions 版本，必须先查官方当前支持版本，不能猜测。

## 二、验证要求

修改过程中持续检查，不要等到最后一次性排错。

至少执行：

```powershell
& 'C:\Users\louis\mysoft\env\pyenv\codex\Scripts\python.exe' docs/examples/getting_started_au111.py --output docs/_build/example-au111

& 'C:\Users\louis\mysoft\env\pyenv\codex\Scripts\python.exe' -m sphinx -E -b html -W --keep-going docs/source docs/_build/html

& 'C:\Users\louis\mysoft\env\pyenv\codex\Scripts\python.exe' -m sphinx -b linkcheck -W --keep-going docs/source docs/_build/linkcheck

& 'C:\Users\louis\mysoft\env\pyenv\codex\Scripts\python.exe' -m compileall mymetal

git diff --check
```

还必须：

1. 实际打开构建后的首页、Getting Started 和两个 Tutorial；
2. 检查桌面导航、代码显示、交叉链接和图片；
3. 检查约 390 px 移动端没有页面级横向溢出；
4. 检查搜索 `surface energy`、`OUTCAR`、`SLURM`、`POSCAR` 能找到 canonical 页面；
5. 检查生成图片不是空白；
6. 检查没有发布 POTCAR、用户路径、账号或集群 secret；
7. 最后再次运行 `git status`，确认 `setup.py` 未被修改或暂存。

目标 CentOS HPC 若可用，还应执行：

```bash
module avail
bash -n <generated-slurm-script>
```

只有在维护者明确指定测试 partition 并授权后，才可执行真实 `sbatch`。

## 三、工作方式

- 先检查当前真实文件和 API，再修改；
- 使用 `apply_patch` 编辑文件；
- 不大规模重建整个 `docs/`；
- 不递归删除任何目录；
- 不覆盖用户已有改动；
- 不因为工作量大而只提交空框架；
- 页面必须有实际内容后才进入导航；
- 遇到错误应定位根因并继续修复；
- 对无法验证的科研结论明确标记，不得猜测；
- 保持在当前 `main` 分支，不创建临时发布分支；
- 默认不 commit、不 push、不部署。

## 四、完成时报告

请给出：

1. 实际修改和新增的文件；
2. 每个用户路径得到的改善；
3. 示例、Sphinx、linkcheck、compileall 的实际结果；
4. 未完成项及具体 blocker；
5. 当前 `git status`；
6. 是否满足 `ai-guide/docs_readability_improvement_plan.md` 第 15 节验收标准。

本轮默认先完成本地改版和验证，不要 commit、push 或部署，等待用户检查。

---

## 可选：上线授权提示词

只有在用户明确希望同一个 session 直接发布线上网站时，才把下面一段作为追加指令发送：

```text
现在授权继续上线。

本地全部验证通过后，使用 $p-git-commit skill 按功能边界提交文档改版。只暂存本轮相关文件、ai-guide/docs_readability_improvement_plan.md 和 ai-guide/next_task_docs_website_upgrade.md，不得暂存 setup.py。

正常 push 到 origin/main，不得 force push。推送前先 fetch 并确认本地 main 与 origin/main 没有未处理分歧。

随后检查最新 GitHub Pages workflow 的 build 和 deploy。实际访问在线文档并验收首页、Getting Started、Tutorial、搜索和约 390 px 移动端。只有 workflow 成功且线上页面确认更新后才算完成。

完成时报告 commit hash、push 分支、Actions run URL、build/deploy 结论和实际验收的 Pages URL。
```
