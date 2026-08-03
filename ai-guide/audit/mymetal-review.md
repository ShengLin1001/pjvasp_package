# mymetal 文档网站升级（Stage 0-6）独立审阅报告

**审阅对象**：`pjvasp_package` 仓库 `docs/mymetal-website-upgrade` 分支（已 push 到 origin）
**审阅日期**：2026-08-03
**审阅者**：独立审查 subagent（非原执行 agent）
**审阅方式**：真实 checkout 目标分支，用 `codex` venv python 实跑脚本/构建，读真实文件核实
**结论**：**产出质量高，无 blocker，可 merge**。发现 2 个 minor 问题 + 1 个统计口径备注。

---

## 审阅范围

在 `docs/mymetal-website-upgrade` 分支上 **git tracked** 的全部 Stage 0-6 产出（排除工作区 untracked 文件污染）：

| 类别 | 数量 | 路径 |
|------|------|------|
| 示例脚本 | 23 个 `.py` | `docs/examples/` |
| 教程页面 | 19 个 `.rst`（含 `index.rst`） | `docs/source/tutorials/` |
| API 参考 | 8 个 `.rst` | `docs/source/api/` |
| 生成图片 | 22 个 `.png` + 1 个 `dummy.sh`（脚本副产物） | `docs/source/_static/images/generated/` |
| 自定义 CSS | 1 个 | `docs/source/_static/css/custom.css`（300 行） |
| 图片生成脚本 | 1 个 | `docs/scripts/generate_structure_images.py`（214 行） |
| CI workflow | 1 个 | `.github/workflows/docs.yml`（127 行） |
| 计划/状态 | 7 个 | `ai-guide/multistage-docs/`（MASTER_PLAN / STAGE_STATUS / ASE_AESTHETIC_STUDY / SUBPACKAGE_COVERAGE_MAP / CRON_PROMPT + 7 个 handoff） |

**注**：审阅初期工作区混入了 `docs/utils-website-upgrade` 分支的 untracked 文件（`vasp_workflow_bulk_overview.py/.png`、`vasp_workflow_bulk.rst`、`neb_utils.rst`、`index.rst`/`vasp.rst` 的 modified）。这些**不属于 mymetal 分支产出**，已临时移走后重新构建核实，不计入本次审阅。

## 审查方法

1. **示例脚本可运行性**：用 `C:/Users/louis/mysoft/env/pyenv/codex/Scripts/python.exe` 实跑全部 23 个 tracked 示例脚本（`--output` 到临时目录），记录退出码、耗时、产出 PNG。
2. **图片非空白**：用项目 python + PIL + numpy 打开 22 个 committed PNG，计算 `mean(|rgb/255 - 1|)` 和非白像素占比，阈值 `size>5KB 且 meandev>0.002`。
3. **RST 审美**：抽查 5 篇 tutorial 详读（`cubic_cell_and_stretch` / `miller_index_and_density` / `slurm_script_generation` / `gsfe_models` / `io_extxyz_and_general`），并对全部 18 篇（除 index）做 10 项 ASE 规范要素矩阵检查。
4. **交叉链接**：用正则提取全部 74 个 rst 里的 `:doc:`（157 处）和 `:ref:`（0 处）引用，解析相对/绝对路径，核对目标存在性。
5. **安全约束**：扫全部示例脚本 + generate_structure_images.py + workflow + rst 里的 `sbatch`/`vasp_std`/`POTCAR`/`rm -rf`/`rmtree`/`setup.py` 模式，逐行核查上下文。
6. **Sphinx 构建**：在纯净 mymetal 分支工作区跑 `sphinx -E -b html -W --keep-going`，核实 0 warning。
7. **计划一致性**：读 MASTER_PLAN / STAGE_STATUS / ASE_AESTHETIC_STUDY / SUBPACKAGE_COVERAGE_MAP / stage_6_handoff，核对自报状态与实际产出。
8. **CI 注册完整性**：核对 docs.yml 的示例调用（23 处）和 generate_structure_images.py 的 render 注册（22 处）是否覆盖全部 tracked 产出。

## 发现的问题

### Blocker（0 个）

无。

### Major（0 个）

无。

### Minor（3 个）

#### M1. `outcar_summary.py` 与 `surface_energy.py` CLI 风格不符 docs/examples 约定

- **位置**：`docs/examples/outcar_summary.py`、`docs/examples/surface_energy.py`
- **现象**：这两个脚本用单横杠 `-output` / `-path_case` / `-cases`（`pei_*` CLI 风格），不接受标准 `--output`。其余 21 个 docs/examples 脚本都用标准 argparse `--output`（双横杠）。
- **影响**：违反 `mymetal-docs-examples` skill 明确的"docs/examples 用标准 `--output` 双横杠"约定。直接 `python outcar_summary.py --output X` 会报 `unrecognized arguments`（审阅实测 exit=2）。
- **缓解**：`.github/workflows/docs.yml` 已为这两个脚本写了正确的原生 CLI 调用（`outcar_summary.py mymetal/example/test-post/y_dir -cases 0.997 1.000 -output ...`、`surface_energy.py` 无参），所以 **CI 能跑通**，不影响部署。但用户照搬其他脚本的 `--output` 习惯会失败。
- **来源**：这两个脚本来自更早的 `f180f53` commit（"重构计算材料学文档站点"），早于 Stage 0-6，是遗留脚本被纳入 docs/examples 目录。Stage 0-6 执行 agent 未统一它们的 CLI。
- **建议**：低优先级。要么统一改成 `--output`，要么在对应 tutorial（`outcar_batch.rst` / `surface_energy.rst`）的"运行完整脚本"段明确标注这两个脚本的特殊参数。

#### M2. `surface_energy.rst` 和 `outcar_batch.rst` 无结果图、验证段简略

- **位置**：`docs/source/tutorials/surface_energy.rst`、`docs/source/tutorials/outcar_batch.rst`
- **现象**：这两篇 tutorial 没有 `.. figure::` 段（其余 16 篇都有）。`surface_energy.rst` 的"验证方法"段也较简略。
- **原因**：`surface_energy.py` 依赖外部 fixture（`test-surface-energy/`，含真实 VASP 结果），不生成图；`outcar_batch.rst` 对应的 `outcar_summary.py` 输出是 CSV 表格非图。两者性质上不适合产图。
- **影响**：ASE 审美规范第 5 条"复杂概念配图"在这两篇未落实，但属合理缺失（不是所有教程都需要图）。不影响正确性。
- **建议**：可选——为 `outcar_batch` 加一张 OUTCAR 解析流程示意图；为 `surface_energy` 加一张 slab/bulk 能量对比示意图（可用合成数据）。

#### M3. `periodic_table_and_arkel.rst` 交叉链接偏少

- **位置**：`docs/source/tutorials/periodic_table_and_arkel.rst`
- **现象**：全文 `:doc:` 引用较少（末尾"相关 API"段只有 `:func:` 无 `:doc:`），不像其他 tutorial 有"下一步"段密集交叉链接到相关教程。
- **影响**：ASE 审美规范第 7 条"交叉链接密集"在此篇略弱。不影响正确性。
- **建议**：可选——加"下一步"段链接到 `miller_index_and_density`（同属 universal 子包）或 `../manual/n2p2`（n2p2 也用元素属性）。

### 统计口径备注（非问题）

`stage_6_handoff.md` 自报"新增示例 11 个 / 教程 11 个 / 图片 14 张"，这是 **Stage 1-5 新增**的口径（Stage 0 之前已存在的 12 个示例/8 篇 tutorial/8 张图不计入）。mymetal 分支累计 tracked 总数是 23 示例 / 19 tutorial rst / 22 PNG。handoff 数字与实际无矛盾，只是口径需读者注意。

## 通过项

### ✅ 示例脚本可运行性（21/23 标准 + 2/23 原生 CLI = 23/23 可跑）

- 21 个 `--output` 风格脚本：全部 exit=0，耗时 2.0-7.7s，各产出 1 个 PNG（`periodic_table_and_arkel` 产出 2 个）。
- `outcar_summary.py` + `surface_energy.py`：用原生 CLI + 仓库内 fixture（`mymetal/example/test-post/y_dir/`、`mymetal/example/test-surface-energy/`，均 tracked）可跑通。
- 无 ImportError、无超时、无崩溃。

### ✅ 图片非空白（22/22 通过）

全部 22 个 committed PNG：
- 文件大小 68KB - 946KB，全部 > 5KB 阈值。
- `mean(|rgb/255-1|)` 全部 > 0.002（最低 `vasp_workflow_bulk_overview.png` 0.0446，最高 `periodic_table_heatmap` 0.2935）——均远超空白阈值，有实质内容。
- 分辨率合理（777×1989 到 3788×2677）。

### ✅ Sphinx 构建（0 warning）

纯净 mymetal 分支工作区跑 `sphinx -E -b html -W --keep-going`：
- `build succeeded`，退出码 0，0 warning。
- 50 个源文件全部构建，intersphinx（ase/python/numpy）inventory 正常加载。
- 与 stage_6_handoff 自报"build succeeded，0 warning"一致。

> 注：审阅初期因工作区混入 utils 分支的 `neb_utils.rst`（untracked），曾出现 1 个 `Inline literal start-string without end-string` warning。移走该外来文件后 warning 消失——**该 warning 不属于 mymetal 分支产出**。

### ✅ 交叉链接完整性（157/157 `:doc:` 有效）

- 74 个 rst 文件，157 处 `:doc:` 引用，0 处 broken。
- 22 个 `.. _label:` 定义备用，`:ref:` 未使用（0 处）。
- 相对路径和绝对路径（`/...`）解析全部命中真实文件。

### ✅ 安全约束全部遵守

| 约束 | 核实结果 |
|------|----------|
| 不调用真实 VASP/sbatch | ✅ `slurm_script_generation.py` 的 `CMD="vasp_std"` 是写进生成脚本的字符串常量，不执行；`if_sbatch=False` dry-run；`vasp_universal_overview.py`/`vasp_workflow_bulk_overview.py` 的 vasp_std 出现在图例/box 文字里 |
| 不发布 POTCAR | ✅ 所有 POTCAR 提及都是说明性文字（"需要 POTCAR"/"不分发 POTCAR"），无真实 POTCAR 文件或内容 |
| 不修改 setup.py | ✅ `git diff main docs/mymetal-website-upgrade -- setup.py` 为空；setup.py 完全未动 |
| 无递归删除命令 | ✅ 全仓库扫描 `rm -rf`/`rmdir /s`/`rd /s`/`del /s`/`Remove-Item -Recurse`/`shutil.rmtree` = 0 命中 |
| 不大规模重建 docs/ | ✅ Stage 0-6 增量添加，未覆盖用户已有改动 |
| 中文叙事英文 API | ✅ tutorial 中文叙事，`:func:`/`:class:`/`:mod:` 引用英文模块路径，参数名/代码英文 |

### ✅ ASE 审美规范落实（路径 A：rtd_theme + custom.css）

`custom.css`（300 行）完整落实 ASE_AESTHETIC_STUDY 路径 A：
- 白底 `#ffffff`、正文 `rgb(34,40,50)`、标题 `rgb(41,49,61)`、代码块 `#f3f4f5`+`1px solid #d1d5da`、橙色点缀 `rgb(232,146,23)`。
- 内容区 1120px、正文 16px/行高 1.65、代码 14px。
- admonition 左色条（note 蓝/warning 红/important 橙）、figure 居中带边框、表格浅灰边框无斑马纹。
- 移动端 600px 响应式（task-grid 单列、代码 13px、field-list 单列）。

### ✅ RST 内容结构（18/18 tutorial 达标）

全部 18 篇 tutorial（除 index）均含：
- 开头一句话说清功能 + `:Audience:/:Time:/:Requires:/:Runs VASP: No` 元数据头
- `.. code-block:: console` 运行命令 + `.. literalinclude::` 引用脚本
- "预期输出"段 + "验证方法"段 + "Related API"段（`:func:`/`:class:` 列表）
- 16/18 篇有 `.. figure::` 结果图（outcar_batch/surface_energy 因性质无图，见 M2）

抽查的 5 篇详读质量很高：背景解释清晰、mathjax 公式正确（如 `find_cubic` 面积守恒公式）、"常见错误"段实用、交叉链接到 manual 和相邻 tutorial。

### ✅ CI 注册完整性（100% 覆盖）

- `docs.yml` 注册 **23 个**示例调用，覆盖 mymetal 分支全部 23 个 tracked 示例（含 surface_energy/outcar_summary 的原生 CLI 调用）。
- `generate_structure_images.py` 注册 **22 个** PNG render，覆盖全部 22 个 tracked PNG。
- workflow 含 `compileall`（隐式，通过 pip install）、sphinx `-W` 构建、Pages 部署、linkcheck（定时任务）。

### ✅ 计划一致性

- STAGE_STATUS 自报 Stage 0-6 全部 completed（2026-08-02），与 git log 的 8 个 stage commit 一致。
- MASTER_PLAN 的 9 条安全约束全部遵守（见上）。
- SUBPACKAGE_COVERAGE_MAP 详细列出 6 个子包共 71 个模块的 LOC/函数/类/公开 API，与实际源码结构吻合。
- 每阶段 handoff 文档齐全（stage_0 到 stage_6 共 7 个）。

## 建议

1. **可 merge**：产出质量高，无 blocker/major，Sphinx 构建 0 warning，安全约束全守。建议 merge `docs/mymetal-website-upgrade` 到 main 触发 GitHub Pages 部署。
2. **M1 修不修两可**：`outcar_summary.py`/`surface_energy.py` 的 CLI 风格不影响 CI（workflow 已写对调用），但影响用户一致性体验。若要修，统一改成 `--output` 双横杠；若不修，在对应 tutorial 标注特殊参数。
3. **M2/M3 可选优化**：为 outcar_batch/surface_energy 加流程示意图、为 periodic_table_and_arkel 加"下一步"交叉链接——纯增强，非必需。
4. **部署注意**：如 stage_6_handoff 所述，CI workflow 的 `on.push.branches` 只含 `[main]`，push 到 mymetal 分支不会自动部署，需 merge 到 main。
5. **后续 stage**：mymetal 分支未覆盖的子包模块（如 `electronic_structure/plotter.py` 4586 LOC、`post/oldmain.py` 1285 LOC、`ml/n2p2/workflow.py` 1972 LOC 等大模块）可在未来 stage 补 API 文档。

---

*报告生成方式：独立 subagent 用真实文件和命令核实，未修改 mymetal 分支任何 tracked 文件。审阅过程中临时移走的 utils 分支 untracked 文件已恢复。*
