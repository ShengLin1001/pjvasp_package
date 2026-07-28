# pjvasp_package 文档网站可读性提升方案

> 审阅日期：2026-07-28
> 仓库：<https://github.com/ShengLin1001/pjvasp_package>
> 在线文档：<https://shenglin1001.github.io/pjvasp_package/>
> 本地基线：`main`，commit `267da73db8475ec2d2b64eeb8daf0f74dda3470a`
> 报告性质：审阅与实施规划；本轮未改写正式文档源码，也未运行真实 VASP、LAMMPS 或 SLURM 作业。

本文使用三种证据状态：

- **已确认**：已从当前文件、源码、构建输出、网页或实际运行结果核实。
- **合理推断**：由已确认事实推出，但尚未经过目标 HPC 或真实科研计算验证。
- **建议**：下一阶段应实施的设计或工程动作，不代表当前已有能力。

## 1. 执行摘要

当前文档工程本身是健康的，但内容架构还没有把仓库的实际能力转化为可学习、可验证的用户路径。

**已确认的正面基线**

1. 当前本地环境为 Python 3.10.7、Sphinx 8.1.3、`sphinx_rtd_theme` 3.1.0。严格 HTML 构建命令成功，21 个实际纳入构建的源页面没有 warning：

   ```text
   python -m sphinx -b html -W --keep-going docs/source docs/_build/html
   ```

2. `linkcheck` 严格构建也成功，当前已纳入文档的链接未发现失效项：

   ```text
   python -m sphinx -b linkcheck -W --keep-going docs/source docs/_build/linkcheck
   ```

3. 2026-07-28 的 GitHub Pages workflow 成功构建并部署当前 `main`：
   <https://github.com/ShengLin1001/pjvasp_package/actions/runs/30352865173>
4. 在线搜索可用；桌面侧栏、页面内标题锚点、自动 Previous/Next 均可用。390 px 宽移动端没有整页横向溢出，长代码块可以局部横向滚动。
5. 文档已经按 User Guide、Workflow Guide、API Reference 和 Reference 做了初步分区，方向正确，不需要推倒重建。

**已确认的核心缺口**

1. 首页标题是 `mymetal Manual`，项目仓库名是 `pjvasp_package`；首页没有清楚解释两者关系，也没有直接回答“适合谁、能完成什么任务、下一步点哪里”。证据见 `docs/source/index.rst` 和 `docs/source/conf.py:13-16`。
2. `Quick Start` 只有 `generate_film`、VASP I/O、surface energy、n2p2 四段互不相连的片段，没有完整输入、运行命令、预期输出、图、验证断言或下一步链接。证据见 `docs/source/user_guide/quickstart.rst:1-54`。
3. `Examples` 页面只列出五个目录文字，没有可点击案例、代码、结果或状态说明。仓库实际有 19 个 notebook，但部分 notebook 使用已失效的旧 import，至少一个保存了 `ImportError`，三个已导出的科学 PNG 均为 640×480 的纯白图。证据见 `docs/source/user_guide/examples.rst` 与 `mymetal/example/`。
4. 仓库存在真实的 FCC/HCP、film/stretch、surface energy、heterostructure、OUTCAR、EOS、NEB、GSFE、decohesion、SLURM、LAMMPS、n2p2 等能力，但现有 workflow 页面多为十几到几十行的概述；许多用户只能从源码、shell 脚本或各子目录 README 发现它们。
5. 导航中的 API 页只渲染了约 68 个 function、2 个 class 和 17 个 method，而静态源码盘点发现 80 个非 example 模块中有 426 个公开顶层 function/class。核心入口 `generate_film`、`create_fcc_111`、`create_hcp_basal`、`find_hetero` 和 SLURM 入口没有形成完整、可跳转的签名页。
6. 公开对象 docstring 覆盖率约为 83.1%（354/426），但格式不一致：约 248 个为 Google 风格、3 个为 NumPy 风格、103 个为非结构化文本；72 个公开对象没有 docstring，83 个有返回值的函数没有 `Returns` 段。`mymetal/calculate/calenergy/surfenergy.py:28-53` 还存在单位和中文注释乱码，属于科学含义风险，而不只是视觉问题。
7. 中文 README 已承担大量安装、架构和工作流说明，而正式文档主体为英文，当前没有语言策略和内容唯一来源，造成受众与维护责任不清。
8. `docs/build/` 中有 108 个已跟踪的旧 HTML/doctree 文件，而当前 CI 输出到 `docs/_build/html`；这些旧生成物不再参与部署，且最后更新时间早于当前源码，应该退出源码管理。
9. workflow 只监听 `docs/**`、`mymetal/**` 和 `setup.py`，没有监听 `vasp_utils/**`、`slurm_utils/**`、`lmp_utils/**`、`myvasp/**`、`n2p2_utils/**`。因此脚本接口改变时，文档不会自动重建。

**建议结论**

- **P0 不更换主题**。保留 `sphinx_rtd_theme`，先重写首页、建立一个真正完整的 Au(111) Getting Started、修正命令与单位事实、清理导航责任，并让示例进入 CI。
- **P1 建立 4–6 个完整科研任务教程**，优先使用无需 VASP 运行的结构构建、已有 OUTCAR 后处理、surface-energy 已有结果和 SLURM dry-run。
- **P2 再考虑 Gallery、notebook 集成、自动图片再生成和 API 覆盖率门禁**。当前 notebook 还不适合作为自动文档源，直接引入 `sphinx-gallery` 或 `myst-nb` 会先放大维护问题。

目标不是复制 ASE 的规模或样式，而是采用它最有效的原则：一个连续任务、可观察输出、解释性图片、可下载脚本、相关 API 和明确的下一步。

## 2. 审阅范围与方法

### 2.1 已审阅范围

按任务要求实际检查了：

- `README.md`
- `docs/`、`docs/source/`、`docs/requirements.txt`
- `.github/workflows/docs.yml`
- `mymetal/`、`mymetal/example/`
- `vasp_utils/`、`slurm_utils/`、`lmp_utils/`、`myvasp/`
- 与 n2p2 相关的 `mymetal/ml/n2p2/` 和 `n2p2_utils/`
- 当前在线站点的首页、Quick Start、Examples、API、搜索、桌面侧栏和 390 px 移动端
- ASE 当前官方文档首页、Tutorial index 和一个完整 introductory tutorial

### 2.2 审阅方法

| 方法 | 实际动作 | 用途 |
| --- | --- | --- |
| 文件与导航盘点 | 枚举 `docs/source/**/*.rst`，读取全部 21 个当前构建页；核对 45 个 RST 源文件与 `exclude_patterns` | 区分“存在于仓库”和“真实进入网站” |
| 严格构建 | 运行 HTML `-W --keep-going` | 验证 warning、导入和交叉引用 |
| 链接检查 | 运行 `linkcheck -W --keep-going` | 验证当前已纳入页面的链接 |
| 源码静态盘点 | AST 统计公开 function/class、docstring 与常见段落；定位教程候选真实入口 | 避免依据文件名虚构 API |
| 示例审阅 | 检查 19 个 notebook、输入/输出文件、PNG、OUTCAR、surface-energy 结果与 import | 判断可复现性和素材可用性 |
| 代表性运行 | 实际运行 `generate_film`、`my_read_vasp`/`my_write_vasp`、`cal_surface_energy`、`nnpdata` 和 OUTCAR parser | 验证 Quick Start 及案例入口，不只阅读代码 |
| 在线交互审阅 | 检查首页、搜索、代码块、跨页链接、桌面和移动端 | 区分源码问题与线上表现 |
| ASE 对照 | 使用当前 <https://docs.ase-lib.org/>，审阅首页、gallery index 和 `00-n2cu` 教程 | 提取信息架构与教学方法，不复制内容 |

### 2.3 已验证结果与限制

- `generate_film(symbols="Au", structure="fcc", num_layers=12, my_vacuum=20.0, slice_plane=(1, 1, 1))` 实际返回 `Au12` 的 `ase.Atoms`；这足以支持一个无需 VASP 的 Getting Started。
- 使用 `mymetal/example/test-post/1.000/CONTCAR` 实际完成 `my_read_vasp` 与 `my_write_vasp` 往返写入。
- 使用一个仓库内 OUTCAR 实际得到收敛状态、iteration、elapsed time、energy、stress、volume、pressure 和 maximum force。`PostTime`、`PostData`、`PostData2` 可作为教程入口。
- `post_general()` 在隔离目录的实际调用中因默认 `job_list=None` 路径触发 `TypeError`；因此在代码修复前，不应把它作为文档中的首选批处理 API。
- `find_hetero` 当前在本次隔离环境中无法导入，原因是缺少可选依赖 `hetbuilder`。异质结构案例只能列为条件性 P1/P2，不能声称已完成运行验证。
- 未运行真实 VASP、LAMMPS、n2p2 训练或 SLURM `sbatch`，也未对仓库中已有能量数值做独立物理复核。报告只把这些文件视为仓库现有演示数据。
- 当前平台实际为 Windows PowerShell；CentOS HPC 的 `module avail`、module 名称、调度器策略和文件系统行为尚未现场验证，必须由目标集群补验。
- 工作树开始时已有用户修改 `setup.py`；本轮保持不动。除本报告外没有修改正式源码或文档。

## 3. 当前文档工程概况

### 3.1 Sphinx、主题与扩展

| 项目 | 已确认事实 | 评价 |
| --- | --- | --- |
| 依赖范围 | `docs/requirements.txt:1-2` 为 `Sphinx>=7.2,<9`、`sphinx-rtd-theme>=2,<4` | 有兼容范围，但不是精确 lock；可重复性优于无约束，仍可能随 minor 版本变化 |
| 本机解析版本 | Sphinx 8.1.3、`sphinx_rtd_theme` 3.1.0、Python 3.10.7 | 当前构建成功 |
| 项目标识 | `docs/source/conf.py:13-16`：project=`mymetal`、release=`1.0.0` | 与仓库 `pjvasp_package` 的关系没有在首页解释 |
| 扩展 | `autodoc`、`napoleon`、`mathjax`、`todo` | 足够支撑当前基础站点；`todo` 暂未形成实际维护流程 |
| 主题 | `sphinx_rtd_theme`，navigation depth 4，sidebar 不折叠 | 桌面和移动端均可用；不是当前主要瓶颈 |
| 模板与静态资源 | 配置了不存在的 `_templates`；没有 `_static`、自定义 CSS 或正式图片目录 | 不会导致构建失败，但说明视觉和页面组件尚未建立 |
| API 导入 | `autodoc_mock_imports` mock 了 `hetbuilder`、`myvasp`、`ovito`、`pymatgen`、`torch` 等 | 保证构建，但掩盖“文档可构建”和“用户可运行”的差异；页面需标明可选依赖 |

### 3.2 文档入口与 toctree

当前入口 `docs/source/index.rst` 有三个 numbered toctree：

```text
mymetal Manual
├── User Guide
│   ├── Overview
│   ├── Install
│   ├── Quick Start
│   ├── Examples
│   └── Troubleshooting
├── Workflow Guide
│   ├── Workflows
│   ├── VASP
│   ├── LAMMPS
│   ├── SLURM
│   └── n2p2
└── Reference
    ├── API Reference
    │   ├── Build
    │   ├── Calculate
    │   ├── I/O
    │   ├── Post-processing
    │   ├── Universal
    │   └── Machine learning
    ├── Scripts
    ├── Dependencies
    └── Development
```

这一结构已经有“指南”和“参考”的意识，但 Tutorial、How-to 和 Workflow 的职责仍混在短页面中；首页的整个导航树承担了入口卡片的作用，用户要先理解仓库内部模块，才能找到任务。

### 3.3 构建与 GitHub Pages

`.github/workflows/docs.yml` 当前流程是：

```text
push main
→ checkout
→ Python 3.10
→ install docs/requirements.txt
→ install myalloy helper from Git
→ sphinx -b html -W --keep-going
→ upload Pages artifact
→ deploy-pages
```

结论：

- **已确认**：本地严格 HTML 构建、linkcheck 和最近一次线上 Pages 部署均成功。
- **已确认**：CI 使用 `docs/_build/html`，不是已跟踪的 `docs/build/`。
- **已确认**：workflow path filter 仅包含 `docs/**`、`mymetal/**`、`setup.py`（`.github/workflows/docs.yml:6-10`），脚本目录变化不会自动触发文档构建。
- **已确认**：最近一次日志提示多项 action 的 Node 20 runtime 即将废弃；涉及 `actions/checkout@v4`、`actions/setup-python@v5`、`actions/configure-pages@v5`、`actions/upload-pages-artifact@v3`、`actions/deploy-pages@v4`。这不是内容 P0，但应安排一次小型 CI 维护。
- **合理推断**：从 Git master 安装 `myalloy_package` 会使未来构建结果随上游变化；应固定已验证 commit 或 tag。

### 3.4 页面、孤立源与生成物

- `docs/source/` 有 45 个 RST 文件；当前 Sphinx 构建实际读取 21 个。
- `conf.py:29-40` 排除了旧 `installation.rst`、`quickstart.rst`、`workflows.rst` 和 `mymetal*.rst` sphinx-apidoc 页面，共 24 个不可见源文件。
- 严格构建没有 orphan warning，因为这些文件被明确 exclude；但它们仍形成双份内容和维护噪声。
- `docs/build/` 有 108 个 Git 已跟踪生成文件，最后内容早于当前文档源码；当前部署不使用它。
- 建议在迁移确认后由维护者**逐个明确移除或手工清理**旧生成物，并把 `docs/build/` 加入 ignore；不要在自动化任务中递归批量删除。

### 3.5 API 文档与 docstring

当前 API 是“部分 autodoc + 部分手工入口名”的混合模式：

- `docs/source/api/calculate.rst` 等页面对部分模块使用 `.. automodule:: ... :members:`。
- `docs/source/api/build.rst` 因 runtime `myvasp` 等依赖，手工写出 `create_fcc_111`、`create_hcp_basal`、`generate_film`、`find_hetero` 等名称，但没有自动签名、参数、单位和返回值。
- `mymetal.build.film.stretch` 与 `mymetal.build.film.findhetero` 的 automodule 没有 `:members:`，所以在线页不能展开核心函数。
- `mymetal.slurm.submit.pei_slurm_univ_submit` 没有进入当前 API 导航。

静态 docstring 审计结果：

| 指标 | 结果 | 说明 |
| --- | ---: | --- |
| 公开顶层 function/class | 426 | 80 个非 example Python 模块 |
| 有 docstring | 354（83.1%） | 覆盖尚可，但质量不等于覆盖 |
| 无 docstring | 72 | `mymetal/build/bulk/create.py` 的 11 个公开入口全部缺失 |
| Google 风格 | 约 248 | 当前主流，可作为统一方向 |
| NumPy 风格 | 3 | 与主流不一致 |
| 非结构化 | 103 | autodoc 可显示，但参数/返回值难稳定排版 |
| 有返回值但无 `Returns` | 83 | 科学结果类型和单位无法可靠理解 |

应优先修教程会直接链接到的 15–25 个入口，而不是立即把 426 个对象全部纳入 API。

### 3.6 README 与语言

`README.md` 为约 401 行中文说明，包含项目架构、安装、功能、workflow、示例和限制；正式 Sphinx 文档主体为英文。两者语义上重复安装、目录结构、Quick Start 和工作流说明，但没有声明哪一处是 canonical source。

建议采用当前维护成本最低的策略：

- README 保持中文项目入口，只保留定位、最短安装、一个最短示例和文档链接。
- Sphinx 的叙事型内容以中文为主，Python 名称、CLI、文件名和标准术语保留英文。
- docstring/API 参数继续使用技术英文，必要时在 Tutorial 中解释中文含义。
- 暂不维护中英文双份页面；只有当英文用户需求和维护者资源明确后，再引入语言切换。

## 4. 当前信息架构

### 4.1 逐页内容盘点

符号：`✓`=具备，`△`=部分具备，`—`=缺失或不适用。“完整代码”指页面自身可按顺序运行，而不是只出现一个函数调用。

| 当前页面 | 目标受众与任务 | 前置条件 | 完整代码 | 输入/目录 | 预期输出与结果含义 | 图/表/终端 | 页面内下一步 | 重复、过时或其他证据 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `index.rst` | 所有用户；网站入口 | — | — | — | — | — | 仅自动 toctree/Next | `mymetal Manual` 未解释与 `pjvasp_package` 的关系；线上首页是纯文字 |
| `user_guide/overview.rst` | 新用户；认识模块 | — | — | 目录名 △ | 模块用途摘要 △ | — | 仅自动导航 | 与 README 的项目架构重复；未转化为用户任务 |
| `user_guide/install.rst` | 新用户；安装 package | Python/pip △ | 命令 △ | 环境/路径 — | 安装验证 — | — | 仅自动导航 | 未解释 shell 脚本不由 `pip` console entry point 安装；HPC module 只有泛化说明 |
| `user_guide/quickstart.rst` | 新用户；尝试四项 API | 文件依赖未说明 | — | `CONTCAR`、`OUTCAR` 只写文件名 | 无预期输出、断言或物理解释 | 4 个 code block；无图/输出 | — | 四段孤立片段；`:9-54` |
| `user_guide/examples.rst` | 普通用户；寻找案例 | — | — | 只列 5 个目录 | — | — | — | 目录不可点击；实际 example tree 有 19 个 notebook |
| `user_guide/troubleshooting.rst` | 安装失败用户 | △ | 命令 △ | — | 只给一般处理 | — | — | 仅 myvasp/myalloy、`module avail` 和 docs warnings；没有 VASP/SLURM/OUTCAR 场景 |
| `manual/workflows.rst` | HPC 用户；理解总流程 | VASP/SLURM 隐含 | — | 简短目录树 ✓ | 流程产物 △；无验收 | 目录树 | — | 仅约 29 行；未涵盖中断恢复、归档和失败 |
| `manual/vasp.rst` | VASP 用户；发现工具 | VASP 隐含 | — | 目录/脚本 △ | — | — | — | 仅约 22 行；真实 EOS/NEB/GSFE/decohesion 等脚本未形成任务页 |
| `manual/lammps.rst` | LAMMPS 用户；发现模板 | LAMMPS 隐含 | — | 模板名 △ | — | — | — | 仅约 13 行；未解释 stretch、GSFE、Cij 的输入输出 |
| `manual/slurm.rst` | HPC 用户；提交任务 | SLURM 隐含 | — | — | — | — | — | 仅约 14 行；远少于 `slurm_utils/README.md` 和 architecture 文档中的真实能力 |
| `manual/n2p2.rst` | ML 用户；发现 n2p2 | n2p2 隐含 | — | 输入名 △ | — | — | — | 仅约 20 行；未连接 `nnpdata`、symmetry-function 工具与数据验证 |
| `api.rst` | 开发者；进入 API | Python 与可选依赖说明 △ | — | — | — | — | 到 6 个模块页 ✓ | 说明 mocked dependency，但未按对象标记 |
| `api/build.rst` | 结构建模用户/开发者 | ASE、myvasp/hetbuilder 部分隐含 | — | — | 函数用途短句 △ | — | — | 核心函数只列名称，无签名；在线约 16 个函数对象，未展开 `generate_film` 等 |
| `api/calculate.rst` | 计算属性用户 | △ | 最小示例 — | — | docstring 决定；单位存在乱码 | — | — | surface energy 的 `Å²`/`J/m²` 文本在源码中乱码 |
| `api/io.rst` | VASP I/O 用户 | VASP 文件知识隐含 | — | 文件名 △ | 返回 tuple/scale 解释不足 | — | — | 没有链接到可运行 I/O how-to |
| `api/post.rst` | 后处理用户 | OUTCAR/目录隐含 | — | △ | 部分 parser 字段可见；整体流程 — | — | — | `post_general` 实际默认路径调用失败；文档未区分稳定入口 |
| `api/universal.rst` | 开发者；通用工具 | 按模块隐含 | — | — | 由 docstring 决定 | — | — | 函数集合多，缺少任务入口和教程反向链接 |
| `api/ml.rst` | n2p2 用户/开发者 | n2p2 与数据格式隐含 | — | △ | 对象接口 △ | — | — | 没有从 OUTCAR 到 `input.data` 的完整例子 |
| `reference/scripts.rst` | 命令行/HPC 用户 | PATH、shell、集群环境未说明 | — | 只列目录 | — | — | — | `setup.py` 无 `console_scripts`；页面未解释实际调用方式和 `pei_*` 名称 |
| `reference/dependencies.rst` | 安装维护者 | — | — | — | 依赖分类 △ | — | — | mocked/build/runtime/optional 三种依赖边界仍不清 |
| `reference/development.rst` | 贡献者 | Git/Sphinx 隐含 | 命令 △ | — | 构建结果 △ | — | — | 缺少示例、图片、交叉引用和 docs CI 的贡献验收规范 |

### 4.2 示例与素材盘点

| 现有素材 | 已确认状态 | 可用于文档的条件 |
| --- | --- | --- |
| `test-generate-bulk` | 有 notebook 和 POSCAR | 提炼为独立脚本；明确 lattice 参数与结构断言 |
| `test-stretch` | 有 notebook 和一组 POSCAR | 去除研究过程噪声；生成应变前后图与目录树 |
| `test-surface-energy` | 有 FCC/HCP 结构、能量文本、notebook、EOS 输出 | 先确认数据来源与单位；禁止把结果包装成新的科研结论 |
| `test-post` | 有 7 组 OUTCAR/INCAR/KPOINTS/POSCAR/CONTCAR；parser 已实测 | 选 1–2 组小型、可授权 fixture，输出确定性表格 |
| `test-n2p2-sfparams` | notebook 执行慢并保存 NumPy deprecation warning | 只抽取小数据和快速路径；不直接接入 CI notebook 执行 |
| `test-hetbuilder-fixatom` | notebook 使用旧 import，保存过 `ImportError`；本机缺 `hetbuilder` | 先修 import 和依赖说明，再决定是否发布 |
| `test-stack` | 只有一个未执行 code cell | 不能作为教程证据 |
| 三个导出 PNG | 均为 640×480 纯白图 | 不应发布；从真实 notebook/数据重新生成 |
| notebook 内嵌 EOS 图 | 图像非空 | 可作为再生成参考；需由脚本生成并由领域维护者确认曲线含义 |
| `POTCAR` | example tree 中存在 `PAW_PBE Au 04Oct2007` 文件 | 不得复制到下载包或文档发布物；教程要求用户自行提供合法 POTCAR |

## 5. 主要问题及证据

| 问题 | 证据 | 用户影响 | 根本原因 | 建议 | 优先级 | 工作量 | 涉及文件 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 首页不能快速建立项目认知，且 `pjvasp_package`/`mymetal` 名称关系不明 | `docs/source/index.rst:1`；`conf.py:13-16`；在线首页 | 新用户无法在几十秒内判断软件用途和入口 | 产品定位与首页内容问题 | 首页明确“Python package + VASP/HPC workflow scripts”，按任务提供 4–6 个入口，并放一个结果预览 | P0 | 小 | `docs/source/index.rst`、`README.md` |
| Quick Start 不是连续、可验证的第一条成功路径 | `user_guide/quickstart.rst:1-54` 仅 4 段代码；无输出、图、链接 | 用户不知道运行目录、输入从哪里来、成功是什么样 | 教学内容问题 | 改为一个 Au(111) slab → 写 POSCAR → 渲染 → 断言的 10–15 分钟任务；其他片段移到 how-to/tutorial | P0 | 中 | `user_guide/quickstart.rst` 或新 `getting_started/`、`docs/examples/` |
| Examples 页只是目录索引，现有 notebook 也不能直接作为文档 | `user_guide/examples.rst`；19 个 notebook；旧 import、保存的 ImportError、纯白 PNG | 用户从页面无法读案例，也容易运行到过时代码 | 示例缺乏发布标准和 CI | 建立 curated case index；每个案例绑定脚本、输入、输出、图片和状态；研究 notebook 不直接发布 | P0 | 中 | `user_guide/examples.rst`、`mymetal/example/`、新 `docs/examples/` |
| 脚本/workflow 文档与真实接口脱节 | `vasp_utils/README.md:48-108` 使用不存在的 `yin_*`；仓库实际为 `pei_vasp_run_*`；manual 页仅 13–29 行 | 用户复制命令即失败，EOS/NEB/GSFE 等能力难发现 | 多处独立文档长期漂移 | 用源码脚本名生成/维护一张 canonical command table；Sphinx workflow 页为唯一正式说明，子目录 README 只链接过去 | P0 | 中 | `manual/*.rst`、`reference/scripts.rst`、`vasp_utils/README.md`、`slurm_utils/README.md` |
| 安装页没有解释脚本如何进入 PATH | `setup.py:30-31` 只有 packages/install requirements，无 `entry_points`/`console_scripts`；`install.rst` 未说明 | Python import 成功，但 `pei_*` 命令找不到 | package API 与 shell 工具分发方式未建模 | 明确三种调用：Python module、仓库内脚本显式路径、维护者配置 PATH；不要声称 pip 会安装命令 | P0 | 小 | `user_guide/install.rst`、`reference/scripts.rst` |
| 科学单位和 API 合同存在乱码/缺失 | `mymetal/calculate/calenergy/surfenergy.py:28-53`；72 个公开对象无 docstring；83 个 value-returning 函数无 `Returns` | 用户可能误读面积/表面能单位，无法确认返回类型 | 源码编码、docstring 规范和 API QA 问题 | 先修 surface-energy 单位并由领域维护者核对；统一 Google 风格；优先覆盖教程入口 | P0 | 小至中 | `surfenergy.py`、`build/bulk/create.py`、`post/newmain.py`、API RST |
| Tutorial、How-to、Workflow、API 没有形成角色分工和交叉链接 | 所有 21 个页面几乎没有显式 `:doc:`/`:func:` 链接；Quick Start 到 API 无入口 | 用户在“学任务”和“查参数”之间来回搜索 | 信息架构与交叉引用问题 | 建立四种页面模板；每个教程末尾链接稳定 API，每个 API 增加 Related tutorials | P0 | 中 | `index.rst`、`api/*.rst`、新 `tutorials/`、`howto/`、`workflows/` |
| README 与 Sphinx 内容重复，语言受众不明确 | 中文 `README.md` 约 401 行；Sphinx 主体英文；安装/架构/workflow 重复 | 内容容易不一致，用户不知道看哪一份 | 内容治理问题 | 采用“中文叙事文档 + 英文标识/API”；README 只做 gateway，Sphinx 为 canonical detail | P0 | 小 | `README.md`、`docs/source/index.rst`、`docs/source/reference/development.rst` |
| 文档构建未覆盖脚本接口变更，也没有示例 smoke test | `.github/workflows/docs.yml:6-10,46`；只 build HTML | 页面可以构建，但示例和真实命令可失效 | CI 只验证 Sphinx 语法 | 扩大 path filter；运行一个无 VASP 的 Getting Started smoke；独立 linkcheck；检查关键脚本名存在 | P0 | 中 | `.github/workflows/docs.yml`、`docs/examples/` |
| 旧生成 HTML 和排除源长期留在仓库 | 108 个 tracked `docs/build` 文件；45 RST 中仅 21 构建；`conf.py:29-40` | 搜索仓库时出现多份真相，维护者可能误改 | 历史构建方式残留 | 将 `docs/build/` 退出跟踪并 ignore；旧 RST 逐页迁移后再明确移除，不自动递归清理 | P1 | 小 | `.gitignore`、`docs/build/`、旧根级 RST |
| 核心 API 在 reference 中不可发现 | 源码 426 个公开顶层对象；在线 API 仅约 68 functions；build 页未展开核心函数；SLURM API 缺席 | 高级用户必须 grep 源码 | autodoc 配置与依赖耦合 | 建立 curated public API allowlist；对可选依赖模块使用手写 autosummary/stub 或安全导入，不全量暴露内部函数 | P1 | 中 | `api/*.rst`、`conf.py`、目标 docstring |
| 缺少有解释意义且可追踪的图片 | 21 个当前页面没有 `image`/`figure`；example 导出 PNG 是空白 | 结构、真空层、应变和 workflow 难以理解 | 无图片目录、生成脚本和审核规则 | 先生成 4–6 张确定性结构/流程图；每图有 alt、图注、数据来源和脚本 | P1 | 中 | 新 `_static/images/generated/`、`docs/scripts/` |
| 专有/大型科研输入的发布边界不清 | example 含 VASP `POTCAR`；`test-post` 有约 40 MB OUTCAR | 可能造成授权和仓库体积风险 | 示例数据治理缺失 | 下载包排除 POTCAR；只纳入最小、可授权 fixture；记录 provenance 和 checksum | P1 | 中 | `mymetal/example/`、新 `docs/examples/data/`、Development guide |

## 6. 与 ASE 文档的针对性比较

说明：任务给出的旧入口 <https://wiki.fysik.dtu.dk/ase/> 当前引导到 ASE 社区页；本次对照使用 2026-07-28 可访问的官方站点 <https://docs.ase-lib.org/>。重点参考其教学组织，不复制文字、素材或目录规模。

| 比较维度 | pjvasp_package 当前状态 | ASE 的有效做法 | 是否适合本项目 | 建议采用方式 |
| --- | --- | --- | --- | --- |
| 首页回答“是什么” | 一段包描述，标题只写 `mymetal Manual` | 首页直接解释 ASE 是什么，并用 Atoms/calculator/optimization 例子体现角色 | 直接适合 | 用一句定位解释 `pjvasp_package` 是仓库，`mymetal` 是 Python 核心；列出结构构建、VASP/HPC workflow、后处理、ML 数据四类任务 |
| 首页回答“适合谁” | 未显式区分用户 | 入口文字和导航分别服务新用户、calculator 用户、开发者 | 直接适合 | 首页列三类用户：结构建模初学者、VASP/HPC 用户、API/ML 高级用户 |
| 首页从哪里开始 | 依赖整棵 toctree | 首页直接展示最小完整代码和可观察结果 | 直接适合 | 首页保留 8–12 行 Au(111) 预览，并有明确 `Start here` 链接 |
| 首个有效结果 | Quick Start 四段孤立调用，无输出 | ASE 首页 H2 例子展示完整优化过程、终端输出和最终能量 | 直接适合 | Getting Started 必须生成结构、POSCAR、图片并做断言；不运行 VASP |
| Getting Started 连续性 | 页面内部没有故事 | introductory tutorial 从科研问题、建模、约束、计算到结果解释连续展开 | 直接适合 | 一个页面完成“构建 Au(111) slab → 检查真空/PBC → 写 POSCAR → 可视化” |
| Tutorials 的组织 | 无独立 Tutorials；Examples 只列目录 | 以任务命名，分类并提供缩略图 | 缩小后适合 | 首期只做 5–8 个精选案例；不复制 ASE 的 41 项 gallery 规模 |
| How-to 与 API 分离 | Quick Start、workflow、API 职责重叠 | Tutorial 教过程，Modules 查对象，CLI/GUI 单独查用法 | 直接适合 | Tutorial、How-to、Workflow、API 四类页面各有模板与相互链接 |
| 可观察输出 | 代码后没有文本/文件/图 | 教程展示 terminal output、能量值和结构图 | 直接适合 | 每教程至少有一种可验收输出：assert、目录树、表格、结构图或曲线 |
| 图片作用 | 当前导航页零图片 | 结构图直接解释吸附体系和几何关系 | 直接适合 | 只为结构、真空层、应变、流程和曲线配图；不放装饰性 banner |
| 下载与复现 | 只给 example 目录名 | 可下载 `.py`、`.ipynb` 和 zip | 缩小后适合 | 每个核心 tutorial 至少提供一个 `.py`；暂不提供全站 notebook zip |
| 页面内学习路径 | 主要依赖自动 Previous/Next | 教程中主动链接 Installation、Atoms、Constraints、Optimization 等相关概念/API | 直接适合 | 在步骤旁边链接相关 API，结尾固定 Next step/Related API |
| Gallery | 无 | Sphinx-Gallery 自动缩略图、下载和执行信息 | 条件性适合 | 等 3–5 个纯 Python 示例稳定并能在 CI 快速执行后再评估 `sphinx-gallery` |
| Tips/FAQ/Troubleshooting | 只有一页很短 Troubleshooting | 独立 Tips、FAQ，并从相关任务页可达 | 缩小后适合 | 先做 8–12 个真实错误条目，按安装、文件、HPC、后处理分类；无需复制完整体系 |
| Note/Warning/Tip | 当前页面没有 admonition | 教程用 note 提醒先安装、用说明块解释细节 | 直接适合 | 仅在 prerequisites、POTCAR 授权、单位、HPC 副作用和耗时处使用 |
| 搜索与导航 | 搜索可用，结果受内容稀薄影响 | 主题导航和语义标题让搜索命中任务页 | 直接适合 | 优先改标题和页面正文；不因搜索外观更换主题 |
| 视觉主题与生态 logo | RTD 主题稳定但首页留白较多 | ASE 当前站点更现代且展示庞大 calculator 生态 | 多数不适合 | 不复制主题和 logo wall；只用小型 CSS 改善首页宽度、任务入口和图注 |

### 6.1 可以直接借鉴的原则

1. 首页在一个视口内完成定位、用户、能力和下一步。
2. 第一个案例形成连续故事，并展示真实可观察输出。
3. Tutorial 教“为什么和怎样完成任务”，API 只负责“对象合同”。
4. 页面主动提供相关 API、下一步和可下载脚本。
5. 图片必须解释几何、过程或结果，而不是装饰。

### 6.2 需要缩小后采用

1. Gallery：先手写 5–8 张案例卡，案例稳定后再自动生成。
2. notebook：只为确有交互价值的高级教程提供，不把 19 个研究 notebook 全部接入。
3. Tips/FAQ：从真实 issue、运行错误和 HPC 经验积累 8–12 条，而不是先造空目录。
4. API：只公开稳定、用户会直接调用的入口，不追求一次展示 426 个对象。

### 6.3 当前不适合采用

1. 复制 ASE 的目录、文案、主题或大规模 calculator logo 墙。
2. 为“现代感”立即迁移到复杂前端或更换主题。
3. 直接执行需要 VASP license、集群、长时间 notebook 或私有数据的 gallery。
4. 同时维护中文和英文两套完整教程。

## 7. 目标用户与典型使用路径

### 7.1 目标用户

| 用户 | 核心问题 | 首要入口 | 成功标准 |
| --- | --- | --- | --- |
| 结构建模初学者 | “如何用 Python 构建并检查一个 VASP 结构？” | Getting Started | 在本机生成 Au(111) slab、POSCAR 和结构图，无需 VASP |
| 普通 VASP/HPC 用户 | “如何批量准备、提交、检查和汇总任务？” | Tutorials + Workflow Guides | 能理解目录合同，先 dry-run，再在集群提交并定位失败 |
| 后处理用户 | “已有 OUTCAR，怎样得到可核对的能量/压力/力表？” | OUTCAR tutorial + Post API | 从固定 fixture 生成确定性表格并知道字段单位 |
| 机器学习势用户 | “如何把 VASP 结果转换成 n2p2 数据？” | n2p2 tutorial/workflow | 生成 `input.data`，核对结构数、元素和必要字段 |
| 高级用户/开发者 | “函数签名、单位、异常和扩展点在哪里？” | API Reference + Development | 两次导航内到达稳定 API，并能从 API 返回相关 tutorial |

### 7.2 推荐新用户路径

```text
Home
→ Getting Started: Build and verify an Au(111) slab
→ Tutorials: strained films 或 surface energy
→ Related API
```

首页应直接链接完整 Getting Started，因此从首页到第一个完整示例只需一次主要导航操作，满足“最多两次”的验收要求。

### 7.3 推荐高级用户路径

```text
Home
→ Workflow Guides
→ VASP / SLURM / LAMMPS / n2p2

Home
→ API Reference
→ Build / I/O / Post / ML / SLURM
```

每个 workflow 页面应在关键步骤旁链接命令/API；每个 API 页的 Related tutorials 应把高级用户带回可运行场景。

### 7.4 内容所有权

- README：项目 gateway，不承担完整教程。
- Getting Started：唯一的新用户连续路径。
- Tutorial：端到端科研任务，含背景与结果解释。
- How-to：一个具体操作，不展开完整科研背景。
- Workflow：跨目录、外部工具和 HPC 生命周期。
- API：签名、类型、单位、异常、最小示例。
- `mymetal/example/`：原始研究/回归素材，不默认等同于发布教程。
- `docs/examples/`：经过整理、可下载、可在 CI 验证的文档示例。

## 8. 推荐的信息架构

### 8.1 目标导航树

```text
Home
├── Installation
├── Getting Started
│   └── Build, write and verify an Au(111) slab
├── Tutorials
│   ├── Compare FCC and HCP reference cells
│   ├── Generate a strained-film series
│   ├── Calculate surface energy from VASP results
│   ├── Summarize a batch of OUTCAR files
│   ├── Prepare and dry-run a SLURM batch
│   ├── Convert OUTCAR to an n2p2 dataset
│   └── Build a lattice-matched heterostructure [conditional]
├── Workflow Guides
│   ├── VASP directory lifecycle
│   ├── SLURM submission and recovery
│   ├── LAMMPS templates
│   ├── n2p2 data workflow
│   └── Advanced workflows: EOS, NEB, GSFE, decohesion and HOEC
├── How-to Guides
│   ├── Read and write POSCAR/CONTCAR
│   ├── Create a strain list
│   ├── Inspect convergence and forces
│   └── Run repository scripts from PATH or an explicit path
├── Examples and Gallery
├── Reference
│   ├── Python API
│   │   ├── Build
│   │   ├── Calculate
│   │   ├── I/O
│   │   ├── Post-processing
│   │   ├── Machine learning
│   │   ├── SLURM
│   │   └── Universal utilities
│   ├── Command-line and Script Reference
│   ├── File and Directory Contracts
│   └── Dependencies
├── Troubleshooting
└── Development
```

物理目录最多三层，侧栏展示最多两级。Advanced workflows 先做一张总览页；只有当某一流程有可发布 fixture 和维护者时，才拆成独立 Tutorial。

### 8.2 一级栏目解决的问题

| 一级栏目 | 用户问题 | 内容边界 |
| --- | --- | --- |
| Home | 这是什么、适合谁、能做什么、从哪开始？ | 定位、任务入口、一个结果预览；不重复完整安装 |
| Installation | 怎样在本地或 HPC 得到可验证安装？ | Python/package、可选依赖、脚本 PATH、`module avail`、验证命令 |
| Getting Started | 怎样在无需 VASP 的情况下得到第一个有效结果？ | 单一连续任务，不做功能大全 |
| Tutorials | 怎样完成一个真实科研任务？ | 背景、完整步骤、输出、意义、常见错误 |
| Workflow Guides | 多个脚本、目录、外部程序怎样协同？ | 生命周期、目录合同、集群、恢复、归档 |
| How-to Guides | 怎样解决一个具体操作问题？ | 短、目标单一、可直接复制；不讲长背景 |
| Examples and Gallery | 有哪些经过验证的结果可以浏览和下载？ | 精选案例卡、状态、预览和下载，不暴露原始 notebook 噪声 |
| Reference | 精确的接口和合同是什么？ | API、CLI、参数、单位、异常、文件格式 |
| Troubleshooting | 失败时如何定位？ | 真实错误信息、诊断命令、恢复路径 |
| Development | 如何贡献且不让文档漂移？ | 构建、示例、图片、docstring、CI 规则 |

### 8.3 现有页面处置

| 现有页面/目录 | 处置 | 目标位置与说明 |
| --- | --- | --- |
| `index.rst` | 重写 | 保留文件；改为任务型首页 |
| `user_guide/overview.rst` | 合并 | 精简内容进入 Home/About；模块清单放 Reference |
| `user_guide/install.rst` | 保留并重命名导航标题 | 作为 `Installation` canonical page |
| `user_guide/quickstart.rst` | 拆分/替换 | 用完整 Au(111) Getting Started 替代；I/O、surface energy、n2p2 片段迁移到对应页面 |
| `user_guide/examples.rst` | 重写 | 变为 curated Examples and Gallery |
| `user_guide/troubleshooting.rst` | 保留并扩展 | 按安装、文件、HPC、后处理分区 |
| `manual/workflows.rst` | 合并/重命名 | `workflows/index.rst`，只做整体地图 |
| `manual/vasp.rst` | 保留并扩展 | `workflows/vasp.rst` |
| `manual/lammps.rst` | 保留并扩展 | `workflows/lammps.rst` |
| `manual/slurm.rst` | 合并 | 复用 `slurm_utils/README.md` 与 architecture 文档的已验证内容，形成 `workflows/slurm.rst` |
| `manual/n2p2.rst` | 保留并扩展 | `workflows/n2p2.rst` |
| `api.rst`、`api/*.rst` | 保留，逐步校正 | 先加稳定入口、单位、相关教程，不立即全量重生成 |
| `reference/scripts.rst` | 重写 | 以真实 `pei_*` 脚本为 canonical 表格，注明调用方式 |
| `reference/dependencies.rst` | 保留并分类 | 区分 build、runtime、optional、external executables |
| `reference/development.rst` | 保留并扩展 | 增加文档贡献和 CI 验收 |
| 排除的 24 个 RST | 逐页处置 | 内容迁移确认后逐个明确移除；不要继续维护双份页面 |
| `docs/build/` | 退出源码管理 | 由维护者明确清理，部署只使用 `docs/_build/html` |

为降低第一轮改动风险，可以先只改变 toctree 和用户可见标题，保留大部分物理路径；确认链接稳定后再做目录重命名。

### 8.4 推荐 `index.rst` 导航草案

```rst
pjvasp_package documentation
=============================

``pjvasp_package`` combines the ``mymetal`` Python package with reusable
VASP, SLURM, LAMMPS and n2p2 workflow scripts for computational materials
science.

Start here
----------

* :doc:`Build and verify an Au(111) slab <getting_started/au111_slab>`
* :doc:`Choose a task-oriented tutorial <tutorials/index>`
* :doc:`Run a VASP/SLURM workflow <workflows/index>`
* :doc:`Look up a Python or script interface <reference/index>`

.. toctree::
   :maxdepth: 2
   :caption: Start

   installation
   getting_started/index

.. toctree::
   :maxdepth: 2
   :caption: Learn by task

   tutorials/index
   workflows/index
   howto/index
   gallery/index

.. toctree::
   :maxdepth: 2
   :caption: Reference and support

   reference/index
   troubleshooting/index
   development
```

首页正文还应有一张实际 Au(111) slab 结果图和四个任务入口；不要在首页直接展开整棵 API 树。

## 9. 高价值教程与案例规划

### 9.1 选择原则

优先选择同时满足以下条件的案例：

1. 对应明确的科研任务，而不是孤立函数；
2. 入口已从当前源码确认；
3. 能在不运行真实 VASP 的条件下完成全部或主要步骤；
4. 有可观察输出和可自动验证条件；
5. 能连接两个以上模块，帮助用户理解完整工作流；
6. 数据来源和发布许可可说明。

首期建议 7 个案例，其中前 6 个可以构成核心教程体系；异质结构案例需先解决可选依赖和旧 import。EOS、NEB、GSFE、decohesion 已确认存在真实脚本和 post-processing 代码，但目前缺少经过整理、可授权、可快速复现的教学 fixture，第一轮只进入 Advanced Workflow Guide，不承诺为完整 Tutorial。

### 9.2 案例 1：构建、写出并验证 FCC Au(111) slab

| 字段 | 内容 |
| --- | --- |
| 案例名称 | 构建、写出并验证 FCC Au(111) slab |
| 目标用户 | 初学者 |
| 学习目标 | 认识 `ase.Atoms`、slab、Miller plane、layer、vacuum 和 PBC；把结构写成 VASP POSCAR；用断言判断成功 |
| 前置条件 | Python 3.10、`pjvasp_package` editable install、ASE、Matplotlib；不需要 VASP |
| 实际入口 | `mymetal.build.film.stretch.generate_film`；`mymetal.io.vasp.my_write_vasp` |
| 输入 | 一个可下载脚本；参数 `symbols="Au"`、`structure="fcc"`、`num_layers=12`、`my_vacuum=20.0`、`slice_plane=(1, 1, 1)` |
| 执行步骤 | 创建输出目录 → 调用 `generate_film` → 打印 formula/cell/PBC → 检查 12 层结构 → 写 `POSCAR` → 用 ASE/Matplotlib 渲染 |
| 预期输出 | `Au12` 的 `ase.Atoms`、终端摘要、`POSCAR`、`au111_slab.png`；这些结果已用当前源码代表性运行确认 |
| 推荐图片 | slab 侧视图，标出原子层、真空区和 z 方向；另可放小型顶视图 |
| 验证方式 | `len(slab) == 12`；formula 为 `Au12`；输出文件存在且可由 `my_read_vasp` 重新读回；图中真空区可见 |
| 后续链接 | `generate_film` API、POSCAR I/O how-to、strained-film tutorial、ASE `Atoms` intersphinx |
| 实施难度 | 小 |

这是唯一的 P0 Getting Started。首页代码预览可以来自同一个脚本，避免两份实现。

### 9.3 案例 2：比较 FCC(111) 与 HCP basal 参考晶胞

| 字段 | 内容 |
| --- | --- |
| 案例名称 | 比较 FCC(111) 与 HCP basal 参考晶胞 |
| 目标用户 | 初学者、普通用户 |
| 学习目标 | 理解两个参考取向的晶胞矢量、原子数、层间堆垛和周期边界；学会比较而不混淆“结构类型”和“切面” |
| 前置条件 | Python、ASE；明确提供经过说明的 `a`，HCP 还需提供 `c` |
| 实际入口 | `mymetal.build.bulk.create.create_fcc_111`；`mymetal.build.bulk.create.create_hcp_basal` |
| 输入 | 一个脚本和由领域维护者确认的演示 lattice parameters；可参考 `mymetal/example/test-generate-bulk/` 的现有 POSCAR，但不直接假定其物理值为标准 |
| 执行步骤 | 显式设置 lattice parameters → 创建两种 `Atoms` → 打印 cell/volume/atom count → 写出 POSCAR → 并排渲染顶视图和侧视图 |
| 预期输出 | 两个 `Atoms`、两份 POSCAR、结构对照表、`fcc_hcp_cells.png` |
| 推荐图片 | 同一比例尺的 FCC/HCP 顶视图和侧视图，颜色一致，标注 a/c 与 stacking direction |
| 验证方式 | 检查 cell determinant 为正、atom count 与 size 一致、写入后可读回；教程明确提示当前函数默认 `a=None`/`c=None` 不能直接作为有效输入 |
| 后续链接 | bulk create API、film/slab tutorial、ASE cell/Atoms reference |
| 实施难度 | 小至中 |

### 9.4 案例 3：生成一系列 strained film 与批量目录

| 字段 | 内容 |
| --- | --- |
| 案例名称 | 为薄膜生成一系列应变构型和 VASP 目录 |
| 目标用户 | 普通 VASP 用户 |
| 学习目标 | 从一个基准结构生成明确的 strain-factor series；理解 cell deformation、是否缩放原子和目录命名的关系 |
| 前置条件 | 已完成 Getting Started；需要一个基准 `Atoms` 或 POSCAR；只生成输入时不需要 VASP |
| 实际入口 | `mymetal.build.film.stretch.stretch_list_along_direction_to_cell`；`mymetal.io.vasp.my_write_vasp`；素材 `mymetal/example/test-stretch/` |
| 输入 | 基准 POSCAR、显式 `stretch_factor_list`、方向列表、输出根目录；教程只选一个方向和 5–7 个 factor，避免参数爆炸 |
| 执行步骤 | 读基准结构 → 定义 factor → 生成 `Atoms` list → 按 factor 建目录并写 POSCAR → 输出 manifest CSV/文本 |
| 预期输出 | `strain_*/POSCAR` 目录、factor/cell/volume manifest、应变前后结构图 |
| 推荐图片 | baseline 与最大拉/压应变的 cell overlay；另用文本目录树显示批量输出 |
| 验证方式 | 输出结构数等于 factor 数；factor=1.0 的 cell 与基准一致；目标方向 cell length 单调变化；每个 POSCAR 可读回 |
| 后续链接 | stretch API、VASP directory workflow、SLURM dry-run tutorial、stretch post-processing |
| 实施难度 | 中 |

目录树应使用可复制的 text block，而不是截图；结构变化才使用图片。

### 9.5 案例 4：从已有 VASP 结果计算表面能

| 字段 | 内容 |
| --- | --- |
| 案例名称 | 从 bulk 与 relaxed slab 结果计算 surface energy |
| 目标用户 | 普通用户 |
| 学习目标 | 理解 bulk reference、surface atom count、surface area、双表面 factor 和单位转换；避免把函数调用当作物理模型本身 |
| 前置条件 | 已有 bulk/slab energy、结构和明确单位；教程使用仓库已有、经确认可发布的 fixture，不运行 VASP |
| 实际入口 | `mymetal.build.film.extrfilm.cal_area`；`mymetal.calculate.calenergy.surfenergy.cal_surface_energy`；素材 `mymetal/example/test-surface-energy/` |
| 输入 | bulk energy、bulk atom count、relaxed slab energy、slab atom count、area、`factor`、`energy_unit`；附输入文件树 |
| 执行步骤 | 读取结构 → 计算面积 → 从 fixture/小型 metadata 读取能量 → 显式写出公式 → 调用 API → 生成结果表 |
| 预期输出 | 一行输入/单位/结果表、可下载 CSV；若展示现有 FCC/HCP 曲线，只使用仓库真实数据并标注 provenance |
| 推荐图片 | slab 侧视图与面积矢量示意；可选表面能对结构或 strain 的曲线，前提是领域维护者确认数据含义 |
| 验证方式 | 手工公式与 API 结果在给定 tolerance 内一致；`eV/Å²` 与 `J/m²` 转换经过独立单元检查；输入面积大于 0 |
| 后续链接 | `cal_area`、`cal_surface_energy` API、surface-energy VASP workflow、单位说明 |
| 实施难度 | 中 |

发布前必须先修复 `surfenergy.py` docstring 中的乱码单位，并由计算材料学维护者确认公式、`factor` 语义和示例数据。不得从本报告推导或伪造新的能量值。

### 9.6 案例 5：批量读取 OUTCAR 并汇总收敛、能量、压力和力

| 字段 | 内容 |
| --- | --- |
| 案例名称 | 把一批 OUTCAR 汇总成可检查的结果表 |
| 目标用户 | 普通用户、高级用户 |
| 学习目标 | 理解 OUTCAR parser 的输入目录、各字段来源、收敛判定和缺失文件处理；得到可进一步绘图的表格 |
| 前置条件 | 1–2 个可发布 OUTCAR fixture；无需 VASP executable |
| 实际入口 | `mymetal.post.newmain.PostTime`、`PostData`、`PostData2`；shell 对照入口 `vasp_utils/vasp_universal/pei_vasp_univ_post`；素材 `mymetal/example/test-post/` |
| 输入 | `cases/<factor>/OUTCAR` 目录；固定 manifest；不要直接把约 40 MB 全量 fixture 复制进文档发布物 |
| 执行步骤 | 扫描明确的 case list → 对每个目录读取 time/data/data2 → 规范化列名和单位 → 写 CSV/RST table → 演示一个缺失 OUTCAR 的可解释错误 |
| 预期输出 | case、converged、iteration、elapsed time、energy、stress、volume、pressure、maximum force 表；本轮已从单个 fixture 实际解析出这些字段 |
| 推荐图片 | 不优先用截图；正文使用可搜索的表格和 terminal code block。只有后续做 convergence curve 时再生成 SVG |
| 验证方式 | 对一个固定 OUTCAR 把关键字段与原文行交叉核对；行数等于 manifest case 数；缺失文件必须被标记而不是静默跳过 |
| 后续链接 | Post API、VASP post-processing workflow、plot workflow API、Troubleshooting |
| 实施难度 | 中 |

当前不应把 `post_general()` 作为教程入口：其默认调用路径在本轮实测中触发 `TypeError`。教程先使用已验证的 parser class；代码修复后再评估统一入口。

### 9.7 案例 6：组织并 dry-run 一个 SLURM 批量任务

| 字段 | 内容 |
| --- | --- |
| 案例名称 | 为多个计算目录生成 SLURM 脚本，并在提交前检查拓扑 |
| 目标用户 | 普通 HPC 用户、高级用户 |
| 学习目标 | 理解 `parallel`、`each_subdir`、`single_alloc` 三种模式、chunk、parent/worker、wall time 与 dry-run |
| 前置条件 | Python；一个包含 3 个 toy case 的目录；生成脚本不需要 SLURM，真正提交需要目标集群和 `sbatch` |
| 实际入口 | `slurm_utils/slurm_universal/pei_slurm_univ_submit.py` → `mymetal.slurm.submit.pei_slurm_univ_submit` |
| 输入 | toy `dir/case001..003`、显式 mode/chunks/command/module profile；第一遍不传 `-if_sbatch` |
| 执行步骤 | 显示 preset → dry-run 生成脚本 → 检查每目录 `sub_slurm_univ.sh` 与 `slurm/` parent/worker → 解释资源 → 只在 HPC 小节演示如何显式提交 |
| 预期输出 | 生成的 shell scripts 和目录树；dry-run 不调用 runner/launcher、不产生 job ID。现有 `slurm_utils/README.md:39,68-78,97-110` 与 architecture 文档确认该行为 |
| 推荐图片 | 三种 mode 的简洁 SVG 拓扑图；目录树用 text；SBATCH header 用 code block |
| 验证方式 | 不传 `-if_sbatch` 时没有提交副作用；脚本数、case list 和 chunk 分配可核对；shell 可做 `bash -n`（目标 CentOS 再验证） |
| 后续链接 | SLURM workflow、script reference、recovery/troubleshooting、VASP runner |
| 实施难度 | 中 |

教程不得在 CI 中调用 `sbatch`。真实提交步骤应放在 clearly marked HPC-only section，并以 `.. warning::` 说明副作用。

### 9.8 案例 7：把 VASP OUTCAR 转换为 n2p2 dataset

| 字段 | 内容 |
| --- | --- |
| 案例名称 | 从一个 VASP OUTCAR 生成并检查 n2p2 `input.data` |
| 目标用户 | 普通用户、高级 ML 用户 |
| 学习目标 | 理解 structure count、element、energy/force/lattice 字段和 tag；学会检查转换结果，不涉及训练 |
| 前置条件 | 一个可发布 OUTCAR fixture；Python/ASE；生成数据不需要 n2p2 executable |
| 实际入口 | `mymetal.ml.n2p2.dataset.nnpdata.load_from_outcar`、`nnpdata.write`；对称函数入口 `mymetal.ml.n2p2.calculate.sf.mysfparams` 只作为后续链接 |
| 输入 | `OUTCAR`、`index`、`tag`、输出路径；素材可从 `test-post` 选择一个最小 fixture |
| 执行步骤 | 创建 `nnpdata` → 加载 OUTCAR → 打印结构数/元素/tag → 写 `input.data` → 读取头尾记录并检查 block |
| 预期输出 | `input.data` 和 terminal summary；当前 Quick Start 调用已经在单个 fixture 上实际写出文件 |
| 推荐图片 | VASP outputs → parser → normalized structures → n2p2 dataset 的 SVG；正文展示一段短格式示例 |
| 验证方式 | 结构数与输入 fixture 一致；每个 block 有 begin/end；元素、atom count、energy/force 字段存在；重复运行结果确定 |
| 后续链接 | `nnpdata` API、n2p2 workflow、symmetry-function guide、dataset troubleshooting |
| 实施难度 | 中 |

### 9.9 条件性案例：构建 lattice-matched heterostructure

| 字段 | 内容 |
| --- | --- |
| 案例名称 | 搜索匹配超胞并构建 heterostructure |
| 目标用户 | 高级用户 |
| 学习目标 | 从 bottom/top primitive cells 搜索匹配、构建 supercells、stack 并解释 mismatch |
| 前置条件 | `hetbuilder`、ASE、经过确认的两份结构；需先修现有 notebook 的旧 import |
| 实际入口 | `mymetal.build.film.findhetero.find_hetero`、`build_supercells`、`stack_atoms` |
| 输入 | 两个 `ase.Atoms`/POSCAR、面积/角度/mismatch 搜索参数、stack distance/vacuum |
| 执行步骤 | 读结构 → 搜索候选 → 选择一项 → 构建 supercell → stack → 写结构 → 可视化 |
| 预期输出 | candidate table、bottom/top supercells、stacked POSCAR、结构图 |
| 推荐图片 | bottom/top primitive、匹配 supercell 和 stacked side view 三联图 |
| 验证方式 | 本地依赖可导入；candidate mismatch 与输出 cell 可独立计算；无重叠原子；POSCAR 可读回 |
| 后续链接 | heterostructure API、dependency guide、fix-atom how-to |
| 实施难度 | 大 |

本轮本机因缺少 `hetbuilder` 未能运行该入口，因此不能进入第一轮验收。

## 10. 配图与可视化规划

### 10.1 配图清单

| 对应页面 | 图片表达的信息 | 真实来源/生成入口 | 格式与建议尺寸 | 图注与 alt | 生成脚本 | CI 可重现 |
| --- | --- | --- | --- | --- | --- | --- |
| Home / Getting Started | Au(111) slab、12 层、z 向真空 | `generate_film` 的同一 tutorial script | PNG 或 WebP，1600×900；侧视图主体占 70–80% 宽度，保留真空 | 必须；alt 如“Au(111) 12-layer slab with vacuum along z” | `docs/scripts/generate_structure_images.py` | 是，无 VASP |
| FCC/HCP tutorial | 两种参考 cell 的堆垛与 cell vectors | `create_fcc_111`、`create_hcp_basal` | PNG，1800×900 两列；相同缩放和配色 | 必须，图注解释比较维度 | 同上 | 是，需固定 lattice input |
| Strained-film tutorial | baseline、压缩、拉伸三个 cell 的变化 | `stretch_list_along_direction_to_cell` | PNG，1800×700；三联图；避免过大空白 | 必须；alt 指明 strain factors 与方向 | 同上 | 是 |
| Surface-energy tutorial | slab 面积、双表面 factor 与真空 | `cal_area` + fixture structure | SVG 示意图 1200×700；结构可另用 PNG | 必须；不能用图形替代公式和单位文本 | `docs/scripts/generate_workflow_diagrams.py` 或受版本控制 SVG 源 | 是 |
| Surface-energy/EOS 结果 | 真实数据点、fit 或 surface-energy comparison | `test-surface-energy` 的已审核数据；不能使用当前纯白 PNG | SVG 优先、另导出 1600×1000 PNG；坐标轴和单位完整 | 必须，注明数据来源、仅作演示、fit 方法 | 独立 `docs/scripts/plot_surface_energy_example.py` | 条件性；数据需先审核 |
| OUTCAR tutorial | 多 case 的结果和异常状态 | `PostTime`/`PostData`/`PostData2` 输出 | 优先原生 RST/HTML table，不制作截图 | 表格 caption；列头含单位，状态不可只靠颜色 | `docs/examples/outcar_summary.py` 生成 CSV；Sphinx include | 是，需小 fixture |
| SLURM workflow | 三种 mode 的 parent/worker/job topology | `slurm_utils/docs/submission-architecture-and-modes.md` | SVG，1600×900；每种 mode 独立 panel | 必须；alt 描述作业关系 | 简单脚本或可审阅的 SVG 源；不依赖浏览器截图 | 是 |
| n2p2 tutorial | OUTCAR 到 `input.data` 的数据流和检查点 | `nnpdata.load_from_outcar`/`write` | SVG，1400×650 | 必须；alt 描述转换阶段 | `docs/scripts/generate_workflow_diagrams.py` | 是 |
| Heterostructure tutorial | primitive → matched supercells → stacked structure | `find_hetero`、`build_supercells`、`stack_atoms` | PNG，1800×900 三联图 | 必须，标注 bottom/top、mismatch、vacuum | `docs/scripts/generate_heterostructure_images.py` | 条件性；需 `hetbuilder` |
| Advanced workflow guide | EOS、NEB、GSFE、decohesion 的输入→计算→post map | 对应 `pei_vasp_run_*`、`mymetal.post.*` 与真实目录约定 | SVG，1600×1000；只画流程，不造结果曲线 | 必须；alt 列出阶段 | 手工可审阅 SVG 或一个小型生成脚本 | 是 |

### 10.2 图片规则

1. 所有结构图使用统一元素颜色、相机方向、背景和比例尺；同页比较图必须保持一致。
2. 所有曲线必须显示 quantity、unit、数据点与拟合/插值的区别；caption 说明数据来源和是否为仓库 fixture。
3. 每张图片同时保存生成脚本、输入文件清单和输出路径。只提交最终图片与小型输入，不提交临时帧。
4. 纯 Python/ASE/Matplotlib 图片应可在 CI 重现；涉及真实 VASP 的图片使用经过审核的固定结果，不在 CI 重算。
5. 终端输出优先用 `.. code-block:: text`，结果表优先用 RST/CSV table。截图只在解释图形界面或浏览器操作时使用。
6. PNG/WebP 用于原子结构等复杂 raster；SVG 用于流程图、示意图和曲线。建议网页显示宽度 700–1000 CSS px，并保留 2× 像素密度。
7. `alt` 不能写“figure”或重复 caption；应描述用户从图中需要获得的信息。
8. 不发布 POTCAR、许可不明确的外部数据、用户路径、账号、集群信息或未经确认的科研结果。

建议目录：

```text
docs/source/_static/images/generated/   # 进入网站的审核后图片
docs/examples/                          # 可下载、可运行示例
docs/examples/data/                     # 最小且可授权的固定输入
docs/scripts/                           # 图片与结果表生成脚本
```

## 11. 页面模板

### 11.1 Tutorial 页面模板

固定顺序：

1. 目标与科研背景；
2. 预计时间；
3. 前置条件和外部依赖；
4. 最终结果预览；
5. 输入文件和起始目录树；
6. 完整代码/命令；
7. 分步骤解释关键参数；
8. 预期 terminal、文件、表格或图；
9. 结果的物理/计算含义；
10. 验证断言；
11. 常见错误；
12. 下一步；
13. 相关 API；
14. 下载脚本与数据来源。

推荐 reStructuredText 骨架：

```rst
.. _tutorial-au111-slab:

Build and verify an Au(111) slab
================================

:Audience: Beginners
:Time: 10--15 minutes
:Requires: Python 3.10, pjvasp_package, ASE, Matplotlib
:Runs VASP: No

Goal
----

Build a 12-layer Au(111) slab, add vacuum, write ``POSCAR``, render the
structure, and verify the result.

Final result
------------

.. figure:: /_static/images/generated/au111_slab.png
   :width: 800px
   :alt: Au(111) 12-layer slab with vacuum along the z direction

   The generated slab. The empty region along z is the vacuum layer.

Prerequisites
-------------

Complete :doc:`../installation` and verify that ``import mymetal`` succeeds.

.. note::

   This tutorial does not run VASP and does not require a POTCAR.

Input and starting directory
----------------------------

.. code-block:: text

   au111-tutorial/
   └── build_au111.py

Run the example
---------------

Download :download:`build_au111.py <../../examples/getting_started_au111.py>`,
then run:

.. code-block:: console

   $ python build_au111.py --output results

The complete source is:

.. literalinclude:: ../../examples/getting_started_au111.py
   :language: python
   :linenos:

Expected output
---------------

.. code-block:: text

   formula: Au12
   atoms: 12
   wrote: results/POSCAR
   wrote: results/au111_slab.png

Verify the result
-----------------

Explain the atom-count, cell, PBC, vacuum and POSCAR round-trip checks here.

Common problems
---------------

``ModuleNotFoundError``
   Link to the installation diagnosis.

The structure has no visible vacuum
   Explain the viewing direction and cell extent.

Next steps
----------

Continue with :doc:`../tutorials/strained_films`.

Related API
-----------

* :func:`mymetal.build.film.stretch.generate_film`
* :func:`mymetal.io.vasp.my_write_vasp`
```

### 11.2 How-to 页面模板

How-to 只回答一个具体问题，建议控制为：

```text
标题：How to read a CONTCAR and preserve its scale
问题
最短可用命令/代码
输入假设
参数选择
预期输出
一个验证检查
两个最常见错误
相关 API / 完整 Tutorial
```

不重复解释完整 surface-energy 理论或 SLURM topology；这些内容链接到 Tutorial/Workflow。

### 11.3 API Reference 页面模板

每个稳定公开入口至少包含：

- 一句话用途与适用边界；
- 完整签名；
- 参数名、类型、允许值、默认值；
- 每个物理量的单位；
- 返回类型、shape/tuple 结构和单位；
- 可预期异常与错误条件；
- 10–20 行的最小示例及其简短输出；
- `Related tutorials`；
- 可选依赖标记；
- 对产生文件或副作用的函数明确列出 side effects。

示例骨架：

```rst
``cal_surface_energy``
----------------------

.. autofunction:: mymetal.calculate.calenergy.surfenergy.cal_surface_energy

Units
~~~~~

``bulk_energy`` and ``relaxed_surface_energy`` use ... .
``area`` uses ... . The return value is ... .

Minimal example
~~~~~~~~~~~~~~~

.. code-block:: python

   # Use reviewed, documented values here.

Related tutorials
~~~~~~~~~~~~~~~~~

See :doc:`../tutorials/surface_energy`.
```

单位文本必须由领域维护者核对，不能只依赖 type annotation。`napoleon` 已支持当前主流 Google-style docstring，因此没有理由同时引入第二种 docstring 风格。

### 11.4 Workflow 页面模板

Workflow 页固定回答：

1. 整体目的和不适用场景；
2. 组件关系：Python API、shell script、VASP/LAMMPS/n2p2、SLURM；
3. 起始输入；
4. 目录和文件合同；
5. dry-run 或 preview；
6. 真正提交/运行步骤；
7. 运行中状态和日志；
8. 中断、恢复、重提和幂等性；
9. 常见失败、诊断命令和退出码；
10. 输出、归档和清理边界；
11. 相关 Tutorial、How-to、CLI/API。

涉及集群的页面还必须列出：

- 目标 scheduler 和已验证集群；
- 所需 `module avail`/`module load`（按集群给实例，不写成通用真理）；
- launcher、partition、node/core、wall time；
- 是否会调用 `sbatch`、清理或覆盖文件；
- dry-run 默认值；
- 如何判断“提交失败”“作业失败”“计算未收敛”三种不同状态。

### 11.5 是否迁移到 MyST Markdown

**当前不建议迁移。**

收益：

- Markdown 对部分贡献者更熟悉；
- fenced directive、表格和 notebook 生态更顺手；
- 与未来 `myst-nb` 有自然衔接。

成本与风险：

- 当前 21 页 RST 已无 warning 构建；
- 迁移会产生大量无内容收益的 diff，并增加交叉引用/指令转换验证；
- 现有 notebook 本身尚未清理，先迁移格式不会提高可复现性；
- 混用 RST/MyST 会增加维护规则。

只有当维护者明确决定新内容主要由 Markdown 作者贡献，且至少 5 个教程已经稳定后，再做一个页面的 MyST prototype。当前保持 RST 是最小、可维护方案。

## 12. Sphinx 技术改进建议

| 工具或方案 | 解决的问题 | 收益 | 引入成本 | 维护风险 | 是否推荐 |
| --- | --- | --- | --- | --- | --- |
| 保留 `sphinx_rtd_theme` + 小型 CSS | 首页内容宽度、任务入口、图注/表格层次 | 低风险，移动端基础已可用；避免主题迁移 | 小；`conf.py` + 1 个 CSS | CSS 选择器随主题版本小幅变化 | **P0 推荐**；CSS 控制在必要规则内 |
| 切换 Furo | 更现代的排版、dark mode | 视觉改善 | 中；全站截图和导航回归 | 与现有品牌/布局再适配，不能解决内容空洞 | **暂不推荐** |
| 切换 PyData Sphinx Theme | 更强顶栏、版本/生态导航 | 适合大型多产品文档 | 中至大 | 对本项目规模过重；容易模仿 ASE 外观而忽略内容 | **不推荐当前采用** |
| `sphinx-design` | 首页 card/grid、tab、dropdown | 可以清晰展示 4–6 个任务入口 | 小至中；新增依赖与语法 | 容易过度组件化 | **P1 条件推荐**；原生 RST + CSS 不够时再加 |
| `sphinx-copybutton` | 复制 console code 时去掉 prompt | 改善命令体验，成本低 | 小；一项依赖/config | 很低 | **P1 推荐** |
| `sphinx-gallery` | 执行 `.py`、缩略图、下载、gallery | 对稳定纯 Python 示例价值高 | 中至大；重构案例、运行时间和缓存 | 当前示例过时/慢/含外部依赖，会造成 CI 脆弱 | **P2 条件推荐**；先稳定 3–5 个脚本 |
| `nbsphinx` | 直接渲染现有 notebook | 快速展示输出 | 中 | 保存状态、旧 import、执行顺序、环境和大输出难管 | **当前不推荐** |
| `myst-nb` | notebook + MyST authoring | 长期交互式教程能力 | 中至大 | 同时引入格式迁移和 notebook 执行复杂度 | **当前不推荐** |
| `sphinxcontrib-mermaid` | 文本维护流程图 | 修改 topology 方便 | 小至中；新增 extension/前端渲染 | 版本与渲染差异；当前图数量很少 | **暂不引入**；首批使用审核后的 SVG |
| `sphinx.ext.intersphinx` | 链接 ASE、Python、NumPy 等官方 API | 避免重复解释外部对象 | 小；built-in extension + inventory config | 外部 inventory 网络波动 | **P1 推荐**；可在 CI 缓存/容错 |
| `sphinx.ext.viewcode` | 从 API 跳到源码 | 高级用户定位实现更快 | 很小；built-in | 可能暴露内部实现但仓库公开 | **P1 推荐** |
| `linkcheck` | 检测外部/内部链接失效 | 当前已通过，可防未来漂移 | 小 | 外网站点偶发超时 | **推荐**独立 scheduled/manual job；不要阻塞每次 Pages deploy 的瞬时网络错误 |
| `doctest` | 验证文档中的短 Python 片段 | 防止最小示例失效 | 中；需要清理随机输出和依赖 | 对打印多、文件型 workflow 不合适 | **定向推荐**；只覆盖纯函数/小片段 |
| 示例代码 CI | 验证完整 Getting Started 和 fixture-based tutorials | 直接保证用户路径可运行 | 中 | 示例数据、运行时间需要控制 | **P0 推荐**；首期 1 个 smoke，之后增至 4–6 个 |
| 图片自动生成 | 防止结构图和代码漂移 | 来源可追踪、统一风格 | 中 | Matplotlib/ASE 版本造成像素差异 | **P1 推荐**；固定输入，不要求每次像素级对比 |
| API/docstring 覆盖检查 | 防止公开入口无参数/Returns/单位 | 提升 reference 完整度 | 中 | 一次全局 gate 会产生大量历史债务 | **P2 渐进推荐**；只 gate curated public API |
| MyST Markdown | 降低 Markdown 作者门槛 | 长期可能有益 | 中至大 | 无当前内容收益、产生迁移 diff | **当前不推荐** |

推荐的最小 `conf.py` 技术增量顺序：

1. 增加 `_static` 与一个 `custom.css`；
2. 增加 `intersphinx`、`viewcode`；
3. 若采用 copybutton，再把 `sphinx-copybutton` 写入 `docs/requirements.txt`；
4. 所有 extension 在同一个 PR 中由严格 HTML build 验证；
5. 暂不添加 gallery/notebook/mermaid 依赖。

## 13. 分阶段实施路线图

工作量定义：小=约半天至 1 天；中=约 2–4 天；大=需要跨模块整理、科研数据审核或一周以上。它是 issue 拆分尺度，不是承诺工期。

| 任务 | 优先级 | 具体文件 | 依赖项 | 预计工作量 | 完成标准 |
| --- | --- | --- | --- | --- | --- |
| 重写首页定位、任务入口和名称关系 | P0 | `docs/source/index.rst`、`README.md` | 确认一句官方定位 | 小 | 首页首屏回答用途/用户/任务/入口；Getting Started 一次点击可达 |
| 明确中文叙事 + 英文 API 语言策略 | P0 | `README.md`、`reference/development.rst` | 维护者确认 | 小 | README 与 docs 各自责任写入贡献规范；不再维护重复长段 |
| 建立完整 Au(111) Getting Started | P0 | 新 `getting_started/au111_slab.rst`、`docs/examples/getting_started_au111.py` | `generate_film`、ASE/Matplotlib | 中 | 新环境无需 VASP 可运行；生成 POSCAR/PNG；断言通过；页面含全部模板段落 |
| 修正文档中的真实脚本名和调用方式 | P0 | `install.rst`、`manual/vasp.rst`、`reference/scripts.rst`、`vasp_utils/README.md` | 脚本清单与 PATH 决策 | 中 | 文档不再出现无兼容说明的 `yin_*`；每个示例命令指向存在文件 |
| 修复 surface-energy 单位/乱码并补 API 合同 | P0 | `surfenergy.py`、`api/calculate.rst` | 领域维护者核对单位和 factor | 小 | 页面正确显示单位、参数、返回值和异常；转换有小型检查 |
| 重写 Examples index，标记案例状态 | P0 | `user_guide/examples.rst` 或新 `gallery/index.rst` | tutorial 清单 | 小 | 每个入口可点击，标记 verified/conditional/HPC-only；不链接已知失效 notebook |
| 扩大 docs workflow path filter，并加 Getting Started smoke | P0 | `.github/workflows/docs.yml`、`docs/examples/` | 稳定示例脚本 | 中 | 脚本目录变更触发 CI；HTML build + smoke 通过；CI 不调用 VASP/SLURM |
| 建立 surface-energy tutorial | P1 | `tutorials/surface_energy.rst`、小型 fixture、绘图脚本 | 单位修复、数据授权/审核 | 中 | 公式/API/手算一致；有输入树、结果表、图和 provenance |
| 建立 OUTCAR batch tutorial | P1 | `tutorials/outcar_batch.rst`、`docs/examples/outcar_summary.py`、小 fixture | 确认稳定 parser | 中 | 固定 fixture 输出确定 CSV/RST table；缺失文件有明确错误 |
| 建立 strained-film tutorial | P1 | `tutorials/strained_films.rst`、示例脚本、图片 | 基准结构审核 | 中 | 5–7 个结构与 manifest；factor=1.0/单调性检查通过 |
| 整合 SLURM workflow 与 dry-run tutorial | P1 | `workflows/slurm.rst`、`tutorials/slurm_dry_run.rst` | 目标集群补验 | 中 | 三种 mode、生成文件、提交副作用、恢复路径均有证据；dry-run 可在无 SLURM 环境检查 |
| 建立 n2p2 dataset tutorial | P1 | `tutorials/n2p2_dataset.rst`、小 OUTCAR fixture | 数据发布确认 | 中 | `input.data` 结构数、元素和 block 可自动核对 |
| 给教程入口补 docstring 与反向链接 | P1 | `api/*.rst`、15–25 个目标函数 docstring | Tutorial URL 稳定 | 中 | 每个核心 tutorial 的 Related API 全部可点；每个目标 API 有 Related tutorials |
| 建立 4–6 张可追踪图片 | P1 | `_static/images/generated/`、`docs/scripts/` | 统一视觉参数 | 中 | 每图有脚本、输入、caption、alt；无纯白/空图 |
| 扩展 Troubleshooting | P1 | `troubleshooting/index.rst` | 真实错误和集群信息 | 中 | 至少覆盖 install、PATH、optional dependency、POSCAR/OUTCAR、SLURM、non-convergence 6 类 |
| 将旧 `docs/build/` 退出跟踪 | P1 | `.gitignore`、`docs/build/` | 维护者确认不再使用 | 小 | Git 不再跟踪生成 HTML；Pages 仍从 `_build` 成功部署；清理由维护者明确执行 |
| 增加 `intersphinx`、`viewcode`、copybutton | P1 | `conf.py`、`docs/requirements.txt` | 内容结构稳定 | 小 | 严格构建无 warning；外部对象链接和源码链接可用 |
| 建立 curated `.py` Gallery | P2 | `docs/examples/`、`gallery/`、可选 `sphinx-gallery` | 至少 3–5 个稳定快速脚本 | 大 | 全部案例在限定时长内执行；缩略图、下载、失败日志可用 |
| 评估一个 notebook prototype | P2 | 单个新 notebook、可选 `myst-nb` | 真实交互需求 | 中 | 干净 kernel 可从头执行；不依赖私有路径/集群；证明比 `.py` 页面更有价值 |
| 建立 curated API coverage gate | P2 | CI、public API allowlist | 核心 docstring 已补齐 | 中 | 新增/修改 public API 必须包含 Args/Returns/units/related tutorial；历史内部 API 不阻塞 |
| 固定 docs 上游依赖并升级 GitHub actions | P2 | `docs/requirements.txt`、workflow | 确认兼容版本/commit | 小 | clean runner 构建可重复；Node runtime warning 消失 |
| 文档版本策略 | P2 | Sphinx/Pages 配置、贡献规范 | 出现多个受支持 release 后 | 大 | 仅在确有两个并行支持版本时实施；当前不提前建设 |

## 14. 第一轮具体修改清单

第一轮目标是用有限改动形成“首页 → 可运行结果 → 两个真实进阶任务 → 精确参考”的闭环。建议按以下 12 项实施：

1. **重写 `docs/source/index.rst`**：说明 `pjvasp_package` 与 `mymetal` 的关系，加入 Getting Started、Tutorials、Workflow、Reference 四个任务入口和一张 Au(111) 结果预览。
2. **缩短 `README.md` 的正式教程内容**：保留定位、最短安装、一个 8–12 行示例、文档入口；把详细 workflow 指向 Sphinx。
3. **扩展 `docs/source/user_guide/install.rst`**：区分 Python package、可选依赖、外部 executable、shell script PATH；增加本地/HPC 验证命令和 `module avail` 提示。
4. **新增 `docs/source/getting_started/index.rst` 与 `au111_slab.rst`**，并新增唯一来源脚本 `docs/examples/getting_started_au111.py`。
5. **新增 `docs/scripts/generate_structure_images.py`**，生成首页/Getting Started 共用的 `au111_slab.png`；图片包含 alt、caption 和确定性输入。
6. **把 `user_guide/examples.rst` 改为精选入口页**：第一轮只展示 Getting Started、surface energy、OUTCAR、SLURM dry-run；其余标记“planned”而不是链接失效 notebook。
7. **新增 `tutorials/surface_energy.rst`**：使用审核后的现有 fixture，先修 `surfenergy.py` 单位乱码；提供公式、表格、验证和 provenance。
8. **新增 `tutorials/outcar_batch.rst` 与 `docs/examples/outcar_summary.py`**：使用 1–2 个最小 OUTCAR fixture，先走 `PostTime`/`PostData`/`PostData2`，不宣传当前失效的 `post_general()`。
9. **扩写 `manual/slurm.rst` 与 `reference/scripts.rst`**：复用现有 SLURM architecture 文档；以实际 `pei_*` 文件名为准，第一轮演示 dry-run，不调用 `sbatch`。
10. **为教程涉及的 API 增加精确 reference**：至少覆盖 `generate_film`、`my_read_vasp`、`my_write_vasp`、`cal_area`、`cal_surface_energy`、三个 OUTCAR parser class 和 `pei_slurm_univ_submit`。
11. **最小化视觉调整**：保留 `sphinx_rtd_theme`；新增 `_static/css/custom.css` 只处理首页内容宽度、任务入口、figure 和宽表，不更换主题、不引入前端框架。
12. **加强 `.github/workflows/docs.yml`**：扩展脚本目录 path filter；运行严格 HTML build、Getting Started smoke 和内部链接检查；外部 linkcheck 放独立 scheduled/manual job。确认部署成功后再由维护者处理 tracked `docs/build/`。

### 14.1 第一轮建议文件树

该树基于当前 RST/Sphinx 工程，只增加第一轮实际需要的文件；不创建空 `_templates` 或空 Gallery 层级。

```text
docs/
├── examples/
│   ├── getting_started_au111.py
│   ├── outcar_summary.py
│   └── data/
│       └── outcar_minimal/              # 仅可授权的小 fixture
├── scripts/
│   ├── generate_structure_images.py
│   └── plot_surface_energy_example.py
├── source/
│   ├── index.rst
│   ├── getting_started/
│   │   ├── index.rst
│   │   └── au111_slab.rst
│   ├── tutorials/
│   │   ├── index.rst
│   │   ├── surface_energy.rst
│   │   └── outcar_batch.rst
│   ├── user_guide/
│   │   ├── install.rst
│   │   ├── examples.rst
│   │   └── troubleshooting.rst
│   ├── manual/
│   │   ├── workflows.rst
│   │   ├── vasp.rst
│   │   ├── lammps.rst
│   │   ├── slurm.rst
│   │   └── n2p2.rst
│   ├── api.rst
│   ├── api/
│   │   └── ...                          # 保留现有 6 页，定向补齐
│   ├── reference/
│   │   ├── scripts.rst
│   │   ├── dependencies.rst
│   │   └── development.rst
│   └── _static/
│       ├── css/
│       │   └── custom.css
│       └── images/
│           └── generated/
│               ├── au111_slab.png
│               ├── surface_energy_geometry.svg
│               └── slurm_modes.svg
├── requirements.txt
└── _build/                               # CI/local output，不进入 Git
```

第一轮不必立即把 `manual/`、`user_guide/` 物理重命名为目标目录。先让用户可见导航和内容闭环成立，再逐页迁移路径，避免一次产生大量失效链接。

### 14.2 第一轮验证命令

```bash
python docs/examples/getting_started_au111.py --output docs/_build/example-au111
python docs/examples/outcar_summary.py docs/examples/data/outcar_minimal
python -m sphinx -b html -W --keep-going docs/source docs/_build/html
python -m sphinx -b linkcheck -W --keep-going docs/source docs/_build/linkcheck
```

目标 CentOS HPC 还应执行：

```bash
module avail
bash -n <generated-slurm-script>
```

只有在维护者明确指定测试 partition 和授权后，才执行真实 `sbatch` 验收。

### 14.3 第一轮用户体验验收

邀请一名未熟悉仓库结构的用户完成三个任务，不提供口头补充：

1. 从首页找到并生成 Au(111) POSCAR；
2. 找到 surface-energy 输入各字段的单位；
3. 找到 SLURM dry-run，说明它是否会提交作业。

记录点击路径、完成时间、卡点和用户实际搜索词。验收重点是能否完成任务，不是“页面看起来更现代”。

## 15. 验收标准

### 15.1 信息架构与新用户

- [ ] 首页首个桌面视口内说明项目用途、目标用户、四类能力和开始入口。
- [ ] 首页明确 `pjvasp_package` 是仓库/项目，`mymetal` 是核心 Python package。
- [ ] 用户从首页一次点击到完整 Getting Started，最多两次主要导航到任一核心 Tutorial。
- [ ] Getting Started 在没有 VASP、POTCAR、SLURM 和集群账号时可完成。
- [ ] Tutorial、How-to、Workflow、API 的页面模板和职责写入 Development guide。

### 15.2 教程质量

- [ ] 每个核心 Tutorial 都有目标、时间、前置、输入树、完整步骤、预期输出、意义、验证、常见错误、下一步和 Related API。
- [ ] 所有 Python/CLI 名称均由当前源码或文件存在性检查验证。
- [ ] 适合可视化的核心 Tutorial 至少有一张解释性图片，且有 caption 和有效 alt。
- [ ] 所有数值结果来自标明 provenance 的 fixture；没有伪造 VASP 数据。
- [ ] 所有下载包排除 POTCAR、私有路径、账号、集群 secret 和许可不明数据。

### 15.3 API 与 workflow

- [ ] 第一轮教程涉及的全部稳定 API 有签名、类型、单位、返回值、异常和最小示例。
- [ ] surface-energy 页面正确显示 `Å²`/`J/m²` 等经领域维护者确认的单位，不存在乱码。
- [ ] script reference 中展示的每个 `pei_*` 文件在仓库中存在。
- [ ] 安装页不声称 `pip` 会安装当前未声明的 console scripts。
- [ ] SLURM 文档明确 dry-run、提交副作用、三种 mode、生成文件、中断恢复和状态差异。

### 15.4 构建、链接与移动端

- [ ] `python -m sphinx -b html -W --keep-going ...` 在 clean environment 无 warning。
- [ ] internal `linkcheck` 无失败；外部链接失败能区分暂时网络错误和真实失效。
- [ ] Getting Started 和至少两个 fixture-based 示例在 CI 中运行。
- [ ] `vasp_utils/**`、`slurm_utils/**`、`lmp_utils/**`、`myvasp/**`、`n2p2_utils/**` 的相关变化可触发文档检查。
- [ ] 390 px 移动端没有页面级横向溢出；代码块和宽表在容器内滚动。
- [ ] 搜索 `surface energy`、`OUTCAR`、`SLURM dry run`、`POSCAR` 时，结果中出现对应 canonical task page。
- [ ] `docs/build/` 不再作为部署源或长期跟踪生成物；GitHub Pages 仍从 workflow artifact 部署。

### 15.5 可维护性

- [ ] 每张自动生成图片有源脚本、固定输入和输出路径。
- [ ] README 不复制完整 Tutorial/Workflow。
- [ ] 不复制 ASE 的文字、图片、logo 或 CSS；只采用可解释的信息架构原则。
- [ ] 新 extension 全部进入 `docs/requirements.txt` 并在 CI 安装验证。
- [ ] CI 不运行真实 VASP、LAMMPS、n2p2 training 或 `sbatch`。
- [ ] 每个路线图任务可以独立转成 issue，且有文件、依赖和完成标准。

## 16. 风险与维护建议

### 16.1 科学正确性与“可运行”不是同一件事

Sphinx build、Python script 和图片生成成功，只证明工程链条可运行，不证明 lattice parameter、surface-energy reference、EOS fit 或 NEB/GSFE 曲线具有可发表的物理意义。每个带科研数值的案例需要两类验收：

1. 文档工程维护者确认可复现、单位、输入输出和链接；
2. 领域维护者确认模型、公式、参数、结果解释和数据 provenance。

### 16.2 可选依赖与 mocked import

当前 autodoc mock 能让页面成功构建，即使 `hetbuilder`、`myvasp`、OVITO 等未安装。应在每个相关页面显式标记：

- build-only mock；
- runtime required；
- optional feature；
- external executable；
- 已在哪个平台验证。

不要把“API 页面能渲染”当作“用户环境可运行”的证据。

### 16.3 VASP 数据与许可

- 不发布 POTCAR 或包含其内容的 archive。
- 对 OUTCAR/CONTCAR 等 fixture 记录来源、是否可公开、计算目的和 checksum。
- 教程可引用用户自行准备 POTCAR 的目录位置，但不提供下载。
- 大型 OUTCAR 只保留最小、足够覆盖 parser 的可授权样本；不要复制整套约 40 MB 数据到每个案例。

### 16.4 HPC 可移植性

`module load`、partition、launcher、wall time、文件系统和 VASP executable 都是集群特定配置。推荐：

- 通用页面只说明参数语义和 discovery 方法；
- 真实 preset 单独列为“已验证环境示例”；
- 必要命令缺失时先 `module avail`；
- 所有提交教程默认 dry-run；
- 把“提交成功”“作业完成”“VASP 收敛”作为三个不同状态记录。

### 16.5 示例与源码漂移

文档示例应复用同一个 `.py` 文件完成三件事：

1. Tutorial 的 `literalinclude`；
2. 用户下载；
3. CI 执行。

不要在首页、Tutorial 和 test 中复制三份代码。研究 notebook 可以保留为开发素材，但只有经过清理、从空 kernel 可执行、无用户路径、无保存错误、运行时间可控时，才升级为正式示例。

### 16.6 主题与扩展债务

当前 RTD 主题在桌面、搜索和移动端均可用。主题更换不会自动产生输入、输出、图片、案例或交叉链接。先完成内容 P0/P1；只有出现以下明确需求时才重新评估主题：

- 侧栏无法容纳已经稳定的 IA；
- accessibility audit 发现当前主题无法通过的阻断问题；
- 需要多版本、多语言或大型生态导航；
- 小型 CSS 无法解决核心页面布局。

### 16.7 建议维护节奏

- 每个功能 PR：若更改 public API、脚本名、目录合同或输出格式，必须同步更新 canonical reference。
- 每个教程 PR：必须附严格构建、示例运行和图片 provenance。
- 每月或每次 release：运行外部 linkcheck、搜索词验收和一个移动端 smoke。
- 每 6 个月：复核 optional dependencies、HPC presets、GitHub actions 和已发布 fixture 授权。
- 只有在有稳定新增案例时才扩展 Gallery；不为占位先建空栏目。

最终建议的实施顺序是：

```text
首页定位
→ 一个完整 Getting Started
→ 修正命令/单位/安装事实
→ 两个 fixture-based Tutorials
→ Tutorial ↔ API ↔ Workflow 交叉链接
→ 图片与 CI 机制
→ 条件性 Gallery/notebook
```

这条路径在不更换主题、不运行真实 VASP、不一次性重写全部文档的前提下，能够最直接地改善“看不懂项目、找不到任务、无法判断结果是否正确”三个核心问题。
