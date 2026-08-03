# 多阶段文档网站升级 — vasp_utils / slurm_utils / n2p2_utils / lmp_utils

## 任务背景

前一轮 mymetal package 文档升级（Stage 0-6, 分支 `docs/mymetal-website-upgrade`）已完成，
覆盖了 `mymetal/` 可导入 Python package 的所有子包。但 `vasp_utils/`、`slurm_utils/`、
`n2p2_utils/`、`lmp_utils/` 四个目录型脚本包目前只有极简的 manual 页面（`manual/vasp.rst`、
`manual/slurm.rst`、`manual/lammps.rst`、`manual/n2p2.rst`），缺少：

- 对每个核心脚本的功能说明、参数表、示例
- 流程图、目录树图、脚本生成结果展示图
- 丰富的 dry-run 示例和预期输出
- API 级别的详细参考

本轮升级专注这四个脚本包，为每个核心脚本生成描述、图片、介绍和示例。

## 目标

对 `vasp_utils`、`slurm_utils`、`n2p2_utils`、`lmp_utils` 下的各个子包的各个核心功能
持续生成相应的描述、图片、介绍等网站中用于说明该函数特征的功能，并部署到 GitHub Pages
在线文档。审美参考 ASE 文档（https://docs.ase-lib.org/）。

## 工作分支

`docs/utils-website-upgrade`（已从 main 创建）。所有修改在此分支进行。
main 分支保持不动。最终验证通过后再 merge 回 main 或提 PR。

## 关键约束（无人值守 yolo 模式安全规则）

1. 禁止使用 `rm -rf`、`rmdir /s`、`rd /s`、`del /s`、`Remove-Item -Recurse`
   或任何递归删除命令。需要删除文件时，一次只删一个明确指定的文件路径。
2. 不得修改 `setup.py`（用户原有修改必须保留）。
3. 不得运行真实 VASP、LAMMPS、n2p2 training 或 `sbatch`。
4. 不得发布 POTCAR、用户路径、账号或集群 secret。
5. 示例脚本不得调用真实 VASP/POTCAR；用合成数据、dry-run 模式或仓库内已 tracked fixture。
6. 不大规模重建整个 `docs/`；不覆盖用户已有改动；不覆盖上一轮 mymetal 升级的产出。
7. 每个阶段完成后必须运行 Sphinx 构建确认无回归（`--keep-going`，不强制 `-W` 以免阻塞）。
8. 中文叙事，英文 API/模块/命令/文件名。
9. commit message 遵循仓库风格：`docs(...)` + 简短中文描述。
10. 不删除已有的 `docs/source/_static/images/generated/` 图片。
11. 所有交接文档、计划文件放在 `ai-guide/utils-website-upgrade/` 子目录下。

## 阶段划分

| 阶段 | 内容 | 估时 |
|------|------|------|
| Stage 0 | 分支创建 + 计划文档 + 基线构建 + cron 设置 | 本 session |
| Stage 1 | vasp_utils/vasp_universal 子包 | 1 session |
| Stage 2 | vasp_utils/vasp_workflow_bulk + neb_utils 子包 | 1 session |
| Stage 3 | slurm_utils 子包 | 1 session |
| Stage 4 | lmp_utils 子包 | 1 session |
| Stage 5 | n2p2_utils 子包 | 1 session |
| Stage 6 | 审查收尾：交叉链接、首页、全量构建、commit、push、部署 | 1 session |

## 子包覆盖范围

### Stage 1 — vasp_utils/vasp_universal（单目录工具集）

核心脚本（bash + python 混合）：

| 脚本 | 功能 |
|------|------|
| `pei_vasp_univ_sbatch` | 单个 VASP 目录的 runner：判收敛、断点续算、清理、srun、退出码 |
| `pei_vasp_univ_post` | 汇总 y_dir 计算状态 |
| `pei_vasp_univ_clean_up_full` | 清理大型输出文件 |
| `pei_vasp_univ_clean_up_small` | 清理小型输出文件 |
| `pei_vasp_univ_clean_old_slurm` | 清理旧 slurm 输出 |
| `pei_vasp_univ_clean_outcar` | 清理 OUTCAR |
| `pei_vasp_univ_cp_contcar_cartesian_poscar` | CONTCAR→POSCAR（笛卡尔坐标） |
| `pei_vasp_univ_extract_convergence` | 提取收敛信息 |
| `pei_vasp_univ_extract_energy_components` | 提取能量分量 |
| `pei_vasp_univ_find_and_change` | 查找并替换 INCAR 参数 |
| `pei_vasp_univ_increase_nbands` | 增加 NBANDS |
| `pei_vasp_univ_load_env` | 加载环境 |
| `pei_vasp_univ_monitor_error` | 监控错误关键字 |
| `pei_vasp_univ_monitor_slurm_state` | 监控 slurm 状态 |
| `pei_vasp_univ_resubmit_isif3` | 重新提交 ISIF=3 |
| `pei_vasp_univ_resubmit_isym0` | 重新提交 ISYM=0 |
| `pei_vasp_univ_transfer_normal_to_selective` | 普通→选择性动力学 |
| `pei_vasp_univ_transfer_selective_to_normal` | 选择性→普通动力学 |
| `pei_vasp_univ_get_size_by_distance.py` | 按 k 点距离计算 size |
| `pei_vasp_univ_check_phase_transition` | 检查相变 |
| `pei_vasp_univ_get_struct_infos` | 获取结构信息 |
| `pei_vasp_plot_convergence.py` | 绘制收敛曲线 |

### Stage 2 — vasp_utils/vasp_workflow_bulk + neb_utils

| 脚本 | 功能 |
|------|------|
| `pei_vasp_run_eos.py` | EOS workflow |
| `pei_vasp_run_stretch.py` | stretch workflow |
| `pei_vasp_run_neb` | NEB workflow |
| `pei_vasp_run_convergence` | convergence workflow |
| `pei_vasp_run_cij_energy` | Cij energy workflow |
| `pei_vasp_run_cohesive` | cohesive energy workflow |
| `pei_vasp_run_dos_band` | DOS/band workflow |
| `pei_vasp_run_surface_energy` | surface energy workflow |
| `pei_vasp_run_hoec_energy` | HOEC energy workflow |
| `pei_vasp_run_kpar_ncore` | KPAR/NCORE workflow |
| `pei_vasp_plot_convergence.py` | convergence 绘图 |
| `pei_vasp_plot_hoec_energy.py` | HOEC energy 绘图 |
| `pei_vasp_plot_neb.py` | NEB 绘图 |
| `pei_vasp_plot_stretch.py` | stretch 绘图 |
| `pei_vasp_univ_neb_nebbarrier.py` | NEB barrier 提取 |
| `pei_vasp_univ_neb_nebbarrier_full.py` | NEB barrier 完整提取 |
| `pei_vasp_univ_neb_plot.py` | NEB 绘图 |
| `pei_vasp_univ_neb_select.py` | NEB 路径选择 |

### Stage 3 — slurm_utils

| 脚本 | 功能 |
|------|------|
| `slurm_universal/pei_slurm_univ_submit.py` | 软件无关提交入口 |
| `slurm_universal/pei_slurm_univ_launch_retry` | launcher 重试 |
| `slurm_universal/pei_slurm_univ_monitor_error` | 错误监控 |
| `slurm_universal/pei_slurm_univ_useful_command.sh` | 命令速查 |
| `slurm_universal/pei_slurm_univ_sbatch_retry` | sbatch 重试 |
| `slurm_vasp/pei_slurm_univ_vasp_monitor` | VASP 监控启动器 |
| `slurm_n2p2/pei_slurm_univ_n2p2_train` | n2p2 训练提交 |
| `slurm_n2p2/pei_slurm_univ_n2p2_scaling` | n2p2 scaling 提交 |

### Stage 4 — lmp_utils

| 文件 | 功能 |
|------|------|
| `template/gsfe_model.py` | GSFE 模型构建 |
| `template/stretch_template.in` | stretch 模板 |
| `template/gsfe_template.in` | GSFE 模板 |
| `template/Cij_energy_template.in` | Cij energy 模板 |
| `template/*.mod` | 各模块文件 |
| `post/Cij_energy.py` | Cij energy 后处理 |
| `post/gsfe.py` | GSFE 后处理 |
| `post/stretch.py` | stretch 后处理 |
| `lmp_universal/pei_lmp_run_properties` | properties workflow |

### Stage 5 — n2p2_utils

| 脚本 | 功能 |
|------|------|
| `data_input.py` | 数据输入 |
| `data_read.py` | 数据读取 |
| `sfs_gen_basic_SF.py` | 对称函数生成 |
| `active_sf_0_collect_sf.py` | active learning SF 收集 |
| `active_sf_0_select_feat.py` | active learning 特征选择 |
| `active_sf_0_sub_cal.py` | active learning 子计算 |
| `n2p2_universal/pei_n2p2_univ_run` | n2p2 workflow runner |
| `n2p2_universal/pei_n2p2_univ_load_env` | 环境加载 |
| `n2p2_universal/pei_n2p2_univ_clean_train` | 训练目录清理 |

### Stage 6 — 审查收尾 + 部署

- 全站交叉链接检查
- 首页 index.rst 更新（加新 tutorial 入口、新图廊）
- 全量 Sphinx 构建
- compileall（如有 python 变更）
- git commit（分功能边界）
- git push origin docs/utils-website-upgrade
- 检查 GitHub Actions 构建（注意：CI 只在 push main 时自动部署）
- 建议：merge 到 main 或修改 workflow 触发条件

## 每阶段通用产出

每个子包阶段必须产出：

1. **可运行示例脚本**（`docs/examples/<topic>.py`）：无需真实 VASP/LAMMPS/n2p2/sbatch，
   含确定性断言，生成流程图/目录树图/脚本内容展示图。
2. **RST 文档页**（`docs/source/manual/<topic>.rst` 或 `docs/source/tutorials/<topic>.rst`）：
   包含目标、前置条件、脚本说明、参数表、dry-run 示例、预期输出、结果图、
   验证方法、常见错误、相关交叉链接。
3. **生成的图片**（`docs/source/_static/images/generated/<topic>.png`）：
   非空白，有 caption 和 alt text。图片类型包括：
   - 流程图（脚本执行流程）
   - 目录树图（生成文件布局）
   - 脚本内容展示（渲染的 shell 脚本文本）
   - dry-run 输出展示
   - 架构图（组件关系）
4. **manual 页面更新**：扩展 `docs/source/manual/vasp.rst`、`slurm.rst` 等，
   或新建子页面。
5. **图片生成脚本注册**：在 `docs/scripts/generate_structure_images.py` 注册
   新示例的 render 函数。
6. **CI workflow 更新**：在 `.github/workflows/docs.yml` 的 examples 步骤注册
   新示例脚本。
7. **导航更新**：在相关 `index.rst` 或 `manual/` 的 toctree 注册新页面。

## ASE 审美规范（已确立，沿用上一轮）

参见 `ai-guide/multistage-docs/ASE_AESTHETIC_STUDY.md`。关键点：

- 白底 `#ffffff`，正文色 `rgb(34,40,50)`
- 代码块背景 `#f3f4f5` + `1px solid #d1d5da` 边框
- 链接色调暗，橙色点缀
- 每页开头一句话说清"这个脚本是什么"
- 紧跟最小可运行示例
- 关键参数用表格说明
- 复杂概念配流程图/目录树图
- 交叉链接密集
- 末尾"See also"列出相关模块

## 子agent分工模式

每阶段可并行调用最多 3 个 subagent：

- **Agent A（示例编写者）**：读脚本源码 → 写 `docs/examples/<topic>.py` → 运行验证 →
  产出图片路径。
- **Agent B（文档撰写者）**：根据 Agent A 的脚本和图片，写 RST 页面、manual 页面、
  更新 toctree 和 generate_structure_images.py。
- **Agent C（审查者）**：根据 ASE 审美规范和本计划检查 Agent B 的产出，反馈意见。

串行依赖：Agent A 先行 → Agent B 跟进 → Agent C 审查。但不同子模块的 A 可并行。

## 验证清单（每阶段结束必须执行）

```bash
# 在仓库根目录
'C:/Users/louis/mysoft/env/pyenv/codex/Scripts/python.exe' -m sphinx -E -b html --keep-going docs/source docs/_build/html
'C:/Users/louis/mysoft/env/pyenv/codex/Scripts/python.exe' docs/scripts/generate_structure_images.py
# 检查图片非空白
ls -la docs/source/_static/images/generated/
```

## Python 环境

- Sphinx 构建：`C:/Users/louis/mysoft/env/pyenv/codex/Scripts/python.exe` (3.10.7)
  - 已安装 sphinx 8.1.3, ase 3.28, mymetal, myvasp importable
- 示例脚本运行：同一环境

## 状态追踪

每阶段开始/结束时更新 `ai-guide/utils-website-upgrade/STAGE_STATUS.md`。
每阶段结束时写 `ai-guide/utils-website-upgrade/handoff/stage_N_handoff.md`。
