# Post-Processing 模块网站扩展计划 (2026-08-21)

## 目标

从 `mymetal/post/` 提炼有用且重要的后处理模块，为每个模块在网站中创建
图文并茂的 tutorial 页面，补充真实计算数据图表，扩充 example 案例。

## 扩展范围（按优先级）

### Priority 1: 核心后处理模块（有历史数据，可生成真实图表）

| 模块 | 函数 | 服务器数据路径 | 扩展内容 |
|------|------|----------------|----------|
| `hoec_energy` | `post_hoec_energy` | `construct_dataset/calculate/A11-2/y_hoec_energy/` | 高阶弹性常数 2/3/4 阶，能量-应变拟合曲线，窗口扫描 |
| `Cij_energy` | `post_lammps_Cij_energy` | `pj-test-properties-gold/y_Cij_energy/` | 二阶弹性常数 energy-strain 方法 |
| `gsfe` | `post_gsfe` | `pj-test-properties-gold/y_gsfe/` | 广义层错能 γ-surface |
| `stretch` | `post_stretch`/`post_lammps_stretch` | `construct_dataset/calculate/A11-1/y_stretch/` | 单轴拉伸能量-应变曲线 |
| `convergence` | `post_convergence` | `construct_dataset/calculate/` 收敛测试目录 | ENCUT/k-point 收敛 |
| `relax_convergence` | `my_univ_post_convergence` | VASP OUTCAR ionic relaxation | 离子弛豫收敛轨迹 |
| `neb` | `post_neb` | NEB 计算目录 | NEB 路径 spline 拟合 |
| `E_in_1_2_bulk` | `post_E_in_1_2_bulk` | 双轴变形 E_in_1/2 | 2D 等高线 + profile |
| `kpar_ncore` | `post_kpar_ncore` | KPAR/NCORE benchmark | 并行性能基准 |

### Priority 2: 辅助模块

| 模块 | 内容 |
|------|------|
| `general` | `my_sort`, `my_ployfit`, `get_structure_info` — 后处理通用工具 |
| `newmain` | OUTCAR parser 类族 — 已有 outcar_batch tutorial，补充高级用法 |

## example 命名规范

旧命名 `test-*` → 新规范：`<domain>-<topic>`，全小写连字符分隔。

| 旧名 | 新名 | 说明 |
|------|------|------|
| `test-post` | `post-outcar-parse` | OUTCAR 解析 |
| `test-n2p2-sfparams` | `n2p2-symfunc-params` | n2p2 对称函数参数 |
| `test-surface-energy` | `surface-energy-fitting` | 表面能拟合 |
| `test-hydroxylated` | `surface-hydroxylation` | 羟基化表面 |
| `test-hydroxylated-custom` | `surface-hydroxylation-custom` | 自定义羟基化 |
| `test-hetbuilder-fixatom` | `heterostructure-builder` | 异质结构构建 |
| `test-stretch` | `stretch-film-generation` | 拉伸薄膜生成 |
| `test-generate-bulk` | `bulk-structure-generation` | 体相结构生成 |

新增 example（docs/examples/ 下的 VASP-free 脚本）：
- `post_hoec_energy_demo.py` — HOEC 拟合演示（合成数据）
- `post_cij_energy_demo.py` — Cij 拟合演示（合成数据）
- `post_gsfe_demo.py` — GSFE γ-surface 演示（合成数据）
- `post_stretch_demo.py` — 拉伸曲线演示（合成数据）
- `post_convergence_demo.py` — 收敛测试演示（合成数据）
- `post_relax_convergence_demo.py` — 离子弛豫收敛演示（合成数据）
- `post_kpar_ncore_demo.py` — KPAR/NCORE 基准演示（合成数据）
- `post_E_in_1_2_bulk_demo.py` — E_in_1/2 等高线演示（合成数据）

## 图表风格约束

参考 MASTER_PLAN.md (website-enhancement-2026) 和 p-plot-figure skill：
- 用 `my_plot()` 或 `general_set_all_rcParams()` 入口
- 画完图例后 `general_modify_legend`
- `fig.savefig(path, bbox_inches='tight')`，不叠 `tight_layout()`
- 网站图片用本地静态 PNG（`docs/source/_static/images/generated/`）
- 生成时正确栅格化中文，不依赖浏览器补字体
- 数值不直接显示 `np.float64(...)`，先转 Python `float`/`int`
- 流程图方框按实际多行文字扩容，文字不越界

## 执行流程

1. 本 session 做编排 + 核心 example 脚本 + RST 骨架
2. 长任务（SSH 服务器数据收集 + 真实图表生成）用 cron 链或 subagent
3. 每个模块：example → RST → 真实数据图 → build gate → commit

## 服务器约束

- SSH zcm6 (scg6928@ZC-M6)，登录节点直接运行
- 绝不 sbatch/srun
- 服务器源目录只读；写操作在 `~/mywork/temp/<task>/`
- 不泄露 POTCAR/密码/集群配置
