# Stage 4 Handoff — lmp_utils 子包文档

## 完成时间

2026-08-03（主体由前序 cron commit ffcbb12 完成；本 cron 做状态同步 + figure 补强）

## 原始产出（commit ffcbb12）

### 新增文件

1. **示例脚本**：`docs/examples/lmp_utils_overview.py`（641 行，641 insertions）
   - 多面板综合图：runner 三阶段流程（stretch → Cij → gsfe）、模板对比、
     ``.mod`` 文件依赖关系、GSFE 滑移系、``sed`` 占位符替换机制
   - 覆盖 lmp_utils 全部 20 个脚本（template/* + post/* + lmp_universal/*）
   - 中文字体回退，零 glyph 警告
   - 非空白自检通过（PNG 588034 字节）

2. **RST 文档页**：`docs/source/manual/lmp_utils.rst`（312 行）
   - 整体架构 + runner 流程图（ASCII code-block）
   - `pei_lmp_run_properties` 主 runner 参数表（位置参数：pair_style/pair_coeff/python_path/mass_content/lmp_template_path）
   - 自动重试机制说明（``LMP_MAX_TRY`` / ``LMP_RETRY_SLEEP``）
   - template / post / lmp_universal 三类脚本分节说明
   - 7 处 seealso 交叉链接（lammps / vasp_workflow_bulk / slurm_utils / workflows / gsfe_models / cij_energy_fitting / biaxial_stretch）

3. **图片**：`docs/source/_static/images/generated/lmp_utils_overview.png`（588034 字节）

### 更新文件

1. `docs/source/index.rst`：注册 `manual/lmp_utils` 到 toctree + 首页图廊引用
2. `docs/scripts/generate_structure_images.py`：注册 `render_lmp`（行 224-228）
3. `.github/workflows/docs.yml`：注册 lmp_utils 触发路径（行 10）+ 示例调用（行 80）

## 本 cron 补强（figure 引用）

### 问题

原 `lmp_utils.rst` 页面虽有丰富的 seealso 交叉链接和内容结构，但**缺少
``.. figure::`` 结果图引用**（ASE 审美规范第 5 条"复杂概念配图"未落实）。
``lmp_utils_overview.png`` 已生成且在 index.rst 首页图廊引用，但 manual
页面本身未内嵌该图。与 BLOCKERS.md 记录的 neb_utils 同类问题（minor）一致。

### 修复

在 `lmp_utils.rst` "整体架构"段末尾、"lmp_universal"详细段之前插入 figure 块：

```rst
.. figure:: /_static/images/generated/lmp_utils_overview.png
   :alt: lmp_utils 工作流总览：runner 流程、模板对比、mod 依赖、GSFE 滑移系、sed 替换
   :width: 100%

   ``lmp_utils`` 子包总览图。上图展示 ``pei_lmp_run_properties`` runner 的三阶段
   流程（stretch → Cij → gsfe）、模板与 ``.mod`` 文件的依赖关系、以及 ``sed``
   占位符替换机制。完整脚本见 :doc:`../tutorials/index` 中的相关教程。
```

### 验证

- Sphinx 构建成功（build succeeded），warning 数维持在 23（与 Stage 3 基线一致，
  无 lmp_utils 相关 warning，无 regression）
- 23 个预存 warning 全部是 autodoc import 失败（n2p2/matplotlib 等可选依赖
  未在当前环境安装），与 Stage 4 产出无关——已用 git worktree 在 Stage 3
  commit (99cdfd2) 上构建对比确认

## 安全约束遵守

- 不修改 setup.py ✅
- 不递归删除 ✅
- 不运行真实 VASP/LAMMPS/n2p2/sbatch ✅（示例脚本纯绘图，无真实计算）
- 不发布 POTCAR/secrets ✅
- 不覆盖已有产出 ✅（本 cron 只在 RST 增 figure 块，未删原有内容）

## 给后续阶段的指示

Stage 4 主体已完成。Stage 5（n2p2_utils）和 Stage 6（收尾部署）也均已由前序
cron 完成（commit bad9640 / 7657857），所有产出已 push 到 origin。

剩余可选改进（非阻塞，不影响合并）：
- BLOCKERS.md 记录的 neb_utils 缺独立 figure（minor，open）—— 可仿照本 cron
  的做法补 figure 引用
- mymetal-review.md M1：outcar_summary.py / surface_energy.py 的 CLI 风格
  统一（低优先级）

## 终止条件已满足

STAGE_STATUS.md 所有阶段（Stage 0-6）均已 completed。按 CRON_PROMPT 终止条件，
无需创建下一阶段 cron。工作分支已与 origin 同步（HEAD = origin = f6fed87）。
