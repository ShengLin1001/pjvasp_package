# Stage 1 Handoff — build 子包

## 完成时间
2026-08-02

## 新增文件

1. `docs/examples/gsfe_models.py` — GSFE 模型构建示例（FCC_111/HCP_basal/HCP_prism1w）
2. `docs/examples/cubic_cell_and_stretch.py` — 正交 cell 查找 + 单轴拉伸示例
3. `docs/source/tutorials/gsfe_models.rst` — GSFE 模型教程页
4. `docs/source/tutorials/cubic_cell_and_stretch.rst` — 正交 cell + 拉伸教程页
5. `docs/source/_static/images/generated/gsfe_models.png` — GSFE 三模型侧视图（273KB，非空白）
6. `docs/source/_static/images/generated/cubic_cell_and_stretch.png` — 正交+拉伸对比图（225KB，非空白）

## 修改的文件

1. `docs/source/api/build.rst` — 大幅扩充"其他结构入口"：
   - 新增 create_fcc_111、create_hcp_basal、create_hcp_prism1、create_gsfe_model 的 autofunction + 用途/参数/示例
   - 新增 find_cubic、find_optimal_cell_shape、stretch_along_direction_to_cell 的 autofunction + 用途/参数/示例
   - 每个函数都有 Related tutorials 交叉链接
2. `docs/source/tutorials/index.rst` — toctree 注册 gsfe_models、cubic_cell_and_stretch
3. `docs/scripts/generate_structure_images.py` — 注册两个新示例的 render 函数
4. `.github/workflows/docs.yml` — CI examples 步骤注册两个新脚本

## 验证结果

- `compileall mymetal`：成功
- `sphinx -E -b html -W --keep-going`：**成功**（build succeeded，0 warning）
- `generate_structure_images.py`：成功生成两张新图
- gsfe_models.py 运行：18/12/24 原子，断言通过
- cubic_cell_and_stretch.py 运行：面积守恒 8.9236→8.9236 Å²，断言通过

## subagent 使用

本次用了 2 个并行 leaf subagent：
- Agent A：写 gsfe_models.py（含运行验证）
- Agent B：写 cubic_cell_and_stretch.py（含运行验证）
两者均成功。函数名在集成到 generate_structure_images.py 时需核对（cubic 的
render_figure 签名是 (prim, ortho, stretched, path)，不是 (models, path)）。

## 未完成项

- hydroxyl.py 的 passivate_surface_custom / add_hydro_atoms 未做示例（函数较复杂，
  需要 slab + 邻居分析，留到后续阶段或单独 tutorial）
- workflow/hoec.py 和 workflow/kpar_ncore.py 的目录生成器未做示例（需要 y_full_relax
  目录和 myvasp 依赖，不适合无 VASP 示例）。API 页面已有 automodule 暴露。
- findhetero.py 依赖 hetbuilder（可选），未做示例。

## 下一阶段建议

Stage 2（calculate 子包）可以开始。重点：
- calmath/matrix.py 的 hermite_normal_form
- calmechanics/deformation.py 的 cal_deform_matrix
- calmechanics/strain.py 的 cal_strain_matrix / cal_von_mises_strain（已有 strain_deformation tutorial）
- calmismatch/calhetero.py 的 cal_mismatch
- calqm/kpoints.py 的 get_kpoints_by_size / cal_reciprocal_matrix（已有 kpoints_sampling tutorial）
- electronic_structure/plotter.py 的 plot_brillouin_zone_from_kpath（4586 行，大模块）
- material_science/schmid.py 的 cal_fcc_schmid_factors（已有 schmid_factor tutorial）

## 注意事项（给下一阶段）

- 工作期间发现其他 session 会切换分支到 feat/mymetal-example-expansion。
  cron job 启动时必须确认在 docs/mymetal-website-upgrade 分支，必要时
  git stash + checkout + stash pop 恢复工作。
- codex venv 路径：'C:/Users/louis/mysoft/env/pyenv/codex/Scripts/python.exe'
- RST 中 inline literal 后面紧跟全角括号 `（` 会触发 "Inline literal start-string
  without end-string" 警告。解决：literal 和括号之间加空格。
