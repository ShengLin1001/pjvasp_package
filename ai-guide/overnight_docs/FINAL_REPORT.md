# FINAL REPORT — 无人值守文档深度建设

## Executive summary

今晚对 pjvasp_package 仓库及其 Sphinx/GitHub Pages 文档网站进行了系统性的深度审计与质量提升。该仓库文档已经相当成熟（96 RST 文件、39 个 examples、87.7% docstring 覆盖率），本轮工作聚焦于：(1) 补齐 API autodoc 覆盖缺口；(2) 新增可运行教程；(3) 补充缺失 docstring；(4) 修复科学性错误；(5) 修复代码 bug。

所有验证通过：compileall 0 errors，Sphinx strict build 0 warnings，linkcheck 0 broken。

## Major improvements

### API autodoc 覆盖扩展
- 新增 `mymetal.calculate.calmechanics.hoec` automodule（22 函数 + 1 class）
- 新增 `mymetal.build.workflow.gsfe` automodule（9 函数）
- 新增 `mymetal.build.workflow.general` automodule（3 函数）
- 新增 `mymetal.ml.n2p2.workflow_md` automodule（2 函数 + PeiN2p2MD class）
- 新增 `mymetal.post.general` automodule（4 函数）
- 新增 `mymetal.universal.print.print` automodule（5 函数）
- 新增 `mymetal.universal.print.printafter` automodule（4 函数）
- 新增 `mymetal.universal.search.find` automodule（2 函数）
- 新增 `mymetal.universal.plot.oldplotdos` automodule（5 函数）

### docstring 覆盖率提升
- 从 87.7%（327/373）提升到 94.9%（354/373）
- 补充了 27 个函数的 docstring，涉及 build.bulk.create、universal.plot.general、
  universal.check.type、universal.print.print、build.workflow.general、
  build.film.stretch、universal.atom.density 等模块

## New tutorials

### surface_passivation（表面钝化与悬挂键检测）
- 新增 `docs/examples/surface_passivation.py`（可运行，含 assertions）
- 新增 `docs/source/tutorials/surface_passivation.rst`
- 演示 `add_hydro_atoms`、`get_neighbors`、`find_matching_atom_in_bulk`、
  `find_unsaturated_neighbors`
- 已加入 CI workflow

### hoec_modes（高阶弹性常数模式表）
- 新增 `docs/examples/hoec_modes.py`（可运行，含 assertions）
- 新增 `docs/source/tutorials/hoec_modes.rst`
- 演示 `selftest_hoec`、`get_model`、`get_hoec_modes`、`get_strain_list`、
  `get_deformation_gradient`、`check_symmetry`
- 验证 Wang-Li Table I 一致性
- 已加入 CI workflow

## API/documentation fixes

### 科学性错误修正（P0）
1. **双轴应变公式错误**（`cij_energy_fitting.rst`）：
   文档写 `U/V0 = 0.5*(C11+C12)*eta²`，正确公式为 `(C11+C12)*eta²`。
   脚本使用正确公式，仅 RST 文档有误。
2. **Cauchy-Green 张量标签互换**（`strain.py`）：
   `C = F^T F` 被标为 "left Cauchy-Green"，正确为 "right Cauchy-Green"；
   `B = F F^T` 被标为 "right Cauchy-Green"，正确为 "left Cauchy-Green"。
   数学正确，仅标签互换。
3. **HCP Miller 指数 round-trip 描述**（`miller_index_and_density.rst`）：
   文档称 `[2,-1,0,0] → [1,0,0]`，实际代码返回 `[3,0,0]`（同方向不同大小）。
   已修正描述。

### RST 格式修复
- 修复 `workflow_md.py` docstring 中 6 个导致 Sphinx 构建失败的 RST inline markup 问题

## Tests and reproducibility

- 新增 example `surface_passivation.py`：build Au(111) slab → H passivation →
  dangling bond detection on Si，含 7 个 assertions
- 新增 example `hoec_modes.py`：run self-test → inspect mode tables →
  verify deformation gradient → symmetry detection，含 6 个 assertions
- 两个 example 均加入 `.github/workflows/docs.yml` CI pipeline

## Bugs found

### hydroxyl.py — passivate_surface_custom 真值检查 bug
- 文件：`mymetal/build/film/hydroxyl.py`，line 415
- 问题：`if idx_bulk:` 当 `idx_bulk = -1`（未找到匹配）时为 truthy，
  导致后续 `get_neighbors(-1, ...)` 抛 TypeError
- 修复：`if idx_bulk:` → `if idx_bulk >= 0:`
- 验证：修复前 passivate_surface_custom 在 Si slab 上崩溃，修复后正常运行

## Validation

| 检查 | 命令 | 结果 |
|------|------|------|
| compileall | `python -m compileall mymetal` | 0 errors |
| Sphinx build | `python -m sphinx -E -b html -W --keep-going docs/source docs/_build/html` | 0 warnings, build succeeded |
| linkcheck | `python -m sphinx -b linkcheck -W --keep-going docs/source docs/_build/linkcheck` | 0 broken |
| pytest | `python -m pytest tests/ -q --tb=no` | 51 passed, 16 failed (Windows path only) |
| examples | 运行 surface_passivation, hoec_modes, au111, cij, miller, gsfe, atom_manipulation, surface_energy | 全部 PASS |

## Commits

1. `cdfef38` docs(api): 扩展 API autodoc 覆盖并修复 workflow_md RST 警告
2. `4c0f743` docs(tutorial): 新增表面钝化教程并修复 hydroxyl bug
3. `cd6d976` docs(api): 补充 create/print/check 模块缺失 docstring
4. `6cd67c5` docs(tutorial): 新增 HOEC 弹性常数教程和可运行示例
5. `2991fda` docs(api): 补充 plot.general 模块缺失 docstring
6. `18c183d` docs(api): 补充 workflow.general 和 stretch 模块缺失 docstring
7. `1b0b038` docs(api): 补充 cal_density docstring
8. `edb6ee8` fix(docs): 修正科学性错误：双轴应变公式、Cauchy-Green 标签、Miller 指数 round-trip

## Remaining issues

1. **test_slurm_submit.py 16 failures**：Windows 上 `Path("/work").is_absolute()` 返回 False
   （需要盘符如 `C:/work`）。这些测试在 Linux CI 上通过。这是平台差异，不是代码 bug。
2. **heterostructure tutorial 缺失**：`findhetero` 模块需要可选的 `hetbuilder` 包，
   无法在最小安装下创建可运行 example。已记录在 tutorials/index.rst 的"后续案例"。
3. **19 个函数仍缺 docstring**：主要在 cr.crsum3（4个）、ml.model（3个）、post.oldmain（1个）
   等较小或专门模块。覆盖率已达 94.9%。

## Decisions requiring author review

1. **heterostructure tutorial**：是否要写一个需要 hetbuilder 的教程？
   建议：等 hetbuilder 在 CI 中可用后再加。
2. **workflow_md.py docstring 中移除 `` 内联标记**：为修复 RST 警告，
   移除了一些 ``nnp-dataset`` 等内联字面量标记。是否保留可读性？
3. **file_name 函数 known bug**：`build.film.stretch.file_name` 返回函数对象
   而非字符串（name-shadowing）。已在 docstring 中标记 deprecated。是否删除？
