# BACKLOG — 动态任务清单

## P0 科学正确性
- [ ] 审查文档中 strain/stress/tensor convention 是否与源码一致
- [ ] 审查文档中单位标注（eV/Å/GPa 等）是否完整正确
- [ ] 审查 surface energy / cohesive energy 归一化公式

## P1 缺失文档（implementation > documentation）
- [ ] mymetal.build.bulk.create — 11 个公开函数无 automodule 覆盖（含 vasp_create_* 系列）
- [ ] mymetal.build.bulk.gsfe — create_gsfe_model 无 automodule 覆盖
- [ ] mymetal.build.film.extrfilm — 4 个函数无 automodule 覆盖
- [ ] mymetal.build.film.findcubic — 3 个函数无 automodule 覆盖
- [ ] mymetal.build.film.findhetero — 8 个函数无 automodule 覆盖（heterostructure 构建）
- [ ] mymetal.build.film.findprim — 2 个函数无 automodule 覆盖
- [ ] mymetal.build.film.hydroxyl — 5 个函数无 automodule 覆盖（表面 passivation）
- [ ] mymetal.build.film.stretch — 8 个函数无 automodule 覆盖
- [ ] mymetal.build.workflow.gsfe — 9 个函数无 automodule 覆盖
- [ ] mymetal.build.workflow.general — 3 个函数无 automodule 覆盖
- [ ] mymetal.calculate.calmechanics.hoec — 22 函数 + 1 class，无 automodule 覆盖
- [ ] mymetal.calculate.calmechanics.stretch — 2 个函数无 automodule 覆盖
- [ ] mymetal.calculate.electronic_structure.plotter — 10 函数 + 6 classes 无覆盖
- [ ] mymetal.cr.crplotkcontactgraph — 7 函数 + 2 classes 无覆盖
- [ ] mymetal.cr.crsum3 — 10 函数 + 1 class 无覆盖
- [ ] mymetal.io.extxyz — extxyz_to_atomlist 无覆盖
- [ ] mymetal.io.general — 3 个函数无覆盖
- [ ] mymetal.ml.n2p2.calculate.cur — 3 个函数无覆盖（CUR 特征选择）
- [ ] mymetal.ml.n2p2.workflow — PeiN2p2 class 无覆盖
- [ ] mymetal.ml.n2p2.workflow_md — PeiN2p2MD + 2 函数无覆盖
- [ ] mymetal.post.general — 4 个函数无覆盖
- [ ] mymetal.universal.atom.density — cal_density 无覆盖
- [ ] mymetal.universal.atom.miller — three_index_to_four_index 无覆盖
- [ ] mymetal.universal.check.atoms — 4 个 CNA 函数无覆盖
- [ ] mymetal.universal.plot.* — 多个绘图模块函数无 automodule 覆盖
- [ ] mymetal.universal.print.print — pr/er/warn/fail 无覆盖
- [ ] mymetal.universal.search.find — 2 个函数无覆盖

## P1 docstring 缺失（46 函数）
- [ ] mymetal.build.bulk.create: 7 个函数无 docstring
- [ ] mymetal.universal.plot.general: 7 个函数无 docstring
- [ ] mymetal.cr.crsum3: 4 个函数无 docstring
- [ ] mymetal.universal.check.type: 4 个函数无 docstring
- [ ] mymetal.universal.print.print: 4 个函数无 docstring
- [ ] mymetal.build.workflow.general: 3 个函数无 docstring
- [ ] mymetal.ml.model: 3 个函数无 docstring
- [ ] mymetal.build.film.stretch: 2 个函数无 docstring
- [ ] mymetal.post.newmain: 2 个函数无 docstring
- [ ] 其他: 10 个函数无 docstring

## P2 可复现性
- [ ] 验证所有 39 个 examples 可运行
- [ ] 检查 examples 是否有 assertions
- [ ] 检查 fixture 数据完整性

## P3 可用性
- [ ] heterostructure tutorial（findhetero 代码存在但无 tutorial）
- [ ] hydroxyl/passivation tutorial（hydroxyl 代码存在但无 tutorial）
- [ ] CR (contact resonance) 模块文档
- [ ] electronic_structure plotter 文档
- [ ] 检查 pei_* 脚本 -h 输出质量

## P4 视觉与排版
- [ ] 检查首页信息层级
- [ ] 检查 navigation dead ends
- [ ] 检查 cross-link 完整性
