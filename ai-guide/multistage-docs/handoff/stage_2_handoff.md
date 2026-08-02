# Stage 2 Handoff — calculate 子包

## 完成时间
2026-08-02

## 新增文件

1. `docs/examples/deformation_and_hnf.py` — 变形矩阵 + HNF 示例
2. `docs/examples/reciprocal_lattice.py` — 倒格子矢量 + RK 网格示例
3. `docs/source/tutorials/deformation_and_hnf.rst` — 变形矩阵 + HNF 教程
4. `docs/source/tutorials/reciprocal_lattice.rst` — 倒格子 + RK 教程
5. `docs/source/_static/images/generated/deformation_and_hnf.png` (69KB)
6. `docs/source/_static/images/generated/reciprocal_lattice.png` (176KB)

## 修改的文件

1. `docs/source/api/calculate.rst` — 大幅扩充：
   - 新增 cal_stretch/cal_relative_stretch autofunction
   - 新增 hermite_normal_form autofunction + 示例
   - 新增 cal_mismatch/compare_atoms/cal_stretch_lattice autofunction
   - 新增 cal_reciprocal_matrix/cal_reciprocal_matrix2/get_size_by_distance autofunction (:no-index: 避免与 automodule 重复)
   - 新增 electronic_structure.universal automodule
2. `docs/source/tutorials/index.rst` — 注册 2 个新教程
3. `docs/scripts/generate_structure_images.py` — 注册 2 个新示例
4. `.github/workflows/docs.yml` — CI 注册 2 个新脚本
5. `mymetal/calculate/calmismatch/calhetero.py` — 修复 cal_mismatch docstring 的 Title underline too short 警告

## 验证结果

- compileall：成功
- sphinx -E -b html -W --keep-going：**成功**（0 warning）
- 2 个新示例运行成功，图片非空白
- deformation: F[0,0]=1.05 验证通过
- reciprocal: 两种倒格子方法 allclose 验证通过

## 未覆盖（依赖可选包，只做 API 文档不做可运行示例）

- calmismatch/calhetero 的 filter_results（依赖 hetbuilder）
- calmechanics/stretch（依赖 hetbuilder）
- electronic_structure/plotter（依赖 pymatgen，4586 行大模块，API 已用 automodule 暴露）
- calmechanics/hoec（1186 行，API 已有 automodule）

## 下一阶段

Stage 3（io + post 子包）。
