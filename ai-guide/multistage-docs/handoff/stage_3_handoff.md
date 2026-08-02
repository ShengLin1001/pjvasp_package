# Stage 3 Handoff — io + post 子包

## 完成时间
2026-08-02

## 新增文件
1. docs/examples/io_extxyz_and_general.py — extxyz 轨迹 + general_read/write I/O
2. docs/examples/cij_energy_fitting.py — Cij 弹性常数 energy-strain 拟合（合成数据）
3. docs/source/tutorials/io_extxyz_and_general.rst
4. docs/source/tutorials/cij_energy_fitting.rst
5. docs/source/_static/images/generated/io_extxyz_and_general.png
6. docs/source/_static/images/generated/cij_energy_fitting.png

## 修改的文件
1. docs/source/api/io.rst — 新增 extxyz_to_atomlist/general_read/general_write autofunction + 示例
2. docs/source/api/post.rst — 新增 Cij_energy/convergence/stretch automodule + 用途说明
3. docs/source/tutorials/index.rst — 注册 2 新教程
4. docs/scripts/generate_structure_images.py — 注册 2 新示例
5. .github/workflows/docs.yml — CI 注册 2 新脚本
6. mymetal/io/extxyz.py — 修复 docstring RST 格式（markdown 代码块 → RST literal block）
7. mymetal/io/general.py — 修复 docstring RST 格式（反引号 + 缩进）

## 验证
- sphinx -W：0 warning，build succeeded
- io 示例：5 帧轨迹 round-trip，DataFrame 数值列 allclose
- cij 示例：拟合 C11/C12/C44 与输入偏差 < 0.08%

## 未覆盖
- post/gsfe.py、post/neb.py 的 automodule 因 docstring RST 格式问题暂不渲染（加了 note 说明）
- post/newmain.py 的 PostData 等已有 outcar_batch tutorial
- post/hoec_energy.py、post/relax_convergence.py、post/E_in_1_2_bulk.py 未单独做示例（依赖 VASP 目录结构）
