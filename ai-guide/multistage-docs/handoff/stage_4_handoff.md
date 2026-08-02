# Stage 4 Handoff — universal 子包

## 完成时间
2026-08-02

## 新增文件
1. docs/examples/miller_index_and_density.py
2. docs/examples/periodic_table_and_arkel.py
3. docs/source/tutorials/miller_index_and_density.rst
4. docs/source/tutorials/periodic_table_and_arkel.rst
5. docs/source/_static/images/generated/miller_index_and_density.png
6. docs/source/_static/images/generated/periodic_table_heatmap.png
7. docs/source/_static/images/generated/van_arkel_triangle.png

## 修改的文件
1. docs/source/api/universal.rst — 新增 miller/density/rotate/delatom/check.type/check.atoms(CNA)/patterntrans/matrix.adjust/plotting(plotting+workflow) automodule + autofunction
2. docs/source/tutorials/index.rst — 注册 2 新教程
3. docs/scripts/generate_structure_images.py — 注册 3 新图像
4. .github/workflows/docs.yml — CI 注册 2 新脚本
5. mymetal/universal/plot/plotting.py — 修复 van_arkel_triangle docstring RST 格式

## 验证
- sphinx -W：0 warning，build succeeded
- miller: 3↔4 转换 round-trip 一致，密度 Cu 8.97/Fe 7.85 g/cm³ 正确
- periodic_table: 215KB 非空白
- van_arkel: 172KB 非空白

## 未覆盖
- check/atoms.py 的 CNA 函数依赖 ovito（可选），API 已用 autofunction 暴露但不做可运行示例
- plot/workflow.py 的 my_plot_neb* 因 docstring RST 问题用 exclude-members 排除
- plot/general.py, plot/plot.py, plot/n2p2.py 等大模块未单独做示例
