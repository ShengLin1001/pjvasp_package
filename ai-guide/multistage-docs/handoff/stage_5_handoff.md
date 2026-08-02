# Stage 5 Handoff — ml + slurm + cr 子包

## 完成时间
2026-08-02

## 新增文件
1. docs/examples/slurm_script_generation.py — Slurm 脚本生成 dry-run 示例
2. docs/source/tutorials/slurm_script_generation.rst
3. docs/source/api/cr.rst — CR 分析 API 页面（新增）
4. docs/source/_static/images/generated/slurm_script_generation.png

## 修改的文件
1. docs/source/api/slurm.rst — 新增构建块函数说明 + tutorial 链接
2. docs/source/api/ml.rst — 新增 generate_radial_blocks/angular_blocks/cur_select autofunction + workflow/post/confusionmatrix/plot automodule + note
3. docs/source/api.rst — toctree 注册 api/cr
4. docs/source/tutorials/index.rst — 注册 slurm_script_generation
5. docs/scripts/generate_structure_images.py — 注册 slurm 示例
6. .github/workflows/docs.yml — CI 注册 slurm 示例
7. mymetal/ml/n2p2/calculate/post.py — 修复 module docstring RST 格式

## 验证
- sphinx -W：0 warning，build succeeded
- slurm 示例：生成含 #SBATCH 行的脚本，check_wall_time 验证通过，无 sbatch
- 图片 333KB 非空白

## 未覆盖
- ml/n2p2/workflow.py 的 PeiN2p2 类因 docstring RST 问题用 note 说明（不 automodule）
- cr/crsum3.py 和 cr/crplotkcontactgraph.py 因导入失败用 note 说明
- ml/model.py 和 ml/dataset.py 依赖 torch（mock），automodule 已暴露
