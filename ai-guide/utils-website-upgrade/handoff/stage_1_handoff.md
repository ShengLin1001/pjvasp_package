# Stage 1 Handoff — vasp_utils/vasp_universal 子包

## 完成时间
2026-08-03

## 完成内容

### 新增文件

1. **示例脚本**：`docs/examples/vasp_universal_overview.py` (23661 bytes)
   - 5 面板综合图片：runner 流程图、目录扫描、清理对比、INCAR 操作、退出码约定
   - 支持 `--output` 参数
   - 中文字体回退（Microsoft YaHei + DejaVu Sans）
   - 非空白图片自检（pixel_deviation > 0.002）
   - 不运行真实 VASP/sbatch

2. **RST 文档页**：`docs/source/manual/vasp_universal.rst` (11113 bytes)
   - 覆盖 22 个脚本的详细说明
   - 核心内容：
     - pei_vasp_univ_sbatch 流程图和退出码约定
     - pei_vasp_univ_post 产出文件表
     - 四个清理脚本对比表
     - INCAR 操作（find_and_change, increase_nbands）
     - POSCAR 操作（transfer_*, cp_contcar）
     - 收敛数据提取与绘图
     - 能量分量提取
     - 监控脚本
     - 续算脚本
     - 环境与结构信息
   - 交叉链接到 vasp.rst, slurm.rst, workflows.rst, scripts.rst, kpoints_sampling tutorial

3. **图片**：`docs/source/_static/images/generated/vasp_universal_overview.png` (583676 bytes)
   - 5 面板，2310×1749 像素
   - 所有面板验证非空白

### 更新文件

1. `docs/source/index.rst`：在 toctree 注册 `manual/vasp_universal`
2. `docs/source/manual/vasp.rst`：添加到 vasp_universal 的交叉链接
3. `docs/scripts/generate_structure_images.py`：注册 `render_overview` 函数
4. `.github/workflows/docs.yml`：在 examples 步骤注册新示例脚本

### 验证

- Sphinx 构建：成功，0 warnings，0 errors
- generate_structure_images.py：成功生成所有图片
- 图片非空白：570KB，pixel_deviation=0.087

## 覆盖的脚本（22个）

| 脚本 | 类型 | 文档中的说明 |
|------|------|-------------|
| pei_vasp_univ_sbatch | bash | 流程图 + 退出码 + 完成标记检测 |
| pei_vasp_univ_post | bash | 产出文件表 + 用法 |
| pei_vasp_univ_clean_up_full | bash | 清理对比表 |
| pei_vasp_univ_clean_up_small | bash | 清理对比表 |
| pei_vasp_univ_clean_old_slurm | bash | 清理对比表 |
| pei_vasp_univ_clean_outcar | bash | 清理对比表 |
| pei_vasp_univ_find_and_change | bash | 参数规则 + 用法 |
| pei_vasp_univ_increase_nbands | bash | 用法 |
| pei_vasp_univ_transfer_normal_to_selective | bash | 功能说明 |
| pei_vasp_univ_transfer_selective_to_normal | bash | 功能说明 |
| pei_vasp_univ_cp_contcar_cartesian_poscar | bash | 功能说明 |
| pei_vasp_univ_extract_convergence | bash | 用法 + 产出 |
| pei_vasp_plot_convergence.py | python | 用法 |
| pei_vasp_univ_extract_energy_components | bash | 用法 |
| pei_vasp_univ_monitor_error | bash | 用法 + 参数 |
| pei_vasp_univ_monitor_slurm_state | bash | 功能说明 |
| pei_vasp_univ_resubmit_isif3 | bash | 功能说明 |
| pei_vasp_univ_resubmit_isym0 | bash | 功能说明 |
| pei_vasp_univ_load_env | bash | source 说明 |
| pei_vasp_univ_get_struct_infos | python | 用法 + 参数 |
| pei_vasp_univ_check_phase_transition | python | 退出码表 |
| pei_vasp_univ_get_size_by_distance.py | python | 用法 + 参数 |

## 给 Stage 2 的指示

Stage 2 覆盖 `vasp_utils/vasp_workflow_bulk` + `neb_utils`。这些是目录型 workflow
脚本（EOS, stretch, NEB, convergence, Cij, cohesive, DOS/band, surface energy,
HOEC, KPAR/NCORE）和绘图脚本。

建议：
- 用 subagent 并行处理 workflow 脚本和绘图脚本
- 图片可以展示目录树结构（workflow 生成的 y_dir 布局）和示例输出
- 参考 Stage 1 的 RST 结构和 ASE 审美规范

## 安全约束状态

- 未修改 setup.py ✓
- 未运行真实 VASP/sbatch ✓
- 未使用递归删除 ✓
- 未覆盖已有产出 ✓
- 中文叙事，英文 API/命令名 ✓
