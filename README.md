# pjvasp_package

`pjvasp_package` 是一个面向计算材料科学的综合工作区，核心是 `mymetal` Python 包，外围集成了 VASP、LAMMPS、SLURM/PBS、n2p2 相关脚本、示例数据和 Sphinx 文档。

这个项目不是单一用途的小包，而更像一个研究用工具箱，主要服务于以下几类任务：

- 晶体、薄膜、表面和异质结构的构建
- 面向 VASP 的输入输出处理与后处理
- 体相、表面、拉伸、NEB、GSFE、解理等工作流脚本
- 从 VASP 结果提取机器学习势训练数据
- 针对高通量批量任务的脚本化提交、监控、清理和归档

## 1. 项目概览

仓库当前可以分成四个层次：

1. `mymetal/`
   这是核心 Python 包，负责结构构建、计算辅助、I/O、后处理、绘图、机器学习数据准备等。
2. `vasp_utils/`、`myvasp/`
   这是偏工作流和命令行脚本层的工具，主要围绕 VASP 任务准备、提交和结果整理。
3. `slurm_utils/`、`lmp_utils/`
   这是集群调度与 LAMMPS 输入模板层。
4. `docs/`、`mymetal/example/`、`matlab_code/`
   这是文档、示例和补充分析脚本层。

如果你只关心 Python 接口，优先看 `mymetal/`。

如果你只关心实际跑 VASP 批处理任务，优先看 `vasp_utils/`、`myvasp/`、`slurm_utils/`。

## 2. 仓库适合什么人

这个仓库更适合以下用户：

- 已经在使用 ASE / VASP / LAMMPS 的材料计算研究者
- 需要批量生成结构和批量整理结果的使用者
- 需要把 VASP 输出转成 n2p2 数据集的用户
- 想把“建模 + 提交 + 后处理 + 绘图”放到同一套工作区中的个人研究项目

如果你的目标只是一个干净、最小依赖的公开发行 Python 包，这个仓库目前不是那个定位。

## 3. 核心目录说明

### `mymetal/`

项目主包。按职责大致分为：

- `build/`
  结构生成与几何构造。包括 bulk / film / heterostructure / hydroxylation / stretch 等。
- `calculate/`
  若干材料学计算辅助函数，例如表面能、应变、应力、错配、k 点等。
- `io/`
  VASP 文件读写和部分后处理结果写出。
- `post/`
  面向 VASP 目录批处理的后处理逻辑，例如批量读 `OUTCAR`、统计参数、提取 warning。
- `universal/`
  通用底层工具，包含 atom、plot、search、check、matrix、data 等。
- `ml/`
  机器学习相关工具，尤其是 n2p2 数据集读写与对称函数参数生成。
- `example/`
  示例输入、示例 notebook 和测试数据目录。

### `vasp_utils/`

面向具体 VASP 工作流的脚本集合，包含：

- `vasp_workflow_bulk/`
  例如 `eos`、`convergence`、`cohesive`、`surface_energy`、`neb`、`stretch`。
- `vasp_workflow_planar_defects/`
  例如 `gsfe`、`decohesion`、slip 面相关工作流。
- `vasp_workflow_others/`
  一些专用分析脚本。
- `vasp_universal/`
  通用批处理工具，如批量提交、重投、监控、清理、抽取结果。
- `n2p2_utils/`
  面向神经网络势相关的数据整理和主动学习脚本。

### `myvasp/`

偏个人化的 VASP 辅助脚本集合，包含：

- 清理、复制、grep、修改 INCAR 的 shell 脚本
- `vtstscripts/` 中的 VTST 相关工具副本
- 一些提交脚本、NEB/力分析/数据转换工具

### `slurm_utils/`

集群任务脚本模板目录，包含顺序提交、重投和一些通用命令模板。

### `lmp_utils/`

LAMMPS 输入模板，当前可见模板主要围绕：

- `stretch`
- `gsfe`
- `Cij_energy`

### `docs/`

Sphinx 文档工程和已生成的 HTML 文档。入口通常是：

- `docs/source/`
- `docs/build/html/index.html`

### `yin_github/`

保留的一份上游或历史版本脚本镜像，便于参考比对，不建议作为主入口。

## 4. `mymetal` 的主要能力

### 4.1 结构构建

`mymetal.build` 提供了较多结构生成函数，适合基于 ASE `Atoms` 做二次开发。

典型能力包括：

- 生成 FCC/HCP 体相与表面取向结构
- 从 bulk 切 slab / film
- 对薄膜或晶胞施加单轴或多轴拉伸
- 构造异质结并处理上下层堆叠
- 处理羟基化或薄膜几何调整

比较典型的入口函数：

- `mymetal.build.film.stretch.generate_film`
- `mymetal.build.film.stretch.stretch_list_along_direction_to_cell`
- `mymetal.build.bulk.create.create_fcc_111`
- `mymetal.build.bulk.create.create_hcp_basal`
- `mymetal.build.film.findhetero.find_hetero`

### 4.2 VASP I/O

`mymetal.io.vasp` 对 ASE 的 VASP 读写进行了定制：

- `my_read_vasp()`
  读取 `POSCAR/CONTCAR`，并返回 `Atoms` 与晶格比例因子
- `my_write_vasp()`
  写出 VASP 结构文件，并支持 `lattice_scale_factor`

这对需要保留 VASP 文件原始缩放因子的工作流比较有用。

### 4.3 后处理

`mymetal.post.newmain` 是一套批量后处理工具，围绕目录型计算任务读取 `OUTCAR` 并输出摘要文本。

可以提取的内容包括：

- 收敛状态
- 总能和熵项
- 压力、体积、最大受力
- 输入参数统计
- warning / advice 上下文

这套逻辑默认面向类似 `y_dir/0.997/OUTCAR` 这类批处理目录结构。

### 4.4 材料学计算辅助

`mymetal.calculate` 中已有一些可直接调用的小工具，例如：

- `calenergy.surfenergy.cal_surface_energy`
- `calmechanics.*`
- `calmismatch.calhetero`
- `calqm.kpoints`

这些函数更偏“研究脚本中的积木块”，而不是完整应用框架。

### 4.5 机器学习数据准备

`mymetal.ml.n2p2.dataset.nnpdata` 支持：

- 从现有 `input.data` 读入 n2p2 数据集
- 从 VASP `OUTCAR` 提取多帧结构、能量和力
- 写出 n2p2 所需数据格式

这部分适合把 DFT 计算结果接到神经网络势训练流程上。

## 5. 典型工作流

### 工作流 A：生成晶体或薄膜结构

常见流程是：

1. 用 `mymetal.build.bulk` 或 `mymetal.build.film` 生成初始结构
2. 用 `mymetal.io.vasp.my_write_vasp` 写出 `POSCAR`
3. 配合 `vasp_utils/` 或个人脚本提交任务

示例：

```python
from mymetal.build.film.stretch import generate_film

slab = generate_film(
    symbols="Au",
    structure="fcc",
    num_layers=12,
    my_vacuum=20.0,
    slice_plane=(1, 1, 1),
)
```

### 工作流 B：批量应变或拉伸计算

常见流程是：

1. 先生成参考薄膜或体相结构
2. 使用 `stretch_list_along_direction_to_cell()` 批量生成不同应变构型
3. 将每个构型写入不同子目录
4. 用 `vasp_utils` 或 `myvasp` 的脚本批量提交
5. 用 `mymetal.post` 统一汇总结果

### 工作流 C：表面能计算

常见流程是：

1. 构建 bulk 和 slab
2. 分别完成 VASP 计算
3. 读取能量、原子数和表面积
4. 调用 `cal_surface_energy()`

示例：

```python
from mymetal.calculate.calenergy.surfenergy import cal_surface_energy

gamma = cal_surface_energy(
    bulk_energy=-100.0,
    bulk_atoms_number=4,
    relaxed_surface_energy=-190.0,
    surface_atoms_number=8,
    area=25.0,
    energy_unit="eV",
)
```

### 工作流 D：VASP 输出转 n2p2 数据集

示例：

```python
from mymetal.ml.n2p2.dataset import nnpdata

data = nnpdata()
data.load_from_outcar(outcarfile="./OUTCAR", index=":", tag="train")
data.write(outfile_name="./input.data")
```

### 工作流 E：使用命令行脚本跑批量 VASP

`vasp_utils/README.md` 中已经展示了两类典型场景：

- 基于 `y_full_relax` 的派生工作流
- 基于 `y_src/poscars2/` 的批量结构提交

例如：

```shell
yin_vasp_run_eos
yin_vasp_plot_all -eos
yin_vasp_univ_sub_poscars2
yin_vasp_univ_post
```

这些脚本更依赖你的集群环境、PATH、提交脚本模板和目录约定。

## 6. 目录约定

仓库里多处代码默认使用一种研究工作目录组织方式：

```text
project_folder/
├─ y_full_relax/
├─ y_src/
│  ├─ INCAR
│  ├─ KPOINTS
│  ├─ POTCAR
│  ├─ sub.vasp
│  └─ poscars2/
└─ y_dir/
   ├─ 0.997/
   ├─ 0.998/
   └─ 1.000/
```

其中：

- `y_full_relax/` 常作为参考态目录
- `y_src/` 常用于保存不变的源文件
- `y_dir/` 常用于批量生成和运行的子任务目录

`mymetal.post.newmain` 的许多接口默认就是围绕这种目录结构写的。

## 7. 依赖

根据 `forclear/requirements.txt`，核心依赖大致包括：

- 数值与科学计算：`numpy`、`scipy`、`sympy`、`statsmodels`
- 材料建模：`ase`、`spglib`、`pymatgen`、`mp-api`
- 绘图与可视化：`matplotlib`、`ovito`、`brokenaxes`、`adjustText`
- 文档：`sphinx`、`sphinx_rtd_theme`、`recommonmark`
- 其他工具：`rich`、`typer`、`networkx`

可选依赖包括：

- 异质结构搜索：`hetbuilder`
- 机器学习：`torch`、`torchvision`、`scikit-learn`
- LAMMPS Python 接口：`lammps`、`mpi4py`

## 8. 安装建议

### 8.1 推荐方式：源码工作区方式使用

当前仓库根目录没有看到标准的 `pyproject.toml` 或可直接使用的根级 `setup.py`。现有打包脚本位于 `forclear/setup.py`，但它与当前源码布局并不完全对齐。

因此，最稳妥的使用方式是：

1. 克隆仓库
2. 安装依赖
3. 将仓库根目录加入 `PYTHONPATH`，或在 IDE / notebook 中以源码方式引用

示例：

```shell
pip install -r forclear/requirements.txt
```

然后在当前仓库根目录下运行 Python，或自行把根目录加入 `PYTHONPATH`。

### 8.2 如果你准备整理打包配置

仓库中存在以下打包痕迹：

- `forclear/setup.py`
- `mymetal_pkg.egg-info/`

这说明项目曾经以 `mymetal-pkg` 名称打包过，但当前仓库状态下，建议先补齐根目录打包配置，再使用：

```shell
pip install -e .
```

换句话说，README 这里刻意不把 `pip install -e .` 当成现成保证可用的命令。

## 9. 文档与示例

### Sphinx 文档

仓库自带 Sphinx 工程：

- 源文件：`docs/source/`
- 构建结果：`docs/build/html/`

如果只是快速浏览，可直接打开：

- `docs/build/html/index.html`

### 示例

`mymetal/example/` 下保留了大量研究型示例，包括：

- `test-generate-bulk`
- `test-surface-energy`
- `test-stretch`
- `test-post`
- `test-hydroxylated`
- `test-hetbuilder-fixatom`
- `test-n2p2-sfparams`

如果你想理解项目的真实使用方式，这些示例比 API 文档更重要。

## 10. 当前已知特点与限制

使用这个仓库前，建议先了解它的风格：

- 这是研究驱动的工具仓库，不是经过严格产品化收敛的公共库
- 代码中保留了较多历史脚本、镜像目录和个人工作流痕迹
- 某些模块依赖固定目录命名，例如 `y_dir/`
- 一些功能默认面向 VASP 文件组织方式，而不是纯函数式接口
- 打包配置目前不够统一，源码方式使用更稳妥
- 部分文档字符串较完整，但也存在历史注释、编码痕迹和中英文混合情况

这并不妨碍它在个人研究工作流中高效使用，但意味着第一次接手时更适合“从示例和已有目录模板反推用法”。

## 11. 建议的阅读顺序

如果你第一次看这个项目，建议按下面顺序：

1. 先看 `mymetal/__init__.py`
2. 再看 `mymetal/build/film/stretch.py`
3. 再看 `mymetal/io/vasp.py`
4. 再看 `mymetal/post/newmain.py`
5. 然后对照 `mymetal/example/` 中相近案例
6. 最后再看 `vasp_utils/README.md` 和对应工作流脚本

## 12. 一个简短总结

一句话概括，这个仓库是：

> 以 `mymetal` 为核心、围绕 VASP 材料计算任务组织起来的一套个人研究型工具链。

它的强项不在于“单一模块极致精简”，而在于把下面这些事情放在同一工作区里：

- 结构构建
- 批量任务准备
- 结果后处理
- 训练数据抽取
- 集群脚本模板

如果你准备继续维护这个仓库，最值得优先整理的通常是三件事：

1. 根目录打包配置
2. 面向 `mymetal` 的最小可运行示例
3. `vasp_utils` 与 `mymetal` 之间的职责边界
