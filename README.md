# pjvasp_package

`pjvasp_package` 是面向计算材料学的研究工具仓库。`mymetal` 是可导入的核心
Python package；`vasp_utils/`、`slurm_utils/`、`lmp_utils/`、`myvasp/` 与
`n2p2_utils/` 提供 VASP/HPC 等目录型 workflow 脚本。

适合的任务包括：

- 用 ASE/Python 构建 bulk、film、surface 与 heterostructure；
- 读写 POSCAR/CONTCAR，并汇总 OUTCAR 的能量、压力、体积、力与收敛状态；
- 组织 VASP、LAMMPS、SLURM 与 n2p2 工作流；
- 把可复用结构/后处理逻辑组合进个人科研脚本。

详细 Tutorial、Workflow 与 API 以 Sphinx 文档为唯一正式来源：

- 在线文档：<https://shenglin1001.github.io/pjvasp_package/>
- 本地入口：`docs/source/index.rst`
- 第一个完整案例：`docs/source/getting_started/au111_slab.rst`

## 最短安装

建议使用 Python 3.10 的隔离环境：

```shell
python -m pip install -r requirements.txt
python -m pip install -e .
python -c "import mymetal; print('mymetal: ok')"
```

部分模块还依赖 companion repository：

```shell
python -m pip install "git+https://github.com/ShengLin1001/myalloy_package.git@master"
```

安装、optional dependency、external executable 与 HPC module 的区别见
[Installation](https://shenglin1001.github.io/pjvasp_package/user_guide/install.html)。

## 一个无需 VASP 的结果

```python
from mymetal.build.film.stretch import generate_film

slab = generate_film(
    symbols="Au",
    structure="fcc",
    num_layers=12,
    my_vacuum=20.0,
    slice_plane=(1, 1, 1),
)
assert slab.get_chemical_formula() == "Au12"
```

完整脚本还会写出并重新读取 `POSCAR`、生成结构图并执行确定性断言：

```shell
python docs/examples/getting_started_au111.py --output docs/_build/example-au111
```

它不运行 VASP，不需要 POTCAR，也不会调用 `sbatch`。

## 仓库组成

| 路径 | 责任 |
| --- | --- |
| `mymetal/` | 结构构建、计算辅助、I/O、后处理、绘图与 ML 数据 |
| `vasp_utils/`、`myvasp/` | VASP workflow 与历史辅助脚本 |
| `slurm_utils/` | 软件无关的 SLURM 生成/提交入口与 runner |
| `lmp_utils/` | LAMMPS 模板和 workflow 脚本 |
| `n2p2_utils/` | n2p2 scaling/training 辅助脚本 |
| `docs/` | Sphinx canonical documentation、可运行示例与图片生成脚本 |
| `mymetal/example/` | 原始研究/回归素材，不等同于已发布 Tutorial |

## `pei_*` 脚本怎样调用

`setup.py` 当前没有 `console_scripts`，因此 `pip install -e .` 不会自动把
`pei_*` 命令安装到 PATH。可以显式运行仓库内文件：

```shell
python slurm_utils/slurm_universal/pei_slurm_univ_submit.py -h
bash vasp_utils/vasp_workflow_bulk/pei_vasp_run_neb -h
```

在 CentOS HPC 上也可以把经过审核的脚本目录加入 PATH。缺少 VASP、LAMMPS、
SLURM 或 n2p2 命令时，先执行 `module avail`，再按目标集群配置 module、
partition、launcher 和 executable。

所有提交教程默认 dry-run；不传 `-if_sbatch` 时
`pei_slurm_univ_submit.py` 只生成脚本。

## 本地构建文档

```shell
python -m pip install -r docs/requirements.txt
python -m pip install --no-deps "git+https://github.com/ShengLin1001/myalloy_package.git@master"
python -m sphinx -E -b html -W --keep-going docs/source docs/_build/html
```

然后打开 `docs/_build/html/index.html`。贡献与验证规则见
`docs/source/reference/development.rst`。

## 内容与安全边界

- README 是中文项目入口；Sphinx 维护完整 Tutorial、Workflow 与 Reference。
- 叙事内容使用中文，API、模块、命令、文件名与标准术语保留英文。
- 不发布 POTCAR，也不把未确认许可的 OUTCAR/结构打包为下载资产。
- 不把 VASP 计算完成、SLURM 提交成功和材料学结果已验证混为同一状态。
- `mymetal/example/` 中的 notebook 只有经过清理、复算和许可确认后才进入正式文档。
