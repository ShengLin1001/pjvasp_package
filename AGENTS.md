# 仓库指南

## 项目结构与模块组织

本仓库是一个以 `mymetal` Python 包为核心的计算材料学工具包。核心代码位于 `mymetal/` 目录下，其中包含用于结构构建的子包 `build/`、用于计算流程的 `calculate/`、用于 VASP 输入/输出的 `io/`、用于后处理的 `post/`、用于机器学习工具的 `ml/`，以及通用辅助函数 `universal/`。工作流脚本和模板单独存放：`vasp_utils/` 和 `myvasp/` 用于 VASP，`slurm_utils/` 用于调度器脚本，`lmp_utils/` 用于 LAMMPS 模板，`matlab_code/` 用于 MATLAB 分析。Sphinx 文档源文件位于 `docs/source/`；生成的 HTML 文档位于 `docs/build/`。

## 构建、测试与开发命令

- `python -m pip install -e .`：以可编辑模式安装该包。
- `python -m pip install -r requirements.txt`：安装较完整的运行时和文档依赖。
- `python -m pip install -r setup/strict-requirements.txt`：安装固定版本依赖，用于构建可复现环境。
- `cd docs && make html`：构建 Sphinx 文档，并输出到 `docs/build/html/`。
- `python -m compileall mymetal`：对主包进行快速语法检查。

在当前 CentOS HPC 平台上，如果所需命令缺失或版本过旧，在手动安装前应先检查 `module avail`。

## 代码风格与命名规范

尽可能使用兼容 Python 3.10 的代码。遵循 PEP 8：使用四空格缩进，函数和变量采用 `snake_case`，类名采用 `CamelCase`，模块名应具有描述性。工作流脚本应尽量靠近其对应领域，例如 VASP 批处理辅助脚本应放在 `vasp_utils/vasp_universal/` 中。当逻辑需要共享时，优先在 `mymetal/` 中实现可复用函数，而不是编写大型一次性脚本。

## 测试指南

目前尚无专门的测试套件。对于 Python 代码修改，至少应运行 `python -m compileall mymetal`；在实际可行时，应在受影响的工作流附近添加小型可执行示例。未来测试文件应命名为 `test_*.py`，并放置在 `tests/` 中，或放在需要 VASP/LAMMPS 固定输入文件的示例旁边。

## Commit 与 Pull Request 指南

近期提交使用简洁前缀，例如 `docs:` 和 `feat(DFT):`，通常附带简短中文描述。请保持这种风格：说明变更类型、可选作用域以及具体结果。Pull request 应描述受影响的工作流，列出验证命令，说明所需的 HPC module 或外部工具，并在涉及后处理或文档变更时提供修改前后的输出片段。

## Agent 专用指令

不要批量删除文件或目录。严禁使用递归删除命令，例如 `rm -rf`、`rmdir /s`、`rd /s`、`del /s` 或 `Remove-Item -Recurse`。如果确实需要删除，只能一次删除一个明确指定的文件路径；如果需要批量清理，应停止操作，并要求用户手动删除相关文件。