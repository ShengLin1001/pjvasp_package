# 多阶段文档网站升级 — 总计划

## 目标

对 `mymetal` package 的各个子包持续生成功能描述、可运行示例、结构图/结果图，
并部署到 GitHub Pages 在线文档。审美参考 ASE 文档（https://docs.ase-lib.org/）。

## 工作分支

`docs/mymetal-website-upgrade`（已从 main 创建）。所有修改在此分支进行。
main 分支保持不动。最终验证通过后再 merge 回 main 或提 PR。

## 关键约束（无人值守 yolo 模式安全规则）

1. 禁止使用 `rm -rf`、`rmdir /s`、`rd /s`、`del /s`、`Remove-Item -Recurse`
   或任何递归删除命令。需要删除文件时，一次只删一个明确指定的文件路径。
2. 不得修改 `setup.py`（用户原有修改必须保留）。
3. 不得运行真实 VASP、LAMMPS、n2p2 training 或 `sbatch`。
4. 不得发布 POTCAR、用户路径、账号或集群 secret。
5. 示例脚本不得调用真实 VASP/POTCAR；用合成数据或仓库内已 tracked fixture。
6. 不大规模重建整个 `docs/`；不覆盖用户已有改动。
7. 每个阶段完成后必须运行 Sphinx 构建 + compileall 确认无回归。
8. 中文叙事，英文 API/模块/命令名称。
9. commit message 遵循仓库风格：`docs(...)`、`feat(...)` + 简短中文描述。

## 阶段划分

| 阶段 | 内容 | 估时 |
|------|------|------|
| Stage 0 | 主题/CSS 升级到 PyData 风格 + sphinx 配置 + 基线构建测试 | 1 轮 |
| Stage 1 | build 子包（bulk + film + workflow） | 1-2 轮 |
| Stage 2 | calculate 子包（7 个模块） | 1-2 轮 |
| Stage 3 | io + post 子包 | 1 轮 |
| Stage 4 | universal 子包（atom/check/data/math/matrix/plot） | 1 轮 |
| Stage 5 | ml + slurm + cr 子包 | 1 轮 |
| Stage 6 | 审查收尾：交叉链接、首页、移动端、linkcheck、全量构建、commit、push、部署 | 1 轮 |

## 每阶段通用产出

每个子包阶段必须产出：

1. **可运行示例脚本**（`docs/examples/<topic>.py`）：无需 VASP，含确定性断言，
   生成结构图或结果图。
2. **RST 文档页**（`docs/source/tutorials/<topic>.rst` 或 `docs/source/api/<sub>.rst`）：
   包含目标、前置条件、`literalinclude` 引用同一脚本、预期输出、结果图、
   验证方法、常见错误、相关 API 交叉链接。
3. **生成的图片**（`docs/source/_static/images/generated/<topic>.png`）：
   非空白，有 caption 和 alt text。
4. **API reference 补充**：在 `docs/source/api/<sub>.rst` 用 `autofunction`/
   `automodule` 暴露真实函数签名、参数、返回值、单位、异常、最小示例。
5. **图片生成脚本注册**：在 `docs/scripts/generate_structure_images.py` 注册
   新示例的 render 函数。
6. **CI workflow 更新**：在 `.github/workflows/docs.yml` 的 examples 步骤注册
   新示例脚本。
7. **导航更新**：在相关 `index.rst` 或 `api.rst` 的 toctree 注册新页面。

## 验证清单（每阶段结束必须执行）

```bash
# 在仓库根目录
python -m compileall mymetal
python -m sphinx -E -b html -W --keep-going docs/source docs/_build/html
python docs/scripts/generate_structure_images.py
# 检查图片非空白
ls -la docs/source/_static/images/generated/
```

可选（linkcheck 较慢，每 2-3 阶段跑一次）：
```bash
python -m sphinx -b linkcheck -W --keep-going docs/source docs/_build/linkcheck
```

## 子agent分工模式

每阶段可并行调用最多 3 个 subagent：

- **Agent A（示例编写者）**：读源码 → 写 `docs/examples/<topic>.py` → 运行验证 →
  产出图片路径。
- **Agent B（文档撰写者）**：根据 Agent A 的脚本和图片，写 RST 页面、API reference、
  更新 toctree 和 generate_structure_images.py。
- **Agent C（审查者）**：根据 ASE 审美规范和本计划检查 Agent B 的产出，反馈意见。

串行依赖：Agent A 先行 → Agent B 跟进 → Agent C 审查。但不同子模块的 A 可并行。

## 状态追踪

每阶段开始/结束时更新 `ai-guide/multistage-docs/STAGE_STATUS.md`。
每阶段结束时写 `ai-guide/multistage-docs/handoff/stage_N_handoff.md`。
