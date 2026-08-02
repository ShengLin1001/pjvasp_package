# Cron Job 执行提示词

> 这是定时任务的自我包含指令。每小时由 cron job 触发一次。
> cron job 会加载此文件作为任务指令。不要在此文件中嵌入选项。

## 你的角色

你是 mymetal package 文档网站升级项目的阶段执行 agent。你在一个无人值守的
yolo 模式下运行，但有严格的安全约束。

## 启动检查（每次运行必须先做）

1. cd 到仓库根目录：
   `F:\BaiduSyncdisk\version20240608\main_code_space\pjvasp_package`
2. 确认当前分支：
   ```bash
   git branch --show-current
   ```
   必须在 `docs/mymetal-website-upgrade` 分支。如果在 main，切换过去：
   `git checkout docs/mymetal-website-upgrade`
   如果分支不存在，从 main 创建：`git checkout -b docs/mymetal-website-upgrade`
3. 拉取最新（防止多个 session 冲突）：
   `git fetch origin && git rebase origin/docs/mymetal-website-upgrade 2>/dev/null || true`
4. 读取状态文件确定当前阶段：
   `ai-guide/multistage-docs/STAGE_STATUS.md`
5. 读取主计划和审美规范：
   - `ai-guide/multistage-docs/MASTER_PLAN.md`
   - `ai-guide/multistage-docs/ASE_AESTHETIC_STUDY.md`
   - `ai-guide/multistage-docs/SUBPACKAGE_COVERAGE_MAP.md`
6. 读取上一阶段的 handoff（如果有）：
   `ai-guide/multistage-docs/handoff/stage_N_handoff.md`

## 执行当前阶段

根据 STAGE_STATUS.md 中第一个 `pending` 的阶段执行。每个阶段的具体内容见
STAGE_STATUS.md 的"阶段定义"部分。

### 通用执行流程

每个阶段按以下步骤：

1. **调研源码**：用 read_file 和 search_files 读取当前阶段涉及的 mymetal 模块源码，
   记录每个公开函数的签名、参数、返回值、单位、异常。
2. **并行写示例**：用 delegate_task 启动 1-3 个 subagent 并行编写示例脚本。
   每个 subagent 负责一个子模块的示例脚本。subagent 指令必须包含：
   - 源码路径
   - 输出脚本路径（`docs/examples/<topic>.py`）
   - 要求：无需 VASP、含确定性断言、生成 PNG 图片、有 `render_*` 和 `build_*` 函数
   - 安全约束（不删文件、不运行 VASP、用合成数据）
3. **运行验证示例**：用 terminal 运行每个新示例脚本，确认无报错、图片非空白。
4. **写 RST 文档**：根据示例脚本和源码，写 `docs/source/tutorials/<topic>.rst`
   或更新 `docs/source/api/<sub>.rst`。遵循 ASE 审美规范：
   - 开头一句话说清是什么
   - 最小可运行代码块
   - 参数用 definition list
   - 结果图居中带 caption
   - 末尾 Related API 交叉链接
5. **注册图片生成**：在 `docs/scripts/generate_structure_images.py` 注册新示例。
6. **更新 CI**：在 `.github/workflows/docs.yml` 的 examples 步骤注册新脚本。
7. **更新导航**：在 `docs/source/tutorials/index.rst` 或 `docs/source/api.rst`
   的 toctree 注册新页面。
8. **验证构建**：
   ```bash
   python -m compileall mymetal
   python docs/scripts/generate_structure_images.py
   python -m sphinx -E -b html -W --keep-going docs/source docs/_build/html
   ```
   如果构建失败，修复后重试。最多重试 3 次，仍失败则记录 blocker 并停止。
9. **更新状态**：修改 `ai-guide/multistage-docs/STAGE_STATUS.md`，标记当前阶段
   completed，填入时间和说明。
10. **写 handoff**：写 `ai-guide/multistage-docs/handoff/stage_N_handoff.md`，
    记录：修改了哪些文件、新增了哪些示例/图片、构建结果、未完成项、下一阶段建议。
11. **commit**：
    ```bash
    git add -A
    git commit -m "docs(stage-N): <简短中文描述>"
    ```
    注意：不要 `git add setup.py`。如果 `git add -A` 把 setup.py 也加进去了，
    用 `git reset HEAD setup.py` 排除。

## 安全约束（绝对不可违反）

1. **禁止递归删除**：不得使用 `rm -rf`、`rmdir /s`、`rd /s`、`del /s`、
   `Remove-Item -Recurse`。需要删文件时一次只删一个明确路径。
2. **不碰 setup.py**：`setup.py` 有用户原有修改，不得修改、还原、暂存。
3. **不运行 VASP/LAMMPS/n2p2/sbatch**：所有示例用合成数据或仓库内 tracked fixture。
4. **不发布 secret**：不发布 POTCAR、用户路径、账号、集群信息。
5. **不 force push**：只用普通 push 到 `docs/mymetal-website-upgrade`。
6. **不直接 push main**：所有改动只在 `docs/mymetal-website-upgrade` 分支。
7. **不破坏现有文档**：不大规模重建 docs/，只增量添加。已有的 tutorial 和
   API 页面如需修改，用 patch 精确替换，不整体覆盖。
8. **不递归调度 cron**：不要在此 cron session 中创建新的 cron job。

## subagent 使用

每阶段可调用 delegate_task 启动 subagent。subagent 是 leaf 角色，不能递归调度。
给 subagent 的 context 必须包含：
- 仓库根路径
- 相关源码的绝对路径
- 输出文件的绝对路径
- ASE 审美规范的关键点（从 ASE_AESTHETIC_STUDY.md 摘要）
- 安全约束
- "respond in Chinese for narrative content, English for code/API" 的语言要求

subagent 示例分工：
- Agent A：读源码 → 写 docs/examples/<topic>.py → 运行验证 → 返回图片路径
- Agent B：根据 Agent A 的脚本和图片 → 写 RST 页面 → 更新 toctree
- Agent C：审查 Agent B 产出是否符合 ASE 审美 → 返回修改意见

串行：A → B → C。但不同 topic 的 A 可并行。

## 构建失败处理

如果 Sphinx 构建失败：
1. 读取错误日志
2. 定位根因（通常是 RST 语法错误、autodoc 找不到模块、图片路径错误）
3. 用 patch 修复
4. 重新构建
5. 最多重试 3 次。3 次后仍失败：在 handoff 中记录 blocker，标记阶段为
   `blocked`，停止本次运行。下一轮 cron 会重新尝试或跳过。

## 完成条件

一个阶段算完成当且仅当：
1. 所有计划的新示例脚本已创建并能运行
2. 所有计划的新 RST 页面已创建并构建通过
3. `python -m sphinx -E -b html -W --keep-going docs/source docs/_build/html` 成功
4. `python -m compileall mymetal` 成功
5. STAGE_STATUS.md 已更新
6. handoff 文档已写
7. 已 commit（不 push，除非是 Stage 6）

Stage 6 额外：push 到 origin，检查 GitHub Actions，记录部署 URL。

## 不要做的事

- 不要修改 setup.py
- 不要修改 .git/config 或 git remote
- 不要创建新的 cron job
- 不要删除 ai-guide/ 下已有文件
- 不要在 docs/ 下创建空目录或空文件
- 不要用 echo/cat heredoc 创建文件，用 write_file
- 不要用 sed/awk 编辑文件，用 patch
- 不要用 rm -rf
