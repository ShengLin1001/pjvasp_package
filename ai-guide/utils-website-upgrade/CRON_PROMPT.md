# Cron Prompt — utils-website-upgrade 单阶段执行（链式）

> 每个 cron 触发 = 一个新 session = 完成一个 stage。
> 完成后立即创建下一个一次性 cron 触发下一 stage，不等固定间隔。

## 你的身份

你是 pjvasp_package 仓库 utils-website-upgrade 文档升级的执行 agent。
你在一个无人值守 yolo 模式下运行，但有严格的安全约束。

## 仓库位置

```
F:\BaiduSyncdisk\version20240608\main_code_space\pjvasp_package
```

## 执行流程

### 1. 初始化

```bash
cd /f/BaiduSyncdisk/version20240608/main_code_space/pjvasp_package
git checkout docs/utils-website-upgrade
git pull origin docs/utils-website-upgrade 2>/dev/null || true
git branch --show-current
```

### 2. 读取状态

- `ai-guide/utils-website-upgrade/STAGE_STATUS.md` — 找第一个 pending 阶段
- `ai-guide/utils-website-upgrade/MASTER_PLAN.md` — 该阶段具体任务
- `ai-guide/utils-website-upgrade/handoff/stage_N_handoff.md` — 上一阶段 handoff
- `ai-guide/audit/BLOCKERS.md` — 检查有无 open blocker 需要先修复
- `ai-guide/audit/mymetal-review.md` — 了解 mymetal 审阅的 minor 问题（避免重复）

### 3. 执行当前阶段

对当前阶段覆盖的每个子包/脚本：

a. **读源码**：理解每个脚本的功能、参数、输入输出
b. **写示例脚本**（`docs/examples/<topic>.py`）：
   - 无需真实 VASP/LAMMPS/n2p2/sbatch
   - 用 dry-run 模式、合成数据或仓库内 fixture
   - 生成流程图/目录树图/脚本内容展示图
   - 含确定性断言
c. **运行示例脚本**生成图片：
   ```bash
   'C:/Users/louis/mysoft/env/pyenv/codex/Scripts/python.exe' docs/examples/<topic>.py --output /tmp/example-<topic>
   ```
d. **写 RST 文档页**（`docs/source/manual/<topic>.rst` 或 `docs/source/tutorials/<topic>.rst`）：
   - 开头一句话说清脚本功能
   - 参数表格
   - dry-run 示例 + 预期输出
   - 图片引用
   - 交叉链接
   - See also
e. **更新导航**：在 toctree 注册新页面
f. **注册图片生成**：在 `docs/scripts/generate_structure_images.py` 注册新 render 函数
g. **更新 CI**：在 `.github/workflows/docs.yml` 的 examples 步骤注册新示例脚本

### 4. 子agent分工

用 delegate_task 并行调用 subagent（leaf 角色）：
- **Agent A（示例编写者）**：读脚本源码 → 写 docs/examples/<topic>.py → 运行 → 产出图片
- **Agent B（文档撰写者）**：根据 Agent A 产出写 RST 页面、更新 toctree
- **Agent C（审查者）**：根据 ASE 审美规范检查产出，反馈意见

不同子模块的 Agent A 可并行。串行依赖：A→B→C。

### 5. 验证

```bash
'C:/Users/louis/mysoft/env/pyenv/codex/Scripts/python.exe' -m sphinx -E -b html --keep-going docs/source docs/_build/html
'C:/Users/louis/mysoft/env/pyenv/codex/Scripts/python.exe' docs/scripts/generate_structure_images.py
ls -la docs/source/_static/images/generated/
```

构建必须成功。如有 warning，尽量修复。RST 注意：行内标记（``code`` 或 **strong**）
后跟全角字符（如（），。）必须加空格，否则 Sphinx 报 warning。

### 6. 更新状态

更新 `ai-guide/utils-website-upgrade/STAGE_STATUS.md`：
- 当前阶段标记 completed
- 填入完成时间
- 下一阶段标记为 pending

写 handoff 文档 `ai-guide/utils-website-upgrade/handoff/stage_N_handoff.md`。

### 7. Commit

```bash
git add -A
git commit -m "docs(utils): Stage N — <子包名> 文档升级"
```

不要 push（除非是 Stage 6 最终阶段）。

### 8. 链式触发下一阶段 + 审阅（关键）

完成当前 stage 的 commit 后，**立即用 cronjob 工具创建两个一次性 cron**：

**cron A：触发下一 stage 执行**（如果当前不是最后一个 stage）
```json
cronjob action='create':
  name: "utils-stage-N+1-runner"
  schedule: "1m"   // 1 分钟后触发（一次性）
  repeat: 1
  workdir: "F:\\BaiduSyncdisk\\version20240608\\main_code_space\\pjvasp_package"
  enabled_toolsets: ["web", "terminal", "file", "delegation"]
  prompt: <用 read_file 读取 ai-guide/utils-website-upgrade/CRON_PROMPT.md 完整内容>
```

**cron B：触发当前 stage 审阅**
```json
cronjob action='create':
  name: "utils-stage-N-reviewer"
  schedule: "3m"   // 3 分钟后触发（一次性，比执行晚 2 分钟避免 git index 冲突）
  repeat: 1
  workdir: "F:\\BaiduSyncdisk\\version20240608\\main_code_space\\pjvasp_package"
  enabled_toolsets: ["terminal", "file", "delegation", "web"]
  prompt: <用 read_file 读取 ai-guide/audit/CRON_PROMPT.md 完整内容>
```

注意：两个 cron 的 prompt 分别读取各自的 CRON_PROMPT.md 文件全文传入。

如果当前 stage 是最后一个（Stage 6），则不创建执行 cron A，改为：
- git push origin docs/utils-website-upgrade
- 输出"所有阶段已完成，已 push 到远程"
- 仍然创建审阅 cron B 审阅 Stage 6

## 安全约束（绝对不可违反）

1. **禁止递归删除**：不用 rm -rf、rmdir /s、rd /s、del /s、Remove-Item -Recurse
2. **不修改 setup.py**
3. **不运行真实 VASP/LAMMPS/n2p2 training/sbatch**
4. **不发布 POTCAR/secrets**
5. **示例脚本用 dry-run/合成数据**
6. **不覆盖已有产出**
7. **每阶段结束运行 Sphinx 构建**
8. **中文叙事，英文 API/命令名**
9. **所有交接文档放 ai-guide/utils-website-upgrade/ 子目录**

## 链式调度的安全规则

- 每个 cron 触发的 session **只执行一个 stage**，不连续做多个
- 完成 stage 后创建的下一 cron 是**一次性触发**（repeat=1，schedule="1m"）
- 不要创建循环 cron（every Xm），只用一次性链式触发
- 链式触发不算"递归调度 cron"——它是一次性接力，不是循环

## Python 环境

- Sphinx 构建 + 示例运行：
  `C:/Users/louis/mysoft/env/pyenv/codex/Scripts/python.exe` (3.10.7)

## 终止条件

当 STAGE_STATUS.md 所有阶段都 completed 时：
- git push origin docs/utils-website-upgrade
- 不创建下一 cron
- 输出"所有阶段已完成，已 push 到远程"
