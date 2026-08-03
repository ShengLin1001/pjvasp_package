# 文档网站升级审阅框架

> 审阅任务的多 session 协调中枢。每个审阅 session 开始时先读本文件。
> 审阅与执行（cron）分离：执行 cron 推进 Stage，审阅 cron 在阶段完成后介入。

## 1. 任务背景

pjvasp_package 有两轮文档网站升级任务：

| 轮次 | 分支 | 状态 | 执行 cron |
|------|------|------|-----------|
| mymetal | `docs/mymetal-website-upgrade` | Stage 0-6 全部完成，已 push origin | 无（已完成） |
| utils | `docs/utils-website-upgrade` | Stage 0-1 完成，2-6 pending，cron 推进中 | `b18aaee0e93a`（每 30m） |

mymetal 轮次由执行 agent 独立完成全部 7 阶段，**此前无独立审查**。
utils 轮次执行 cron 在跑，审阅需在每个 Stage 完成后跟进。

## 2. 审阅目标

参考 myobsidian HANDOFF.md 的"一人干活一人审查"模式：审查者必须用真实文件/命令
核实执行 agent 的自述，不能只转述 handoff。

每轮审阅产出一份 `REVIEW_REPORT_stageN.md`，记录：
- 审阅范围（哪些文件/脚本/RST/图片）
- 审查方法（实际运行的命令）
- 发现的问题（blocker / major / minor 分级）
- 通过项
- 建议（给执行 agent 或下一审阅 session）

## 3. 审查标准（ASE 审美 + 计划要求）

### 3.1 示例脚本可运行性
- 用 `C:/Users/louis/mysoft/env/pyenv/codex/Scripts/python.exe` 跑每个示例脚本
- 无报错，生成非空白 PNG
- 不调用真实 VASP/sbatch/POTCAR/LAMMPS/n2p2 training
- 含确定性断言

### 3.2 图片非空白
- PNG 文件 > 5KB
- 用 PIL 打开检查 pixel variation（非纯白/纯色）

### 3.3 RST 符合 ASE 审美
- 开头一句话说清功能
- 最小可运行示例紧跟
- 参数用 definition list 或表格
- 结果图居中带 caption + alt text
- 交叉链接密集（:doc:/:ref:）
- 末尾 See also 列相关模块
- 白底风格（CSS 检查）

### 3.4 交叉链接完整性
- Sphinx 构建 0 个 `unknown document` / `unknown reference` warning
- `:doc:` 引用都有对应目标文件

### 3.5 覆盖度
- RST 是否覆盖 MASTER_PLAN 里列的所有脚本/模块
- handoff 自述的产出文件是否真实存在

### 3.6 安全约束
- 不修改 setup.py
- 无递归删除命令（rm -rf 等）
- 不发布 POTCAR/secrets/用户路径
- 中文叙事，英文 API/命令/模块名

### 3.7 构建状态
- `sphinx -E -b html --keep-going` 成功
- warning 数量（utils 轮次不强制 -W，但记录 warning）
- mymetal 轮次要求 `-W`（强制 0 warning）

## 4. 已知问题（本 session 基线构建时发现）

utils Stage 2 未提交产物 `docs/source/manual/vasp_workflow_bulk.rst` 有 3 个 warning：

1. **toctree 缺失**：`document isn't included in any toctree`（未注册到 index.rst）
2. **unknown document 'neb_utils'**：line 202, 498 的 `:doc:`neb_utils`` 指向不存在的文档（2 处）
3. **inline strong 语法错误**：line 483 `**内部**` 紧跟全角内容无空格分隔

这些是 Stage 2 中间状态，执行 cron 后续会处理。审阅 session 应在 Stage 2 标记
completed 后复查这 3 项是否已修复。

## 5. 审阅 cron 协调规则

- 审阅 cron 与执行 cron `b18aaee0e93a` 错开：执行 cron 每 30m，审阅 cron 每 60m
- 审阅 cron 先读 `ai-guide/utils-website-upgrade/STAGE_STATUS.md`，只审阅
  **最近一个 completed 阶段**（不审阅 pending 或 in-progress）
- 审阅是只读操作：不修改 docs/ 下任何文件，只写 `ai-guide/audit/` 下的报告
- 如发现 blocker，在报告里标记，不自行修复（避免与执行 cron 冲突）；
  blocker 汇总到 `ai-guide/audit/BLOCKERS.md` 供执行 cron 下一轮读取修复

## 6. 文件组织

```
ai-guide/audit/
├── REVIEW_FRAMEWORK.md          # 本文件，审阅框架
├── CHECKLIST.md                 # 每阶段审查清单（可勾选）
├── BLOCKERS.md                  # 跨阶段 blocker 汇总（执行 cron 读）
├── mymetal-review.md            # mymetal Stage 0-6 整体审阅报告
├── utils-stage0-1-review.md     # utils Stage 0-1 审阅报告
└── utils-stageN-review.md       # utils 各阶段审阅报告（按需创建）
```

## 7. 与执行 cron 的通信协议

执行 cron 和审阅 cron 通过文件通信（不直接对话）：

- 审阅 → 执行：`ai-guide/audit/BLOCKERS.md`（审阅写入，执行读取修复）
- 执行 → 审阅：`ai-guide/utils-website-upgrade/handoff/stage_N_handoff.md`
  + `STAGE_STATUS.md`（执行写入，审阅读取对照）

BLOCKERS.md 格式：
```
## [stage-N] <问题简述>
- 严重度: blocker/major/minor
- 文件: <路径>
- 描述: <具体问题>
- 状态: open / fixed-by-executor / verified-by-reviewer
```

## 8. p-git-commit 提交节点

按用户要求，关键节点用 p-git-commit skill 提交并 push：
- 审阅报告写完后：提交 `ai-guide/audit/` 下报告
- blocker 修复验证后：提交验证记录
- commit message 格式：`docs(audit): :memo: <阶段>审阅报告` 或
  `docs(audit): :white_check_mark: <阶段>blocker验证`

提交在 utils 分支进行（审阅报告与 docs 同分支），push 到 origin。
