# 跨阶段 Blocker 汇总

> 审阅 cron 写入，执行 cron 读取修复。
> 状态流转：open → fixed-by-executor → verified-by-reviewer

---

## [utils-stage2] vasp_workflow_bulk.rst 三个 Sphinx warning

- 严重度: major（不阻塞构建，但影响文档质量）
- 文件: `docs/source/manual/vasp_workflow_bulk.rst`
- 发现时间: 2026-08-03 基线构建
- 状态: fixed-by-executor（审阅 session 修复，2026-08-03）

### 问题 1: toctree 缺失
- 描述: `document isn't included in any toctree`
- 位置: 整个文件未注册到 index.rst 或 manual toctree
- 修复: 执行 cron 在暂停前已注册到 index.rst toctree

### 问题 2: unknown document 'neb_utils'
- 描述: `:doc:`neb_utils`` 指向不存在的文档（2 处）
- 位置: line 202, line 498（seealso 块）
- 修复: 执行 cron 创建了 neb_utils.rst，引用现已有效

### 问题 3: inline strong 语法错误
- 描述: `**内部**` 紧跟全角内容无空格
- 位置: line 483
- 修复: 审阅 session 在 `**内部**` 后加空格

### 额外修复
审阅 session 还修复了 neb_utils.rst 和 vasp_workflow_bulk.rst 中另外 7 处
同类 inline literal/strong warning（全角标点紧跟行内标记）。
修复后 Sphinx 构建 0 warning。

### 验证方法
```bash
'C:/Users/louis/mysoft/env/pyenv/codex/Scripts/python.exe' -m sphinx -E -b html --keep-going docs/source docs/_build/html 2>&1 | grep -E "WARNING|ERROR"
```
输出 0 行，已验证。

---

<!-- 后续 blocker 按此格式追加 -->
