# 跨阶段 Blocker 汇总

> 审阅 cron 写入，执行 cron 读取修复。
> 状态流转：open → fixed-by-executor → verified-by-reviewer

---

## [utils-stage2] vasp_workflow_bulk.rst 三个 Sphinx warning

- 严重度: major（不阻塞构建，但影响文档质量）
- 文件: `docs/source/manual/vasp_workflow_bulk.rst`
- 发现时间: 2026-08-03 基线构建
- 状态: open

### 问题 1: toctree 缺失
- 描述: `document isn't included in any toctree`
- 位置: 整个文件未注册到 index.rst 或 manual toctree
- 修复建议: 在 `docs/source/index.rst` 的 manual toctree 加入 `manual/vasp_workflow_bulk`

### 问题 2: unknown document 'neb_utils'
- 描述: `:doc:`neb_utils`` 指向不存在的文档（2 处）
- 位置: line 202, line 498（seealso 块）
- 修复建议: neb_utils 内容尚未独立成页，改为文字说明或移除 :doc: 引用，
  或在 Stage 2 补建 neb_utils.rst 页面

### 问题 3: inline strong 语法错误
- 描述: `**内部**` 紧跟全角内容无空格，RST 解析器报
  `Inline strong start-string without end-string`
- 位置: line 483
- 修复建议: 在 `**内部**` 后加空格，或改用 `` ``内部`` `` 反引号代码标记

### 验证方法
修复后重跑：
```bash
'C:/Users/louis/mysoft/env/pyenv/codex/Scripts/python.exe' -m sphinx -E -b html --keep-going docs/source docs/_build/html 2>&1 | grep vasp_workflow_bulk
```
应输出 0 行。

---

<!-- 后续 blocker 按此格式追加 -->
