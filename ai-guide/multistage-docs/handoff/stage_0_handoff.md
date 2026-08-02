# Stage 0 Handoff — 主题/CSS 升级 + 基线构建

## 完成时间
2026-08-02

## 修改的文件

1. `docs/source/_static/css/custom.css` — 完全重写，实现 ASE 风格色板：
   - 白底 `#ffffff`
   - 正文 `rgb(34,40,50)`，标题 `rgb(41,49,61)`
   - 代码块背景 `rgb(243,244,245)` + `1px solid rgb(209,213,218)` 边框
   - 链接 `rgb(11,92,173)` 蓝
   - admonition: note 蓝、warning 红、important 橙
   - sidebar 浅灰底 `rgb(248,249,250)`
   - task-card 首页网格
   - 移动端 600px 响应式

2. `docs/source/conf.py` — 新增扩展：
   - `sphinx.ext.viewcode`（[source] 链接到源码）
   - `sphinx.ext.intersphinx`（交叉链接 ASE/Python/NumPy）
   - `sphinx_copybutton`（代码块复制按钮，学 ASE）
   - `intersphinx_mapping` 指向 ASE、Python、NumPy

3. `docs/requirements.txt` — 新增 `sphinx-copybutton>=0.5`

## 验证结果

- `sphinx -E -b html -W --keep-going`：**成功**（build succeeded）
- `sphinx_copybutton` 已安装到 codex venv
- ASE 色板在 custom.css 中应用

## 未完成项

- intersphinx 对 ASE 的 `objects.inv` 拉取：ASE 文档站点可能没有标准
  `objects.inv`，intersphinx 会在 linkcheck 时报告但不影响 html 构建。
  后续阶段如需 `:class:\`ase.Atoms\`` 交叉链接，需确认 ASE 的
  intersphinx inventory 可用性。当前不影响功能。

## 下一阶段建议

Stage 1（build 子包）可以开始。重点：
- bulk/create.py 的 create_fcc_111、create_hcp_basal 还没有 tutorial
- film/hydroxyl.py 的 passivate_surface_custom 没有文档
- workflow/hoec.py 和 workflow/kpar_ncore.py 的目录生成器需要示例
- 已有的 generate_film、cal_area、bulk_structures tutorial 保持不动

## 构建命令备忘

```bash
# codex venv 有所有依赖
'C:/Users/louis/mysoft/env/pyenv/codex/Scripts/python.exe' -m sphinx -E -b html -W --keep-going docs/source docs/_build/html
'C:/Users/louis/mysoft/env/pyenv/codex/Scripts/python.exe' -m compileall mymetal
'C:/Users/louis/mysoft/env/pyenv/codex/Scripts/python.exe' docs/scripts/generate_structure_images.py
```
