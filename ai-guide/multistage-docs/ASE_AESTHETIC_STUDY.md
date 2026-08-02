# ASE 文档审美规范

来源：https://docs.ase-lib.org/（2026-08-02 实测提取）

ASE 文档使用 **PyData Sphinx Theme**（`pydata-sphinx-theme`），叠加 sphinx-book-theme
和 sphinx-gallery。整体风格简洁、学术、信息密度高，优先可读性。

## 主题与依赖

- HTML 主题：`pydata-sphinx-theme`（不是 sphinx_rtd_theme）
- 附加：sphinx-gallery（用于示例画廊）、pygments（语法高亮）
- 字体：system-ui sans-serif 正文 + ui-monospace 代码（Cascadia Code / Source Code Pro / Menlo / Consolas）

## 色板（实测 RGB）

| 元素 | 颜色 | 说明 |
|------|------|------|
| 页面背景 | `rgb(255,255,255)` 纯白 | 干净白底 |
| 正文文字 | `rgb(34,40,50)` 深灰偏蓝 | 不是纯黑，降低对比刺眼 |
| 一级标题 h1 | `rgb(41,49,61)` 深灰蓝 | 42px |
| 二级标题 h2 | `rgb(41,49,61)` | 34px |
| 正文链接 | `rgb(20,24,30)` 近黑 | 悬停可能变蓝，主链接低调 |
| 侧栏链接 | `rgb(34,40,50)` | 与正文同色 |
| 代码块背景 | `rgb(243,244,245)` 极浅灰 | 14px monospace |
| 代码块边框 | `1px solid rgb(209,213,218)` | 浅灰边框 |
| 代码文字 | `rgb(34,40,50)` / 行内 `rgb(72,86,107)` | 略偏蓝灰 |
| 行内 code 背景 | `rgb(243,244,245)` | 与块级一致 |
| 搜索按钮 | `rgb(232,146,23)` 橙色 | 唯一亮色点缀 |
| 表格边框 | `rgb(243,244,245)` | 极浅 |

## 排版参数

- 正文字号：16px
- 行高：26.4px（约 1.65）
- h1：42px；h2：34px
- 代码字号：14px
- 侧栏宽度：约 281px
- 内容区宽度：约 1126px

## 布局结构

1. **顶部 header**：白底，左侧 logo + 站名，右侧搜索框（橙色按钮）+ GitLab 图标链接。
   支持 light/dark 切换按钮。
2. **左侧主导航栏**：可折叠树形结构，顶层项如 About / Installation / Tutorials /
   Modules / Tips / Release notes。当前页高亮。
3. **主内容区**：纯白背景，宽度约 1126px。标题 → 段落 → 代码块 → 图片 → 表格
   交替。无装饰性 banner。
4. **右侧 intra-page TOC**（可选）：页内标题跳转。
5. **底部 footer**：简洁，含版权和生成工具。

## 视觉特征要点

1. **极简白底**：几乎全白，唯一亮色是搜索按钮的橙色。
2. **代码块浅灰**：背景 `#f3f4f5`，1px 浅灰边框，不抢眼。
3. **链接低调**：正文链接近黑色，靠下划线或悬停区分，不用亮蓝。
4. **图片居中**：figure 居中显示，带编号和 caption。
5. **表格紧凑**：浅灰边框，无斑马纹。
6. **admonition 轻量**：note/warning 用浅色背景 + 左侧色条，不夸张。
7. **无装饰**：没有大 banner、渐变、阴影、卡片阴影。信息密度优先。
8. **mathjax 公式**：内联和块级公式原生支持。

## 对 mymetal 文档的迁移策略

当前 mymetal 用 `sphinx_rtd_theme`。要达到 ASE 风格有两种路径：

### 路径 A（推荐，低风险）：保留 rtd_theme + 自定义 CSS 模拟 ASE 风格

- 不换主题（避免引入 pydata-sphinx-theme 的依赖和潜在兼容问题）
- 在 `_static/css/custom.css` 中：
  - 白底 `#ffffff`
  - 正文色 `rgb(34,40,50)`
  - 代码块背景 `#f3f4f5` + `1px solid #d1d5da` 边框
  - 链接色调暗
  - 加宽内容区到 1126px
  - 橙色点缀用于关键按钮/标签
- 优点：零新依赖，CI 不变，风险最低
- 缺点：无法完全复现 pydata 主题的侧栏折叠细节

### 路径 B（高保真）：切换到 pydata-sphinx-theme

- `docs/requirements.txt` 加 `pydata-sphinx-theme`
- `conf.py` 设 `html_theme = "pydata_sphinx_theme"`
- 优点：像素级复现 ASE 风格
- 缺点：需重新调 toctree 配置、可能影响现有 custom.css、CI 要测

**决策**：Stage 0 先走路径 A（低风险快速见效），如果时间允许且构建稳定，
可在 Stage 6 收尾时考虑路径 B 升级。当前以"内容丰富度"为第一优先级，
主题保真度为第二。

## 内容风格规范（学 ASE）

1. 每页开头一句话说清"这个对象/函数是什么"。
2. 紧跟一个最小可运行代码块（4-10 行）。
3. 关键参数用 definition list（`:term:` 风格）逐一说明。
4. get/set 方法用表格对照。
5. 复杂概念配图（结构图、流程图）。
6. 数学公式用 mathjax 块级渲染。
7. 交叉链接密集：`Atoms`、`write()` 等都链到对应模块页。
8. 末尾"See also"列出相关模块。
