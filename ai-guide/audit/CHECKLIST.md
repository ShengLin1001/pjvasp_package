# 每阶段审查清单

> 每个审阅 session 对照本清单逐项检查，在 REVIEW_REPORT 里记录结果。

## A. 示例脚本（docs/examples/<topic>.py）

- [ ] A1. 脚本存在且非空
- [ ] A2. 用 codexpy 环境运行无报错：
      `C:/Users/louis/mysoft/env/pyenv/codex/Scripts/python.exe docs/examples/<topic>.py --output /tmp/test`
- [ ] A3. 生成 PNG 图片，文件 > 5KB
- [ ] A4. 图片非空白（PIL 打开，pixel variation > 阈值）
- [ ] A5. 不调用真实 VASP/sbatch/POTCAR/LAMMPS/n2p2 training
- [ ] A6. 含确定性断言（assert 或自检函数）
- [ ] A7. 有 render_* 和 build_* 函数（如 MASTER_PLAN 要求）

## B. RST 文档页（docs/source/manual 或 tutorials/<topic>.rst）

- [ ] B1. 文件存在且非空
- [ ] B2. 开头一句话说清功能
- [ ] B3. 最小可运行示例紧跟（literalinclude 或 code-block）
- [ ] B4. 参数用 definition list 或表格说明
- [ ] B5. 结果图引用（figure/image 带 caption + alt）
- [ ] B6. 交叉链接（:doc:/:ref:）指向真实存在的目标
- [ ] B7. 末尾 See also 列相关模块
- [ ] B8. 中文叙事，英文 API/命令/模块名
- [ ] B9. 无 RST 语法 warning（inline strong/table/section）

## C. 覆盖度

- [ ] C1. RST 覆盖 MASTER_PLAN 该阶段列出的所有脚本/模块
- [ ] C2. handoff 自述的新增文件真实存在（逐个 ls 核对）
- [ ] C3. handoff 自述的更新文件确实被更新（git diff 核对）

## D. 集成

- [ ] D1. 新页面注册到 toctree（index.rst 或 manual 的 toctree）
- [ ] D2. 图片生成脚本注册到 generate_structure_images.py
- [ ] D3. CI workflow 注册新示例脚本
- [ ] D4. 无 "document isn't included in any toctree" warning

## E. 构建与安全

- [ ] E1. Sphinx 构建成功（--keep-going）
- [ ] E2. 记录 warning 数量（mymetal 轮次要求 -W 即 0 warning）
- [ ] E3. setup.py 未被修改（git diff setup.py 为空）
- [ ] E4. 无递归删除命令（grep -r "rm -rf\|rmdir /s\|del /s" docs/examples/）
- [ ] E5. 无 POTCAR/secret/用户路径泄露

## F. 审美（ASE 风格）

- [ ] F1. 图片风格统一（白底、深灰文字、浅灰代码块背景）
- [ ] F2. 表格紧凑无斑马纹
- [ ] F3. admonition 轻量
- [ ] F4. 信息密度高，无装饰性 banner
