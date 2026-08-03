#!/usr/bin/env python3
"""Overview figure for the vasp_utils/vasp_universal sub-package.

This tutorial is VASP-free. It demonstrates the core functionality of the
vasp_utils/vasp_universal sub-package by rendering a single multi-panel
figure that summarises:

  * the pei_vasp_univ_sbatch run/resume flow and its exit-code contract,
  * how pei_vasp_univ_post recursively scans ``y_dir`` trees,
  * the file-type coverage of the four clean-up scripts,
  * the INCAR-editing behaviour of pei_vasp_univ_find_and_change /
    pei_vasp_univ_increase_nbands, and
  * the 0 / 1 / 10 exit-code convention shared across the runner.

No VASP binary, ``sbatch``, or real calculation is invoked.  The figure is
generated with matplotlib only.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np


# ---------------------------------------------------------------------------
# Style constants (task spec: white bg, deep-grey text, light-grey code bg)
# ---------------------------------------------------------------------------
COLOR_TEXT = (34 / 255, 40 / 255, 50 / 255)
COLOR_CODE_BG = "#f3f4f5"
COLOR_PANEL_BG = "#fbfbfc"
COLOR_ARROW = "#4a5568"
COLOR_OK = "#2f855a"
COLOR_WARN = "#c05621"
COLOR_ERR = "#c53030"
COLOR_SKIP = "#6b7280"
COLOR_NODE = "#eef2f7"
COLOR_NODE_BORDER = "#b7c2d0"
COLOR_DECISION = "#fef3c7"
COLOR_DECISION_BORDER = "#d97706"
COLOR_RUN = "#dcfce7"
COLOR_EXIT10 = "#dbeafe"
COLOR_EXIT0 = "#dcfce7"
COLOR_EXIT1 = "#fee2e2"


def _set_rcparams() -> None:
    """Apply a consistent deep-grey-on-white matplotlib style.

    CJK glyphs render via Microsoft YaHei / Noto Sans SC / SimHei (whichever
    matplotlib finds first); Latin glyphs fall back to DejaVu Sans.  This keeps
    the Chinese panel labels and the English script/flag names legible.
    """
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            # CJK + Latin per-glyph fallback: Chinese glyphs resolve to
            # Microsoft YaHei, while ✓ / – / Latin fall back to DejaVu Sans.
            "font.family": ["Microsoft YaHei", "Noto Sans SC", "SimHei", "DejaVu Sans"],
            "axes.unicode_minus": False,
            "text.color": COLOR_TEXT,
            "axes.labelcolor": COLOR_TEXT,
            "axes.edgecolor": COLOR_TEXT,
            "xtick.color": COLOR_TEXT,
            "ytick.color": COLOR_TEXT,
            "axes.titlesize": 13,
            "axes.titleweight": "bold",
            "axes.titlepad": 8,
            "font.size": 9,
        }
    )


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------
def _box(ax, x, y, w, h, text, fc=COLOR_NODE, ec=COLOR_NODE_BORDER, fontsize=8.5,
         fontweight="normal", text_color=None, style="round,pad=0.02"):
    """Draw a labelled box centred at (x, y) and return its bounding artist."""
    box = mpatches.FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle=style, fc=fc, ec=ec, lw=1.1, zorder=3,
    )
    ax.add_patch(box)
    ax.text(
        x, y, text, ha="center", va="center",
        fontsize=fontsize, fontweight=fontweight,
        color=text_color or COLOR_TEXT, zorder=4,
    )
    return box


def _code(ax, x, y, w, h, text, fontsize=8):
    """Draw a light-grey 'code block' box with monospace-ish text.

    Uses a CJK-aware font family list so that any Chinese characters embedded
    in the code samples render correctly alongside the Latin/monospace tokens.
    """
    box = mpatches.FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle="round,pad=0.02", fc=COLOR_CODE_BG, ec="#d1d5db",
        lw=0.8, zorder=3,
    )
    ax.add_patch(box)
    ax.text(
        x, y, text, ha="center", va="center",
        fontsize=fontsize, color=COLOR_TEXT, zorder=4,
        family=["DejaVu Sans Mono", "Microsoft YaHei", "Noto Sans SC", "SimHei", "DejaVu Sans"],
    )
    return box


def _arrow(ax, x1, y1, x2, y2, label=None, color=COLOR_ARROW, lw=1.4,
           label_offset=(0.12, 0.0), label_color=None):
    """Draw an arrow between two points, optionally labelled near the middle."""
    ax.annotate(
        "",
        xy=(x2, y2), xytext=(x1, y1),
        arrowprops=dict(arrowstyle="-|>", color=color, lw=lw,
                        shrinkA=2, shrinkB=2),
        zorder=2,
    )
    if label:
        mx = (x1 + x2) / 2 + label_offset[0]
        my = (y1 + y2) / 2 + label_offset[1]
        ax.text(
            mx, my, label, ha="center", va="center", fontsize=7.5,
            color=label_color or color, fontweight="bold", zorder=5,
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85),
        )


def _turn_off(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


# ---------------------------------------------------------------------------
# Panel 1: pei_vasp_univ_sbatch flow
# ---------------------------------------------------------------------------
def draw_panel_sbatch_flow(ax):
    """Flowchart of pei_vasp_univ_sbatch: cd → check OUTCAR → run/exit."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 11)
    ax.set_title("面板1 · pei_vasp_univ_sbatch 单目录运行/续算流程", loc="left")

    cx = 5.0
    nodes = [
        (cx, 10.2, "cd $job_dir  (realpath)", COLOR_NODE, "start"),
        (cx, 9.0, "读 INCAR 判断 NSW/IBRION\n自动选择 completion_marker\n(relax / static)", COLOR_DECISION, "marker"),
        (cx, 7.6, "is_complete() ?\nOUTCAR 含 marker +\nenergy + timing", COLOR_DECISION, "check1"),
        (1.9, 7.6, "exit 10\n已完成, 跳过", COLOR_EXIT10, "exit10"),
        (cx, 6.1, "CONTCAR 非空 ?\n是 → cp CONTCAR POSCAR\n(检测 EDIFF 过大则修 INCAR)\n否 → 用现有 POSCAR", COLOR_DECISION, "check2"),
        (cx, 4.4, "清理旧输出:\nCHG/CONTCAR/DOSCAR/OUTCAR/\nvasprun.xml ...", COLOR_NODE, "clean"),
        (cx, 3.0, "srun vasp_std  (MY_LAUNCHER)", COLOR_RUN, "run"),
        (cx, 1.6, "is_complete() ?\n再检查 OUTCAR", COLOR_DECISION, "check3"),
        (1.9, 1.6, "exit 0\n本轮完成", COLOR_EXIT0, "exit0"),
        (cx, 0.4, "exit 1\n仍未收敛 / 缺文件", COLOR_EXIT1, "exit1"),
    ]
    for x, y, text, fc, _name in nodes:
        if x == cx:
            _box(ax, x, y, 4.4, 0.95, text, fc=fc, fontsize=8)
        else:
            _box(ax, x, y, 2.0, 0.95, text, fc=fc, fontsize=8.5, fontweight="bold")

    # arrows
    _arrow(ax, cx, 9.72, cx, 9.5)                 # cd -> marker
    _arrow(ax, cx, 8.52, cx, 8.08)                # marker -> check1
    _arrow(ax, 3.8, 7.6, 2.9, 7.6, label="是 Yes", label_offset=(0, 0.18))
    _arrow(ax, cx, 7.12, cx, 6.58, label="否 No", label_offset=(0, 0.18))
    _arrow(ax, cx, 5.62, cx, 4.88)                # check2 -> clean
    _arrow(ax, cx, 3.92, cx, 3.48)                # clean -> run
    _arrow(ax, cx, 2.52, cx, 2.08)                # run -> check3
    _arrow(ax, 3.8, 1.6, 2.9, 1.6, label="是 Yes", label_offset=(0, 0.18),
           label_color=COLOR_OK)
    _arrow(ax, cx, 1.12, cx, 0.88, label="否 No", label_offset=(0, 0.18),
           label_color=COLOR_ERR)

    # legend strip
    legend_items = [
        ("exit 0 = 本轮收敛完成", COLOR_EXIT0),
        ("exit 10 = 早已完成, 跳过", COLOR_EXIT10),
        ("exit 1 = 失败 / 未收敛 / 缺文件", COLOR_EXIT1),
        ("判定框", COLOR_DECISION),
    ]
    for i, (txt, col) in enumerate(legend_items):
        lx = 0.4 + i * 2.45
        ax.add_patch(mpatches.Rectangle((lx, 10.85), 0.28, 0.18, fc=col,
                                        ec=COLOR_NODE_BORDER, lw=0.6, zorder=3))
        ax.text(lx + 0.36, 10.94, txt, ha="left", va="center", fontsize=7,
                color=COLOR_TEXT, zorder=4)
    _turn_off(ax)


# ---------------------------------------------------------------------------
# Panel 2: pei_vasp_univ_post y_dir tree scan
# ---------------------------------------------------------------------------
def draw_panel_post_tree(ax):
    """Directory-tree diagram showing how pei_vasp_univ_post scans y_dir."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 11)
    ax.set_title("面板2 · pei_vasp_univ_post 递归扫描 y_dir 目录树", loc="left")

    # top: root
    _box(ax, 5, 10.3, 4.6, 0.7,
         "运行目录 (. / y_hoec_energy)\nfind ./ -type d -name y_dir",
         fc=COLOR_DECISION, fontsize=8)

    # two y_dir groups
    y_dirs = [(1.8, 8.6, "y_dir/  (group A)"), (8.2, 8.6, "y_dir/  (group B)")]
    for x, y, t in y_dirs:
        _box(ax, x, y, 2.8, 0.6, t, fc=COLOR_NODE, fontsize=8.5, fontweight="bold")
        _arrow(ax, 5, 9.95, x, 8.9, lw=1.2)

    # sub-calc dirs under each y_dir
    sub_a = [
        (0.7, 7.3, "1_0.00\nOUTCAR slurm-*.out"),
        (1.8, 7.3, "2_-0.01\nOUTCAR.gz"),
        (2.9, 7.3, "3_+0.01\n(无 OUTCAR)"),
    ]
    sub_b = [
        (7.1, 7.3, "4_0.00\nOUTCAR"),
        (8.2, 7.3, "5_-0.02\nOUTCAR"),
        (9.3, 7.3, "6_+0.02\nOUTCAR"),
    ]
    for x, y, t in sub_a + sub_b:
        _box(ax, x, y, 1.0, 0.95, t, fc="#f0f7ff", ec="#9db8d8", fontsize=6.6)

    for x in [0.7, 1.8, 2.9]:
        _arrow(ax, 1.8, 8.3, x, 7.78, lw=0.9, color="#7a8aa0")
    for x in [7.1, 8.2, 9.3]:
        _arrow(ax, 8.2, 8.3, x, 7.78, lw=0.9, color="#7a8aa0")

    # aggregation block
    _box(ax, 5, 5.7, 8.6, 1.7,
         "对每个 y_dir 的父目录调用 sub_post():\n"
         "· 解压 OUTCAR.gz → grep energy/EENTRO/stress/volume\n"
         "· 统计 total / un-finished / un-relaxed 作业数\n"
         "· diff INCAR/KPOINTS/POTCAR/sub.vasp 与参考目录\n"
         "· 收集 WARNING / highest band 提示",
         fc=COLOR_CODE_BG, fontsize=7.6)
    _arrow(ax, 1.8, 6.82, 3.8, 6.55, lw=1.0, color="#7a8aa0")
    _arrow(ax, 8.2, 6.82, 6.2, 6.55, lw=1.0, color="#7a8aa0")

    # output files
    outs = [
        (1.2, 3.6, "y_post_data.txt\nenergy / EENTRO / stress"),
        (3.6, 3.6, "y_post_data_2.txt\nvolume / pressure / Fmax"),
        (6.0, 3.6, "y_post_time.txt\niter / relaxed / time / CPU"),
        (8.6, 3.6, "y_post_param.txt\nENCUT/ISMEAR/EDIFF/..."),
    ]
    for x, y, t in outs:
        _box(ax, x, y, 2.2, 1.0, t, fc="#eef6ee", ec="#8fb88f", fontsize=7)
        _arrow(ax, 5, 4.85, x, 4.1, lw=0.9, color="#7a8aa0")

    # verbose extras
    _box(ax, 5, 2.1, 9.0, 1.2,
         "-v 模式额外运行:\n"
         "pei_vasp_univ_extract_convergence  →  y_post_convergence/*.txt\n"
         "pei_vasp_plot_convergence.py        →  收敛曲线 PNG\n"
         "yin_vasp_plot_mag/dos/statistics.py (YIN_GITHUB 旧脚本)",
         fc=COLOR_CODE_BG, fontsize=7.2)
    _arrow(ax, 5, 3.1, 5, 2.7, lw=0.9, color="#7a8aa0")

    # summary line
    _box(ax, 5, 0.7, 9.4, 0.9,
         "收尾: 逐 y_dir 打印 y_post_time 汇总 + 'Your highest band' 警告计数",
         fc=COLOR_DECISION, fontsize=7.8)
    _arrow(ax, 5, 1.5, 5, 1.15, lw=0.9, color="#7a8aa0")

    _turn_off(ax)


# ---------------------------------------------------------------------------
# Panel 3: clean-up scripts comparison table
# ---------------------------------------------------------------------------
def draw_panel_clean_table(ax):
    """Table comparing the four clean-up scripts' file coverage."""
    ax.set_title("面板3 · 清理脚本文件类型对比表", loc="left")
    ax.axis("off")

    headers = ["文件类型", "clean_up_full", "clean_up_small", "clean_old_slurm", "clean_outcar"]
    # rows: file-type group, then per-script behaviour (✓ / —)
    rows = [
        ("VASP 主输出 (OUTCAR/OSZICAR/CONTCAR/DOSCAR...)", "✓", "✓", "—", "—"),
        ("WAVECAR / CHGCAR / PARCHG / EIGENVAL",          "✓", "✓", "—", "—"),
        ("vasprun.xml / REPORT / SUMMARY.*",              "✓", "✓", "—", "—"),
        ("wannier90.* / *.dat / plotfile / p4vasp.log",   "✓", "✓", "—", "—"),
        ("Slurm 日志 (*.e* / *.o* / *.err / *.out)",      "✓", "✓", "仅保留最新", "—"),
        ("y_* 后处理文件 (y_post_*.txt 等)",              "✓", "—", "—", "—"),
        ("evsv/function/ne/nnp/sc/sf/test/train/*.norm",  "✓", "—", "—", "—"),
        ("slurm-*.out 旧副本 (按 mtime)",                 "—", "—", "✓", "—"),
        ("OUTCAR 末尾 \\x00 null 字符",                    "—", "—", "—", "✓ (rstrip)"),
    ]

    col_widths = [0.40, 0.155, 0.155, 0.155, 0.135]
    n_cols = len(headers)
    n_rows = len(rows) + 1
    table = ax.table(
        cellText=rows,
        colLabels=headers,
        colWidths=col_widths,
        cellLoc="center",
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(7.6)
    table.scale(1.0, 1.55)

    # style header
    for j in range(n_cols):
        cell = table[0, j]
        cell.set_facecolor("#2c3e50")
        cell.set_text_props(color="white", fontweight="bold", fontsize=8)
        cell.set_edgecolor("white")

    # style body: first column left-aligned, color ✓/—
    color_yes = COLOR_OK
    color_no = COLOR_SKIP
    color_special = "#b45309"
    for i in range(1, n_rows):
        cell0 = table[i, 0]
        cell0.set_text_props(ha="left", color=COLOR_TEXT, fontsize=7.4)
        cell0.PAD = 0.04
        for j in range(1, n_cols):
            cell = table[i, j]
            txt = rows[i - 1][j]
            if txt == "✓":
                cell.set_text_props(color=color_yes, fontweight="bold")
            elif txt == "—":
                cell.set_text_props(color=color_no)
            else:
                cell.set_text_props(color=color_special, fontweight="bold", fontsize=6.8)
            # zebra
            if i % 2 == 0:
                cell.set_facecolor("#f7f9fb")
            cell0 = table[i, 0]
            if i % 2 == 0:
                cell0.set_facecolor("#f7f9fb")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


# ---------------------------------------------------------------------------
# Panel 4: INCAR operation scripts
# ---------------------------------------------------------------------------
def draw_panel_incar_ops(ax):
    """Explain pei_vasp_univ_find_and_change and pei_vasp_univ_increase_nbands."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 11)
    ax.set_title("面板4 · INCAR 操作脚本", loc="left")

    # ---- find_and_change ----
    _box(ax, 2.6, 10.0, 4.4, 0.7, "pei_vasp_univ_find_and_change",
         fc="#2c3e50", text_color="white", fontweight="bold", fontsize=8.5)

    _box(ax, 2.6, 8.9, 4.4, 1.3,
         "通用 INCAR/POSCAR/KPOINTS 标签编辑器\n"
         "递归 find ./ -type f -name <file>\n"
         "改/加/注释首处匹配行 (VASP 只用首次出现)",
         fc=COLOR_CODE_BG, fontsize=7.2)
    _arrow(ax, 2.6, 9.65, 2.6, 9.55)

    # dedicated cases
    cases = [
        ("INCAR 专用 (-tag value)", [
            "-nbands 18", "-encut 520", "-ediff 1e-6",
            "-isif 3", "-nsw 200", "-algo Normal",
            "-ismear 0", "-lreal Auto", "-ncore 4",
        ]),
        ("注释 / 改写", [
            "-ediff comment      (注释原行)",
            "-algo comment:Normal (改写并注释)",
        ]),
        ("POSCAR / KPOINTS", [
            "-a0 1.0  (POSCAR 第2行 scale)",
            "-cn Mg  (POSCAR 第6行 元素)",
            "-kp 12 12 12  (KPOINTS 第4行)",
        ]),
        ("任意未列出的 tag", [
            "-my_tag value  → 大写写入 INCAR",
            "(查 incar_tag_help_vaspwiki.tsv 补注释)",
        ]),
    ]
    y0 = 7.7
    for title, items in cases:
        _box(ax, 2.6, y0, 4.4, 0.4, title, fc="#e8edf3", ec="#b7c2d0",
             fontsize=7.2, fontweight="bold")
        block = "\n".join(items)
        _code(ax, 2.6, y0 - 0.95, 4.2, 1.5, block, fontsize=6.8)
        y0 -= 2.4

    # ---- increase_nbands ----
    _box(ax, 7.4, 10.0, 4.4, 0.7, "pei_vasp_univ_increase_nbands",
         fc="#2c3e50", text_color="white", fontweight="bold", fontsize=8.5)

    _box(ax, 7.4, 8.9, 4.4, 1.3,
         "从 OUTCAR 读 NELECT / NIONS / NBANDS_old\n"
         "按 NBANDS = ceil(NELECT/2 + 2·NIONS) 计算\n"
         "再调用 find_and_change -nbands <新值>",
         fc=COLOR_CODE_BG, fontsize=7.2)
    _arrow(ax, 7.4, 9.65, 7.4, 9.55)

    # formula visual
    _code(ax, 7.4, 7.4, 4.2, 1.1,
          "awk:  NBANDS = NELECT/2 + 2*NIONS\n"
          "      取 int 后向上取整 (value > int → +1)",
          fontsize=7.0)

    # output echo
    _code(ax, 7.4, 5.9, 4.2, 1.2,
          "$ pei_vasp_univ_increase_nbands OUTCAR\n"
          "NELECT=48\nNIONS=16\n"
          "NBANDS=56      NBANDS_OLD=40",
          fontsize=6.8)

    # typical trigger
    _box(ax, 7.4, 4.3, 4.4, 1.5,
         "触发场景:\n"
         "OUTCAR 出现 'Your highest band is\n"
         "occupied at some k-points' 警告时,\n"
         "自动增加未占据带以加速收敛 / 防止带交叉",
         fc="#fff7ed", ec="#d97706", fontsize=7.2)
    _arrow(ax, 7.4, 5.3, 7.4, 5.05)

    # relation arrow between the two
    _arrow(ax, 5.2, 10.0, 5.2, 8.9, lw=0)
    ax.annotate(
        "increase_nbands 内部\n调用 →", xy=(5.6, 9.4), xytext=(5.6, 9.4),
        ha="center", va="center", fontsize=7, color=COLOR_ARROW,
        fontweight="bold",
    )

    # bottom note
    _box(ax, 5, 0.7, 9.6, 0.9,
         "两个脚本都先校验所有 option/value 对再编辑任何文件, 避免部分写入;"
         " edit 原子性 (mktemp → mv), 保留原文件权限",
         fc=COLOR_PANEL_BG, ec="#d1d5db", fontsize=7.4)
    _turn_off(ax)


# ---------------------------------------------------------------------------
# Panel 5: exit-code convention
# ---------------------------------------------------------------------------
def draw_panel_exit_codes(ax):
    """Table of exit 0 / 1 / 10 meanings across the runner + helpers."""
    ax.set_title("面板5 · 退出码约定 (pei_vasp_univ_sbatch 及配套脚本)", loc="left")
    ax.axis("off")

    headers = ["退出码", "含义", "典型触发", "调度器/上层处理建议"]
    rows = [
        ("0", "本轮计算收敛完成",
         "OUTCAR 含 completion_marker +\nenergy(without entropy) +\nGeneral timing",
         "该目录标记为完成, 不再重投;\n可进入 post / 收敛曲线绘制"),
        ("1", "失败 / 未收敛 / 缺文件",
         "srun 返回非 0; 或 CONTCAR/POSCAR\n均缺失; 或 OUTCAR 仍无收敛标记",
         "需人工排查 (查 slurm-*.out / OUTCAR);\n可改 INCAR 后 resubmit"),
        ("10", "早已完成, 跳过运行",
         "is_complete() 在运行前即为真\n(OUTCAR 已有完整收敛标记)",
         "跳过该目录, 不消耗机时;\npost 时归入 un-finished=0 统计"),
    ]

    col_widths = [0.09, 0.20, 0.38, 0.33]
    n_cols = len(headers)
    n_rows = len(rows) + 1
    table = ax.table(
        cellText=rows,
        colLabels=headers,
        colWidths=col_widths,
        cellLoc="left",
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(7.6)
    table.scale(1.0, 2.3)

    # header
    for j in range(n_cols):
        cell = table[0, j]
        cell.set_facecolor("#2c3e50")
        cell.set_text_props(color="white", fontweight="bold", fontsize=8.5)
        cell.set_edgecolor("white")

    # exit-code colour band
    exit_colors = {1: COLOR_EXIT0, 2: COLOR_EXIT1, 3: COLOR_EXIT10}
    for i in range(1, n_rows):
        code_cell = table[i, 0]
        code_cell.set_facecolor(exit_colors[i])
        code_cell.set_text_props(ha="center", fontweight="bold", fontsize=11,
                                 color=COLOR_TEXT)
        for j in range(1, n_cols):
            cell = table[i, j]
            cell.set_text_props(color=COLOR_TEXT, fontsize=7.4)
            if i % 2 == 0:
                cell.set_facecolor("#f7f9fb")
        if i % 2 == 0:
            code_cell.set_facecolor(exit_colors[i])

    # footer note
    ax.text(0.5, 0.02,
            "注: pei_vasp_univ_check_phase_transition 另有 exit 2/3 (跳过/检查失败); "
            "此处仅列 runner 主流程的 0/1/10。",
            ha="center", va="bottom", fontsize=7, color=COLOR_SKIP,
            transform=ax.transAxes)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


# ---------------------------------------------------------------------------
# Figure assembly
# ---------------------------------------------------------------------------
def render_overview(path_image: Path) -> None:
    """Build the 5-panel overview figure and save it to path_image."""
    _set_rcparams()
    fig = plt.figure(figsize=(16, 12), dpi=150)
    # 3 rows x 2 cols, last cell (panel 5) spans bottom-right
    gs = fig.add_gridspec(
        3, 2, hspace=0.32, wspace=0.12,
        left=0.03, right=0.98, top=0.95, bottom=0.03,
    )
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    ax5 = fig.add_subplot(gs[2, :])

    draw_panel_sbatch_flow(ax1)
    draw_panel_post_tree(ax2)
    draw_panel_clean_table(ax3)
    draw_panel_incar_ops(ax4)
    draw_panel_exit_codes(ax5)

    # suptitle
    fig.suptitle(
        "vasp_utils / vasp_universal 子包功能总览",
        fontsize=16, fontweight="bold", color=COLOR_TEXT, y=0.985,
    )

    fig.savefig(path_image, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # ---- non-blank self-check (mandatory per pjvasp-examples skill) ----
    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError("overview image was not created: " + str(path_image))
    image_rgb = plt.imread(path_image)[..., :3]
    if float(np.mean(np.abs(image_rgb - 1.0))) < 0.002:
        raise AssertionError("overview image is effectively blank: " + str(path_image))


def run_example(path_output: Path) -> dict:
    """Resolve output dir, render the PNG, run deterministic assertions."""
    path_output = Path(path_output).expanduser().resolve()
    path_output.mkdir(parents=True, exist_ok=True)
    path_image = path_output / "vasp_universal_overview.png"

    render_overview(path_image)

    size_kb = path_image.stat().st_size / 1024.0
    image_rgb = plt.imread(path_image)[..., :3]
    deviation = float(np.mean(np.abs(image_rgb - 1.0)))

    summary = {
        "image_path": str(path_image),
        "size_kb": round(size_kb, 1),
        "pixel_deviation": round(deviation, 4),
    }
    print("wrote: " + str(path_image))
    print("size_kb: " + str(round(size_kb, 1)))
    print("pixel_deviation: " + str(round(deviation, 4)))
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render the vasp_utils/vasp_universal overview figure (VASP-free).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("docs/_build/example-vasp-universal"),
        help="Output directory for the rendered PNG (default: docs/_build/example-vasp-universal)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary = run_example(args.output)
    # deterministic assertions
    assert summary["size_kb"] > 10.0, "PNG smaller than 10 KB: " + str(summary)
    assert summary["pixel_deviation"] > 0.002, "PNG effectively blank: " + str(summary)
    print("OK: vasp_universal overview figure verified non-blank.")


if __name__ == "__main__":
    main()
