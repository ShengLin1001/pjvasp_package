#!/usr/bin/env python3
"""Overview figure for the slurm_utils sub-package.

This tutorial is sbatch-free. It demonstrates the core functionality of the
slurm_utils sub-package (slurm_universal / slurm_vasp / slurm_n2p2) by rendering
a single multi-panel figure that summarises:

  * the three ``-mode`` (parallel / each_subdir / single_alloc) of
    pei_slurm_univ_submit.py,
  * how pei_slurm_univ_submit.py recursively discovers y_dir and the file
    layout it generates,
  * the retry-strategy comparison of pei_slurm_univ_launch_retry vs
    pei_slurm_univ_sbatch_retry,
  * the preset registry of pei_slurm_univ_submit.py, and
  * the engine dispatch of pei_slurm_univ_monitor_error and the option list of
    pei_slurm_univ_useful_command.sh.

No ``sbatch``, no real job, no real calculation is invoked.  The figure is
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
COLOR_VASP = "#ede9fe"
COLOR_VASP_BORDER = "#7c3aed"
COLOR_N2P2 = "#fef3c7"
COLOR_N2P2_BORDER = "#d97706"
COLOR_HEADER_BG = "#2c3e50"


def _set_rcparams() -> None:
    """Apply a consistent deep-grey-on-white matplotlib style.

    CJK glyphs render via Microsoft YaHei / SimHei (whichever matplotlib finds
    first); Latin glyphs fall back to DejaVu Sans.  This keeps the Chinese panel
    labels and the English script/flag names legible.
    """
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "font.family": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
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
    """Draw a light-grey 'code block' box with monospace-ish text."""
    box = mpatches.FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle="round,pad=0.02", fc=COLOR_CODE_BG, ec="#d1d5db",
        lw=0.8, zorder=3,
    )
    ax.add_patch(box)
    ax.text(
        x, y, text, ha="center", va="center",
        fontsize=fontsize, color=COLOR_TEXT, zorder=4,
        family=["DejaVu Sans Mono", "Microsoft YaHei", "SimHei", "DejaVu Sans"],
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
# Panel 1: three -mode comparison
# ---------------------------------------------------------------------------
def draw_panel_modes(ax):
    """Flowchart comparing parallel / each_subdir / single_alloc."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 11)
    ax.set_title("面板1 · pei_slurm_univ_submit.py 三种 -mode 对比", loc="left")

    # top: shared entry
    _box(ax, 5, 10.3, 8.0, 0.7,
         "递归发现 K 个 y_dir/<case> 计算目录 (统一排序后按 -chunks 分流)",
         fc=COLOR_DECISION, fontsize=8)
    _arrow(ax, 5, 9.95, 5, 9.7, lw=1.0)
    _box(ax, 5, 9.3, 4.4, 0.7,
         "按 -mode 选择提交策略",
         fc=COLOR_NODE, fontsize=8.5, fontweight="bold")

    # three columns
    col_x = [1.85, 5.0, 8.15]
    col_titles = [
        ("parallel", COLOR_RUN, "每目录一个独立作业, 不等待"),
        ("each_subdir", COLOR_EXIT10, "1 个 shared parent (-n1) 管 K worker"),
        ("single_alloc", COLOR_EXIT0, "K 个计算资源父作业, 各自顺序执行"),
    ]
    for x, (name, fc, desc) in zip(col_x, col_titles):
        _arrow(ax, 5, 8.95, x, 8.7, lw=1.0)
        _box(ax, x, 8.35, 2.9, 0.85, name + "\n" + desc,
             fc=fc, fontsize=7.6, fontweight="bold")

    # parallel column
    _box(ax, 1.85, 7.2, 2.9, 0.85,
         "for each <case>:\n  sbatch sub_slurm_univ.sh",
         fc=COLOR_CODE_BG, fontsize=7.2)
    _arrow(ax, 1.85, 7.93, 1.85, 7.63, lw=0.9)
    _box(ax, 1.85, 5.95, 2.9, 0.95,
         "K 个独立作业并行排队\n互不等待, 各自 -child_wall_time",
         fc=COLOR_NODE, fontsize=7.0)
    _arrow(ax, 1.85, 6.78, 1.85, 6.43, lw=0.9)
    _box(ax, 1.85, 4.7, 2.9, 1.0,
         "-chunks: 忽略\n-chunk_parent_layout: 忽略\n(本就是每目录一作业)",
         fc="#fff7ed", ec="#d97706", fontsize=6.9)
    _arrow(ax, 1.85, 5.48, 1.85, 5.2, lw=0.9)
    _box(ax, 1.85, 3.5, 2.9, 0.85,
         "无 parent.sh\n仅生成 sub_slurm_univ.sh",
         fc=COLOR_PANEL_BG, fontsize=7.0)
    _arrow(ax, 1.85, 4.2, 1.85, 3.93, lw=0.9)

    # each_subdir column
    _box(ax, 5.0, 7.2, 2.9, 0.85,
         "sbatch parent.sh (-n1)\nparent 内并发 K 个 bash worker",
         fc=COLOR_CODE_BG, fontsize=7.0)
    _arrow(ax, 5.0, 7.93, 5.0, 7.63, lw=0.9)
    _box(ax, 5.0, 5.95, 2.9, 0.95,
         "worker 逐个:\n  sbatch --wait sub_slurm_univ.sh",
         fc=COLOR_NODE, fontsize=7.0)
    _arrow(ax, 5.0, 6.78, 5.0, 6.43, lw=0.9)
    _box(ax, 5.0, 4.7, 2.9, 1.0,
         "默认 -chunk_parent_layout auto\n→ shared (1 父作业)\n可改 per_chunk 恢复 K 父作业",
         fc="#fff7ed", ec="#d97706", fontsize=6.9)
    _arrow(ax, 5.0, 5.48, 5.0, 5.2, lw=0.9)
    _box(ax, 5.0, 3.5, 2.9, 0.85,
         "生成 slurm/ 下:\n  chunk001..K.sh + parent.sh",
         fc=COLOR_PANEL_BG, fontsize=7.0)
    _arrow(ax, 5.0, 4.2, 5.0, 3.93, lw=0.9)

    # single_alloc column
    _box(ax, 8.15, 7.2, 2.9, 0.85,
         "sbatch K 个 parent.sh\n每个 chunk 独占计算资源",
         fc=COLOR_CODE_BG, fontsize=7.0)
    _arrow(ax, 8.15, 7.93, 8.15, 7.63, lw=0.9)
    _box(ax, 8.15, 5.95, 2.9, 0.95,
         "parent 内顺序执行 cmd:\n  cd <case> && $cmd",
         fc=COLOR_NODE, fontsize=7.0)
    _arrow(ax, 8.15, 6.78, 8.15, 6.43, lw=0.9)
    _box(ax, 8.15, 4.7, 2.9, 1.0,
         "-chunk_parent_layout auto\n→ per_chunk (K 父作业)\nshared 在本模式被禁止",
         fc="#fff7ed", ec="#d97706", fontsize=6.9)
    _arrow(ax, 8.15, 5.48, 8.15, 5.2, lw=0.9)
    _box(ax, 8.15, 3.5, 2.9, 0.85,
         "生成 slurm/ 下:\n  chunk001..K.sh (无 child)",
         fc=COLOR_PANEL_BG, fontsize=7.0)
    _arrow(ax, 8.15, 4.2, 8.15, 3.93, lw=0.9)

    # bottom legend / dry-run note
    _box(ax, 5, 2.0, 9.6, 1.1,
         "通用约定:\n"
         "· -if_sbatch 默认 False (dry-run, 只生成脚本); 裸写 -if_sbatch 等价 True\n"
         "· -child_wall_time 仅 parallel / each_subdir 的子作业;  -parent_wall_time 仅父作业",
         fc=COLOR_DECISION, fontsize=7.3)
    _arrow(ax, 1.85, 3.08, 3.5, 2.55, lw=0.7, color="#7a8aa0")
    _arrow(ax, 5.0, 3.08, 5.0, 2.55, lw=0.7, color="#7a8aa0")
    _arrow(ax, 8.15, 3.08, 6.5, 2.55, lw=0.7, color="#7a8aa0")

    _box(ax, 5, 0.7, 9.6, 0.75,
         "三模式区别核心: parallel=零编排; each_subdir=轻编排(-n1 守 K worker); "
         "single_alloc=重编排(K 份完整资源, 各自顺序跑)",
         fc=COLOR_PANEL_BG, ec="#d1d5db", fontsize=7.2)

    _turn_off(ax)


# ---------------------------------------------------------------------------
# Panel 2: directory discovery + file layout
# ---------------------------------------------------------------------------
def draw_panel_discovery_layout(ax):
    """Directory tree + generated file layout."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 11)
    ax.set_title("面板2 · 目录发现与生成文件布局", loc="left")

    # discovery (left half)
    _box(ax, 2.5, 10.3, 4.4, 0.7,
         "pei_slurm_univ_submit.py\n-path_root <abs>  -dir_root <rel>",
         fc=COLOR_DECISION, fontsize=7.6)
    _arrow(ax, 2.5, 9.95, 2.5, 9.7, lw=1.0)

    _box(ax, 2.5, 9.2, 4.4, 0.85,
         "find <dir_root> -type d -name y_dir\n(任意深度, 自动去重排序)",
         fc=COLOR_CODE_BG, fontsize=7.2)
    _arrow(ax, 2.5, 8.78, 2.5, 8.5, lw=0.9)

    # tree of y_dir / <case>
    _box(ax, 2.5, 7.95, 4.4, 0.55,
         "y_dir/  (group A, B, ...)", fc=COLOR_NODE,
         fontsize=8, fontweight="bold")
    sub_cases = [
        (0.9, 6.95, "1_0.00"),
        (2.0, 6.95, "2_-0.01"),
        (3.1, 6.95, "3_+0.01"),
        (4.2, 6.95, "..."),
    ]
    for x, y, t in sub_cases:
        _box(ax, x, y, 0.95, 0.6, t, fc="#f0f7ff", ec="#9db8d8", fontsize=6.6)
        _arrow(ax, 2.5, 7.68, x, 7.25, lw=0.7, color="#7a8aa0")

    _box(ax, 2.5, 5.85, 4.4, 0.75,
         "汇总全部 y_dir 的一级子目录\n按 basename 排序后, -chunks K 分流",
         fc=COLOR_DECISION, fontsize=7.0)
    _arrow(ax, 2.5, 6.65, 2.5, 6.23, lw=0.9)

    _box(ax, 2.5, 4.85, 4.4, 0.95,
         "-lsubdir a b c   # 按 basename 过滤\n-dir_root .      # 默认以提交位置为根\n"
         "若 dir_root 本身是 y_dir 也纳入",
         fc=COLOR_CODE_BG, fontsize=6.9)

    # file layout (right half)
    _box(ax, 7.5, 10.3, 4.6, 0.7,
         "生成文件布局", fc=COLOR_NODE, fontsize=8.5, fontweight="bold")
    _arrow(ax, 7.5, 9.95, 7.5, 9.7, lw=1.0)

    _code(ax, 7.5, 7.9, 4.6, 3.4,
          "<path_root>/\n"
          "├─ <dir_root>/**/y_dir/<case>/\n"
          "│    └─ sub_slurm_univ.sh      # 每计算目录一脚本\n"
          "│       (parallel / each_subdir\n"
          "│        子作业)\n"
          "│\n"
          "└─ slurm/\n"
          "   ├─ sub_slurm_each_subdir_\n"
          "   │    chunk001.sh\n"
          "   ├─ sub_slurm_each_subdir_\n"
          "   │    chunk002.sh\n"
          "   ├─ ...\n"
          "   └─ sub_slurm_each_subdir_\n"
          "        parent.sh   # shared 且\n"
          "                     # chunks>1 时生成",
          fontsize=6.6)

    # bottom: file-roles strip
    roles = [
        (1.2, 2.6, "sub_slurm_univ.sh", "每目录计算脚本\n含 module / launcher / cmd"),
        (3.8, 2.6, "chunkNNN.sh", "each_subdir worker\n逐个 sbatch --wait"),
        (6.4, 2.6, "parent.sh", "shared 编排父作业\n-n1, 并发 K worker"),
        (8.9, 2.6, "single_alloc", "chunkNNN 即父作业\n内含 cmd 顺序执行"),
    ]
    for x, y, t, d in roles:
        _box(ax, x, y + 0.35, 2.3, 0.55, t, fc="#eef6ee", ec="#8fb88f",
             fontsize=7.0, fontweight="bold")
        _box(ax, x, y - 0.35, 2.3, 0.65, d, fc=COLOR_PANEL_BG, fontsize=6.6)
        _arrow(ax, x, y + 0.07, x, y - 0.02, lw=0.6, color="#7a8aa0")

    _box(ax, 5, 0.7, 9.6, 0.75,
         "parent.sh 只引用本次生成的 worker 列表, 不会误执行旧残留; "
         "worker 输出到 slurm-<parent-jobid>-chunkNNN.out",
         fc=COLOR_PANEL_BG, ec="#d1d5db", fontsize=7.0)

    _turn_off(ax)


# ---------------------------------------------------------------------------
# Panel 3: retry-strategy comparison table
# ---------------------------------------------------------------------------
def draw_panel_retry_table(ax):
    """Table comparing launch_retry vs sbatch_retry."""
    ax.set_title("面板3 · 重试策略对比 (launch_retry vs sbatch_retry)", loc="left")
    ax.axis("off")

    headers = ["维度", "pei_slurm_univ_launch_retry", "pei_slurm_univ_sbatch_retry"]
    rows = [
        ("包装对象", "srun / mpirun (并行启动器)", "sbatch (作业提交命令)"),
        ("重试的是什么",
         "launcher 自己没起来的暂态失败\n(程序根本没执行)",
         "sbatch 提交这一步没被 slurmctld 接住\n(根本没作业被创建)"),
        ("判定方式",
         "白名单: stderr 命中\nLAUNCH_FAIL_PATTERNS 才重试",
         "黑名单: 未命中 SBATCH_PERMANENT_\nPATTERNS 就重试"),
        ("成功判据",
         "exit 0",
         "exit 0 或输出含\n\"Submitted batch job\""),
        ("程序退出码",
         "原样透传 (非启动失败不重试)",
         "提交成功后原样透传\n(--wait 的作业结果不重投)"),
        ("默认次数 / 间隔",
         "99 次 × 10 s\n(空占分配, 不烧计算)",
         "99 次 × 10 s\n(不占分配, 只花提交进程时间)"),
        ("重试代价",
         "握着整个分配空转\n故必须严格只在启动失败时重试",
         "重试绝对安全\n不产生重复作业, 不占任何计算"),
        ("可叠加",
         "可与 sbatch_retry 叠加\n(两者职责不重叠)",
         "可与 launch_retry 叠加\n(两者职责不重叠)"),
    ]

    col_widths = [0.16, 0.42, 0.42]
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
    table.set_fontsize(7.4)
    table.scale(1.0, 2.0)

    # header
    for j in range(n_cols):
        cell = table[0, j]
        cell.set_facecolor(COLOR_HEADER_BG)
        cell.set_text_props(color="white", fontweight="bold", fontsize=8.2)
        cell.set_edgecolor("white")

    # body: first column bold, color the two script columns
    launch_col = "#f0f7ff"
    sbatch_col = "#fff7ed"
    for i in range(1, n_rows):
        c0 = table[i, 0]
        c0.set_text_props(ha="left", color=COLOR_TEXT, fontweight="bold", fontsize=7.2)
        c0.PAD = 0.04
        c1 = table[i, 1]
        c1.set_text_props(ha="left", color=COLOR_TEXT, fontsize=7.0)
        c1.PAD = 0.04
        c2 = table[i, 2]
        c2.set_text_props(ha="left", color=COLOR_TEXT, fontsize=7.0)
        c2.PAD = 0.04
        if i % 2 == 0:
            c0.set_facecolor("#f7f9fb")
            c1.set_facecolor(launch_col)
            c2.set_facecolor(sbatch_col)
        else:
            c1.set_facecolor("#fbfdff")
            c2.set_facecolor("#fffaf3")

    # footer note
    ax.text(0.5, 0.02,
            "注: launch_retry 的特征串只收 launcher 的「终态」串; sbatch_retry 默认偏向「都重试」, "
            "仅挡永久性配置错误 (Invalid partition / Unable to open file 等)。",
            ha="center", va="bottom", fontsize=6.8, color=COLOR_SKIP,
            transform=ax.transAxes)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


# ---------------------------------------------------------------------------
# Panel 4: preset registry
# ---------------------------------------------------------------------------
def draw_panel_preset_table(ax):
    """Table of submit.py presets."""
    ax.set_title("面板4 · pei_slurm_univ_submit.py preset 注册表", loc="left")
    ax.axis("off")

    headers = ["preset 名", "module_profile_type", "launcher", "cmd", "资源 (N/n)", "用途"]
    rows = [
        ("zcm6_vasp_0", "zcm6_vasp_0", "srun", "pei_vasp_univ_sbatch",
         "1 / 128", "VASP 5.4.4, 经 MY_LAUNCHER 传 srun"),
        ("zcm6_lammps_0", "zcm6_lammps_0", "srun", "lmp -in lmp.in",
         "1 / 24", "LAMMPS (自己编译版)"),
        ("zcm6_lammps_1", "zcm6_lammps_1", "srun", "lmp -in lmp.in",
         "1 / 24", "LAMMPS (nc 编译版)"),
        ("zcm6_n2p2_train_0", "zcm6_n2p2_0", "mpirun", "nnp-train",
         "1 / 24", "n2p2 训练, 推荐 16-32 核"),
        ("zcm6_n2p2_scaling_0", "zcm6_n2p2_0", "mpirun", "nnp-scaling 10000",
         "1 / 1", "n2p2 scaling, 单核估并行扩展性"),
    ]

    col_widths = [0.16, 0.16, 0.09, 0.22, 0.09, 0.28]
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
    table.set_fontsize(7.2)
    table.scale(1.0, 1.9)

    # header
    for j in range(n_cols):
        cell = table[0, j]
        cell.set_facecolor(COLOR_HEADER_BG)
        cell.set_text_props(color="white", fontweight="bold", fontsize=7.8)
        cell.set_edgecolor("white")

    # body
    for i in range(1, n_rows):
        c0 = table[i, 0]
        c0.set_text_props(ha="left", color=COLOR_TEXT, fontweight="bold", fontsize=7.0)
        c0.PAD = 0.04
        for j in range(1, n_cols):
            cell = table[i, j]
            cell.set_text_props(ha="left", color=COLOR_TEXT, fontsize=6.9)
            cell.PAD = 0.04
        if i % 2 == 0:
            for j in range(n_cols):
                table[i, j].set_facecolor("#f7f9fb")

    # footer note
    ax.text(0.5, 0.02,
            "所有 preset 默认 mode=each_subdir, chunks=5, partition=amd_512, dir_root=. ; "
            "命令行显式值始终覆盖 preset。-list_presets / -show_preset NAME 只读查看。",
            ha="center", va="bottom", fontsize=6.8, color=COLOR_SKIP,
            transform=ax.transAxes)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


# ---------------------------------------------------------------------------
# Panel 5: monitor dispatch + useful_command options
# ---------------------------------------------------------------------------
def draw_panel_monitor_commands(ax):
    """Engine dispatch of monitor_error + useful_command.sh options."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 11)
    ax.set_title("面板5 · 监控引擎调度 + useful_command.sh 选项", loc="left")

    # ---- left: monitor_error dispatch ----
    _box(ax, 2.5, 10.3, 4.4, 0.7,
         "pei_slurm_univ_monitor_error",
         fc=COLOR_HEADER_BG, text_color="white", fontweight="bold", fontsize=8.5)
    _arrow(ax, 2.5, 9.95, 2.5, 9.7, lw=1.0)

    _box(ax, 2.5, 9.2, 4.4, 0.85,
         "按 -<engine> 选监控脚本\n生成 1 节点 1 核看守作业\nsbatch 到 \$PEI_MONITOR_DIR",
         fc=COLOR_CODE_BG, fontsize=7.2)
    _arrow(ax, 2.5, 8.78, 2.5, 8.5, lw=0.9)

    # engine dispatch table
    _box(ax, 2.5, 7.95, 4.4, 0.55, "引擎调度 (case 分支)",
         fc=COLOR_NODE, fontsize=7.6, fontweight="bold")

    _code(ax, 2.5, 6.7, 4.4, 1.9,
          "-vasp)\n"
          "  monitor_cmd=pei_vasp_univ_monitor_error\n"
          "  ;;\n"
          "*)\n"
          "  echo \"Option $engine not recognized\"\n"
          "  exit 1\n"
          "  ;;\n"
          "# 新引擎: 加一行 case 即可",
          fontsize=6.8)
    _arrow(ax, 2.5, 7.68, 2.5, 7.65, lw=0.7)

    # generated watchdog script
    _code(ax, 2.5, 4.7, 4.4, 2.5,
          "# 生成 sub_pei_vasp_univ_monitor_error.sh\n"
          "#SBATCH -p amd_512\n"
          "#SBATCH -N 1  -n 1\n"
          "#SBATCH -J pei_monitor\n"
          "srun pei_vasp_univ_monitor_error \\\n"
          "  -skip_ljobid $SLURM_JOB_ID $@\n"
          "# (经 sbatch_retry 提交, 容错暂态繁忙)",
          fontsize=6.8)
    _arrow(ax, 2.5, 6.0, 2.5, 5.95, lw=0.7, color="#7a8aa0")

    # env knobs
    _box(ax, 2.5, 3.0, 4.4, 1.1,
         "环境变量可调:\n"
         "PEI_MONITOR_PARTITION  (默认 amd_512)\n"
         "PEI_MONITOR_JOBNAME    (默认 pei_monitor)\n"
         "PEI_MONITOR_DIR        (默认 ./slurm_monitor)",
         fc="#fff7ed", ec="#d97706", fontsize=6.9)
    _arrow(ax, 2.5, 3.45, 2.5, 3.55, lw=0.7, color="#7a8aa0")

    _box(ax, 2.5, 1.7, 4.4, 0.85,
         "传透选项: -phase_check, -cancel_on_phase_change,\n"
         "-no_phase_check ... (透给引擎监控器)",
         fc=COLOR_PANEL_BG, fontsize=6.8)
    _arrow(ax, 2.5, 2.45, 2.5, 2.13, lw=0.7, color="#7a8aa0")

    # ---- right: useful_command.sh options ----
    _box(ax, 7.5, 10.3, 4.4, 0.7,
         "pei_slurm_univ_useful_command.sh",
         fc=COLOR_HEADER_BG, text_color="white", fontweight="bold", fontsize=8.2)
    _arrow(ax, 7.5, 9.95, 7.5, 9.7, lw=1.0)

    _box(ax, 7.5, 9.2, 4.4, 0.85,
         "无参数 → 打印命令速查表\n(Slurm / VASP / editor / network)\n带参数 → 实际探测网络可达性",
         fc=COLOR_CODE_BG, fontsize=7.0)
    _arrow(ax, 7.5, 8.78, 7.5, 8.5, lw=0.9)

    # options table
    opts = [
        ("-openai", "探测 OpenAI 端点"),
        ("-claude", "探测 Anthropic / Claude 端点"),
        ("-github", "探测 GitHub 端点"),
        ("-pypi", "探测 PyPI 索引 (含清华镜像)"),
        ("-conda", "探测 Anaconda / conda-forge"),
        ("-net", "探测以上全部组"),
        ("-ip", "显示公网 IP"),
        ("-proxy", "打印代理环境变量"),
        ("-monitor", "汇总 RUNNING 作业 / 路径 / 输出"),
        ("-commands", "打印速查表"),
        ("-all", "-net -ip -proxy -commands"),
        ("-timeout SEC", "单请求超时 (默认 8)"),
    ]
    y0 = 7.95
    for opt, desc in opts:
        _box(ax, 6.5, y0, 1.8, 0.42, opt,
             fc="#eef6ee", ec="#8fb88f", fontsize=6.8, fontweight="bold")
        _box(ax, 8.55, y0, 2.9, 0.42, desc,
             fc=COLOR_PANEL_BG, fontsize=6.7)
        y0 -= 0.55

    # bottom: n2p2 + vasp convenience launchers note
    _box(ax, 5, 0.7, 9.6, 0.75,
         "slurm_vasp: pei_slurm_univ_vasp_monitor = sbatch pei_slurm_univ_vasp_monitor.sh  |  "
         "slurm_n2p2: n2p2_train (amd_512, 24核, nnp-train) / n2p2_scaling (amd_512, 1核, nnp-scaling 10000)",
         fc=COLOR_DECISION, fontsize=6.9)

    _turn_off(ax)


# ---------------------------------------------------------------------------
# Figure assembly
# ---------------------------------------------------------------------------
def render_overview(path_image: Path) -> None:
    """Build the 5-panel overview figure and save it to path_image."""
    _set_rcparams()
    fig = plt.figure(figsize=(16, 12), dpi=150)
    gs = fig.add_gridspec(
        3, 2, hspace=0.34, wspace=0.12,
        left=0.03, right=0.98, top=0.95, bottom=0.03,
    )
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    ax5 = fig.add_subplot(gs[2, :])

    draw_panel_modes(ax1)
    draw_panel_discovery_layout(ax2)
    draw_panel_retry_table(ax3)
    draw_panel_preset_table(ax4)
    draw_panel_monitor_commands(ax5)

    fig.suptitle(
        "slurm_utils 子包功能总览",
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
    path_image = path_output / "slurm_utils_overview.png"

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
        description="Render the slurm_utils overview figure (sbatch-free).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("docs/_build/example-slurm-utils"),
        help="Output directory for the rendered PNG (default: docs/_build/example-slurm-utils)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary = run_example(args.output)
    # deterministic assertions
    assert summary["size_kb"] > 10.0, "PNG smaller than 10 KB: " + str(summary)
    assert summary["pixel_deviation"] > 0.002, "PNG effectively blank: " + str(summary)
    print("OK: slurm_utils overview figure verified non-blank.")


if __name__ == "__main__":
    main()
