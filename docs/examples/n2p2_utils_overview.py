#!/usr/bin/env python3
"""Overview figure for the n2p2_utils sub-package.

This tutorial is n2p2-free. It demonstrates the core functionality of the
n2p2_utils sub-package (Python scripts + bash scripts + n2p2_universal/) by
rendering a single multi-panel figure that summarises:

  * the end-to-end n2p2 workflow from VASP data to a trained potential,
  * the active-learning symmetry-function (SF) selection sub-flow,
  * the parameter table of sfs_gen_basic_SF.py (G2/G4, -G/-e/-c/-n/-z),
  * the pei_n2p2_univ_run pipeline (nnp-scaling -> nnp-train) with its
    0 / 1 / 10 exit-code contract and done-file sentinel, and
  * the two-stage clean-up strategy of pei_n2p2_univ_clean_train.

No n2p2 binary, no ``mpirun``, no real training is invoked. The figure is
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
COLOR_DATA = "#dbeafe"
COLOR_DATA_BORDER = "#2563eb"
COLOR_TRAIN = "#ede9fe"
COLOR_TRAIN_BORDER = "#7c3aed"
COLOR_ACTIVE = "#fef3c7"
COLOR_ACTIVE_BORDER = "#d97706"
COLOR_EXIT10 = "#dbeafe"
COLOR_EXIT0 = "#dcfce7"
COLOR_EXIT1 = "#fee2e2"
COLOR_HEADER_BG = "#2c3e50"


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
            # Microsoft YaHei, while Latin / - / digits fall back to DejaVu Sans.
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
# Panel 1: n2p2 end-to-end workflow
# ---------------------------------------------------------------------------
def draw_panel_workflow(ax):
    """Full pipeline: VASP OUTCAR -> input.data -> sflist -> active SF -> train."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 11)
    ax.set_title("面板1 · n2p2 完整工作流 (VASP 数据 -> 训练势函数)", loc="left")

    cx = 5.0

    # Stage A: VASP data -> input.data
    _box(ax, 1.6, 10.1, 2.6, 0.8,
         "VASP OUTCAR / CONTCAR\n(数据集目录)",
         fc=COLOR_DATA, ec=COLOR_DATA_BORDER, fontsize=8)
    _box(ax, cx, 10.1, 3.0, 0.8,
         "data_read.py + data_input.py\n(VASPReader -> input.data)\n+ energy.npy / forces.npy",
         fc=COLOR_CODE_BG, fontsize=7.4)
    _box(ax, 8.4, 10.1, 2.6, 0.8,
         "input.data\n(n2p2 结构+标签)",
         fc=COLOR_NODE, fontsize=8, fontweight="bold")
    _arrow(ax, 2.9, 10.1, 3.5, 10.1)
    _arrow(ax, 6.5, 10.1, 7.1, 10.1)

    # Stage B: SF generation -> input.nn
    _box(ax, 1.6, 8.5, 2.6, 0.9,
         "sfs_gen_basic_SF.py\n-G G2/G4 -e Mg -c 6.2\n-n 5 -z 1,4,16",
         fc=COLOR_CODE_BG, fontsize=7.2)
    _box(ax, cx, 8.5, 2.6, 0.9,
         "sflist\n(对称函数参数行)",
         fc=COLOR_NODE, fontsize=8, fontweight="bold")
    _box(ax, 8.4, 8.5, 2.6, 0.9,
         "input.nn (+ input.nn.all)\n(n2p2 控制文件)",
         fc=COLOR_NODE, fontsize=8, fontweight="bold")
    _arrow(ax, 2.9, 8.5, 3.7, 8.5, label=">> sflist", label_offset=(0, 0.2))
    _arrow(ax, 6.3, 8.5, 7.1, 8.5)
    _arrow(ax, 8.4, 9.7, 8.4, 8.95, lw=1.0, color="#7a8aa0")
    _arrow(ax, 1.6, 9.7, 1.6, 8.95, lw=1.0, color="#7a8aa0")

    # Stage C: active learning SF selection (sub-flow, see panel 2)
    _box(ax, cx, 6.9, 8.6, 1.1,
         "active_sf_0 子流程 (见面板2):\n"
         "sub_cal.py (分批) -> collect_sf.py (收集 feat_atom/feat_av) -> "
         "select_feat.py (CUR/FPS 选 N 个 SF) -> 精选 input.nn",
         fc=COLOR_ACTIVE, ec=COLOR_ACTIVE_BORDER, fontsize=7.6)
    _arrow(ax, cx, 8.05, cx, 7.45, label="input.nn.all", label_offset=(0.7, 0))

    # Stage D: training
    _box(ax, 1.6, 5.1, 3.4, 0.9,
         "data_get_data.bash\n(read -> input -> 复制到\ntrain/ 和 active_sf_0/)",
         fc=COLOR_CODE_BG, fontsize=7.0)
    _box(ax, cx, 5.1, 2.8, 0.9,
         "train_train.bash\nmpirun -np 16 nnp-train",
         fc=COLOR_RUN, fontsize=7.6, fontweight="bold")
    _box(ax, 8.4, 5.1, 3.0, 0.9,
         "learning-curve.out\nscaling.data\nweights.*.out",
         fc=COLOR_NODE, fontsize=7.4, fontweight="bold")
    _arrow(ax, 3.3, 5.1, 3.6, 5.1)
    _arrow(ax, 6.4, 5.1, 6.9, 5.1)
    _arrow(ax, cx, 6.35, cx, 5.55, label="精选 input.nn", label_offset=(0.9, 0))

    # Stage E: results
    _box(ax, 1.6, 3.5, 3.4, 0.9,
         "train_get_result.bash\n复制 input.nn / scaling.data\nlearning-curve.out -> result/post/",
         fc=COLOR_CODE_BG, fontsize=7.0)
    _box(ax, cx, 3.5, 2.8, 0.9,
         "train_final_result.bash\n读 epoch.txt -> 提取\nweights / trainpoints",
         fc=COLOR_CODE_BG, fontsize=7.0)
    _box(ax, 8.4, 3.5, 3.0, 0.9,
         "nnp-data/ (势函数)\n+ test 评估 (LAMMPS)",
         fc=COLOR_TRAIN, ec=COLOR_TRAIN_BORDER, fontsize=7.4, fontweight="bold")
    _arrow(ax, 8.4, 4.65, 8.4, 3.95, lw=1.0, color="#7a8aa0")
    _arrow(ax, 3.3, 3.5, 3.6, 3.5)
    _arrow(ax, 6.4, 3.5, 6.9, 3.5)

    # n2p2_universal runner band
    _box(ax, cx, 1.9, 9.0, 1.1,
         "n2p2_universal/ 统一运行层:\n"
         "pei_n2p2_univ_load_env (source, 加载 eigen/gsl/openmpi module, nnp-* 入 PATH)\n"
         "pei_n2p2_univ_run (nnp-scaling -> nnp-train, 哨兵文件, 退出码 0/10/1, 见面板4)\n"
         "pei_n2p2_univ_clean_train (两阶段清理, 见面板5)",
         fc=COLOR_PANEL_BG, ec="#d1d5db", fontsize=7.4)
    _arrow(ax, cx, 3.05, cx, 2.45, lw=1.0, color="#7a8aa0")

    # bottom summary
    _box(ax, cx, 0.5, 9.4, 0.7,
         "完整链路: VASP 数据 -> input.data -> 对称函数 -> 精选 SF -> nnp-train -> 势函数 + LAMMPS 评估",
         fc=COLOR_DECISION, fontsize=7.8)
    _arrow(ax, cx, 1.35, cx, 0.85, lw=1.0, color="#7a8aa0")

    _turn_off(ax)


# ---------------------------------------------------------------------------
# Panel 2: active learning SF selection sub-flow
# ---------------------------------------------------------------------------
def draw_panel_active_sf(ax):
    """Flowchart of the active-learning SF selection sub-flow."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 11)
    ax.set_title("面板2 · active learning 对称函数选择流程", loc="left")

    cx = 5.0

    # input.nn.all
    _box(ax, cx, 10.3, 4.6, 0.7,
         "input.nn.all  (全部候选 SF)",
         fc=COLOR_NODE, fontsize=8.5, fontweight="bold")

    # sub_cal: split + batch
    _box(ax, cx, 9.0, 7.4, 1.0,
         "active_sf_0_sub_cal.py\n"
         "分离 symfunction_short 行 -> input.nn.no (无 SF)\n"
         "逐条 SF 写入 input.nn, 复制到 work_dir/0000, 0001, ... 提交",
         fc=COLOR_CODE_BG, fontsize=7.2)
    _arrow(ax, cx, 9.95, cx, 9.5, label="分离 + 分批", label_offset=(1.0, 0))

    # work_dir batch outputs
    _box(ax, 1.8, 7.6, 2.4, 0.9,
         "work_dir/0000/\nfunction.data\n(单条 SF 值)",
         fc="#f0f7ff", ec="#9db8d8", fontsize=7.0)
    _box(ax, cx, 7.6, 2.4, 0.9,
         "work_dir/0001/\nfunction.data",
         fc="#f0f7ff", ec="#9db8d8", fontsize=7.0)
    _box(ax, 8.2, 7.6, 2.4, 0.9,
         "work_dir/NNNN/\nfunction.data\n(每条 SF 一批)",
         fc="#f0f7ff", ec="#9db8d8", fontsize=7.0)
    _arrow(ax, 3.0, 8.5, 2.0, 8.05, lw=0.9, color="#7a8aa0")
    _arrow(ax, cx, 8.5, cx, 8.05, lw=0.9, color="#7a8aa0")
    _arrow(ax, 7.0, 8.5, 8.0, 8.05, lw=0.9, color="#7a8aa0")

    # collect_sf
    _box(ax, cx, 6.1, 8.6, 1.1,
         "active_sf_0_collect_sf.py\n"
         "遍历 work_dir/*/function.data -> 每条 SF 的\n"
         "feat_atom (逐原子最大值) + feat_av (逐帧平均值)\n"
         "累计拼接 -> feat_atom.npy / feat_av.npy",
         fc=COLOR_CODE_BG, fontsize=7.2)
    _arrow(ax, 1.8, 7.15, 3.5, 6.65, lw=0.9, color="#7a8aa0")
    _arrow(ax, cx, 7.15, cx, 6.65, lw=0.9, color="#7a8aa0")
    _arrow(ax, 8.2, 7.15, 6.5, 6.65, lw=0.9, color="#7a8aa0")

    # select_feat
    _box(ax, cx, 4.4, 8.6, 1.3,
         "active_sf_0_select_feat.py  (skcosmo)\n"
         "1. 删除恒零列 (0 值占比 > max_0=5%)\n"
         "2. CUR (n_to_select=N) 或 FPS 特征选择\n"
         "3. 把选中的 symfunction_short 行写回 input.nn",
         fc=COLOR_ACTIVE, ec=COLOR_ACTIVE_BORDER, fontsize=7.4)
    _arrow(ax, cx, 5.55, cx, 5.05, label="feat_atom / feat_av", label_offset=(1.6, 0))

    # selected output
    _box(ax, 2.6, 2.8, 3.6, 0.9,
         "input.nn  (精选 N 个 SF)\n+ SFs_all.dat / SFs_drop.dat",
         fc=COLOR_NODE, fontsize=7.6, fontweight="bold")
    _box(ax, 7.4, 2.8, 3.6, 0.9,
         "active_sf_0_select.bash\n(collect_sf.py -> select_feat.py N\nN 默认 36/48)",
         fc=COLOR_CODE_BG, fontsize=7.0)
    _arrow(ax, 3.8, 3.75, 2.6, 3.25, lw=1.0, color="#7a8aa0")
    _arrow(ax, cx, 3.75, 7.4, 3.25, lw=1.0, color="#7a8aa0")

    # downstream
    _box(ax, cx, 1.3, 8.6, 1.1,
         "下游: 精选 input.nn 复制到 train/ -> train_train.bash 训练\n"
         "(active_sf_0_select.bash 末尾: cp input.nn ../train/)",
         fc=COLOR_RUN, fontsize=7.6)
    _arrow(ax, 2.6, 2.35, 3.5, 1.85, lw=1.0, color="#7a8aa0")
    _arrow(ax, 7.4, 2.35, 6.5, 1.85, lw=1.0, color="#7a8aa0")

    # note
    _box(ax, cx, 0.3, 9.4, 0.5,
         "注: select_feat.py 还含 SF._explore() 邻域探索 (eta*/2, rs+-0.1*rc, rc*/1.26), "
         "当前默认注释关闭",
         fc=COLOR_PANEL_BG, ec="#d1d5db", fontsize=7.0)

    _turn_off(ax)


# ---------------------------------------------------------------------------
# Panel 3: sfs_gen_basic_SF.py parameter table + G2 vs G4
# ---------------------------------------------------------------------------
def draw_panel_sf_params(ax):
    """Parameter table for sfs_gen_basic_SF.py and G2 vs G4 comparison."""
    ax.set_title("面板3 · sfs_gen_basic_SF.py 参数表与 G2/G4 对比", loc="left")
    ax.axis("off")

    # ---- left: parameter table ----
    headers = ["参数", "类型/默认", "含义", "示例"]
    rows = [
        ("-G", "str (G2)", "对称函数形式\n(G2=径向, G4=角向)", "-G G2"),
        ("-e", "str (H)", "元素列表, 逗号分隔\n(支持 H..No 全元素校验)", "-e Mg  /  -e Mg,Al"),
        ("-c", "float (12.0)", "截断半径 Rc\n(单条命令一个 cutoff)", "-c 6.2"),
        ("-n", "int (20)", "区间数 N, 生成 2N+1 个\neta/rs 采样点", "-n 5  (11 点)"),
        ("-z", "str (1,4,16)", "角向 SF 的 zeta 列表\n(逗号分隔)", "-z 1,4,16"),
    ]

    col_widths = [0.10, 0.18, 0.46, 0.26]
    n_cols = len(headers)
    n_rows = len(rows) + 1
    table = ax.table(
        cellText=rows,
        colLabels=headers,
        colWidths=col_widths,
        cellLoc="left",
        loc="upper left",
        bbox=[0.0, 0.30, 0.62, 0.68],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(7.6)
    table.scale(1.0, 1.7)

    for j in range(n_cols):
        cell = table[0, j]
        cell.set_facecolor(COLOR_HEADER_BG)
        cell.set_text_props(color="white", fontweight="bold", fontsize=8)
        cell.set_edgecolor("white")

    for i in range(1, n_rows):
        cell0 = table[i, 0]
        cell0.set_text_props(ha="center", color=COLOR_TEXT, fontweight="bold",
                             fontsize=8.5,
                             family=["Microsoft YaHei", "Noto Sans SC", "SimHei", "DejaVu Sans Mono"])
        for j in range(1, n_cols):
            cell = table[i, j]
            cell.set_text_props(color=COLOR_TEXT, fontsize=7.2)
            if j == 3:
                cell.set_text_props(color=COLOR_WARN, fontweight="bold", fontsize=7.0,
                                    family=["Microsoft YaHei", "Noto Sans SC", "SimHei", "DejaVu Sans Mono"])
            if i % 2 == 0:
                cell.set_facecolor("#f7f9fb")
        if i % 2 == 0:
            cell0.set_facecolor("#f7f9fb")

    # ---- right: G2 vs G4 comparison ----
    ax.text(0.66, 0.95, "G2 vs G4 对比", ha="left", va="top",
            fontsize=10, fontweight="bold", color=COLOR_TEXT,
            transform=ax.transAxes)

    # G2 box
    g2_text = (
        "G2 (径向, type 2)\n"
        "symfunction_short <el> 2 <el> <eta> <rs> <rc>\n\n"
        "  G2 = exp(-eta*(Rij-Rs)^2) * fc(Rij)\n\n"
        "  - 描述中心原子与单个邻居的距离\n"
        "  - 2N+1 个 eta, 每个 eta 配 rs=0 和 rs=cutoff/2\n"
        "  - 对每对 (fel, sel) 各生成一组"
    )
    _box_axes(ax, 0.80, 0.72, 0.32, 0.40, g2_text,
              fc="#dbeafe", ec="#2563eb", fontsize=7.2)

    # G4 box (note: source uses G3 label but task spec says G4; n2p2 angular=type 3)
    g4_text = (
        "G4 (角向, type 3)\n"
        "symfunction_short <el> 3 <el> <el> <eta> <lambda> <zeta> <rc> <rs>\n\n"
        "  G4 = 2^(1-zeta) * (1+lambda*cos(theta))**zeta\n"
        "        * exp(-eta*Rij^2) * fc\n\n"
        "  - 描述三体夹角 theta(ijk)\n"
        "  - 每个 eta x 每个 zeta x lambda=+1/-1\n"
        "  - 去重: ang(a,a,b)==ang(a,b,a)"
    )
    _box_axes(ax, 0.80, 0.30, 0.32, 0.40, g4_text,
              fc="#ede9fe", ec="#7c3aed", fontsize=7.0)

    # bottom note
    ax.text(0.02, 0.06,
            "调用示例 (sfs_gen_gen_sf.bash):\n"
            "  basic_SF.py -G G2 -e Mg -c 6.2 -n 5  >> sflist\n"
            "  basic_SF.py -G G4 -e Mg -c 6.2 -n 5 -z 1,4,16  >> sflist",
            ha="left", va="bottom", fontsize=7.0, color=COLOR_TEXT,
            transform=ax.transAxes,
            family=["Microsoft YaHei", "Noto Sans SC", "SimHei", "DejaVu Sans Mono", "DejaVu Sans"],
            bbox=dict(boxstyle="round,pad=0.3", fc=COLOR_CODE_BG, ec="#d1d5db", lw=0.8))

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


def _box_axes(ax, cx, cy, w, h, text, fc=COLOR_NODE, ec=COLOR_NODE_BORDER,
              fontsize=8):
    """Draw a box using axes-fraction coordinates (for table-panel hybrids).

    CJK fonts lead the family list so Chinese glyphs resolve before the
    mono fallback (DejaVu Sans Mono lacks CJK).
    """
    box = mpatches.FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle="round,pad=0.01", fc=fc, ec=ec, lw=1.1,
        transform=ax.transAxes, zorder=3,
    )
    ax.add_patch(box)
    ax.text(cx, cy, text, ha="center", va="center",
            fontsize=fontsize, color=COLOR_TEXT, zorder=4,
            transform=ax.transAxes,
            family=["Microsoft YaHei", "Noto Sans SC", "SimHei", "DejaVu Sans Mono", "DejaVu Sans"])


# ---------------------------------------------------------------------------
# Panel 4: pei_n2p2_univ_run pipeline
# ---------------------------------------------------------------------------
def draw_panel_univ_run(ax):
    """Flowchart of pei_n2p2_univ_run: nnp-scaling -> nnp-train + exit codes."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 11)
    ax.set_title("面板4 · pei_n2p2_univ_run pipeline (nnp-scaling -> nnp-train)", loc="left")

    cx = 5.0

    # entry
    _box(ax, cx, 10.3, 6.6, 0.7,
         "pei_n2p2_univ_run <job_dir> [scaling_bins] [done_file] [done_marker]",
         fc=COLOR_HEADER_BG, text_color="white", fontweight="bold", fontsize=8)
    _code(ax, cx, 9.5, 6.4, 0.55,
          "cd $job_dir  (默认 .)   |   scaling_bins=10000   |   done_file=n2p2.done",
          fontsize=7.0)

    # done_file check
    _box(ax, cx, 8.5, 6.6, 0.9,
         "done_file 存在且含 done_marker ?\n(哨兵文件机制: n2p2 无原生完成标记,\n本 runner 自写 sentinel)",
         fc=COLOR_DECISION, ec=COLOR_DECISION_BORDER, fontsize=7.4)
    _box(ax, 1.6, 8.5, 2.4, 0.9,
         "exit 10\n已完成, 跳过",
         fc=COLOR_EXIT10, fontsize=8.5, fontweight="bold")
    _arrow(ax, 3.4, 8.5, 2.8, 8.5, label="是 Yes", label_offset=(0, 0.2),
           label_color=COLOR_OK)

    # input check
    _box(ax, cx, 7.1, 6.6, 0.8,
         "校验 input.nn + input.data 均存在\n(缺失 -> exit 1)",
         fc=COLOR_NODE, fontsize=7.6)
    _arrow(ax, cx, 8.05, cx, 7.5, label="否 No", label_offset=(0.7, 0))

    # Stage 1: scaling
    _box(ax, 1.6, 5.7, 3.4, 1.0,
         "scaling.data 非空 ?\n是 -> 跳过 nnp-scaling\n否 -> srun nnp-scaling $bins",
         fc=COLOR_DECISION, ec=COLOR_DECISION_BORDER, fontsize=7.2)
    _box(ax, 6.0, 5.7, 3.6, 1.0,
         "Stage 1: nnp-scaling\n(srun nnp-scaling 10000)\n产出 scaling.data\n失败 -> exit 1",
         fc=COLOR_RUN, fontsize=7.4)
    _arrow(ax, 3.3, 5.7, 4.2, 5.7, label="需运行", label_offset=(0, 0.2))
    _arrow(ax, cx, 6.7, 1.6, 6.2, lw=1.0, color="#7a8aa0")
    _arrow(ax, cx, 6.7, 6.0, 6.2, lw=1.0, color="#7a8aa0")

    # Stage 2: training
    _box(ax, cx, 4.0, 6.6, 0.9,
         "Stage 2: nnp-train\nsrun nnp-train  (无 checkpoint 续训,\n残留输出需先 clean 再重跑)",
         fc=COLOR_RUN, fontsize=7.4, fontweight="bold")
    _arrow(ax, 6.0, 5.2, cx, 4.45, lw=1.0, color="#7a8aa0")
    _arrow(ax, 1.6, 5.2, cx, 4.45, lw=1.0, color="#7a8aa0",
           label="已有 scaling.data", label_offset=(-1.8, 0.1))

    # write sentinel + exit 0
    _box(ax, cx, 2.6, 6.6, 0.9,
         '写 done_marker -> done_file\n("n2p2 pipeline finished")',
         fc=COLOR_NODE, fontsize=7.6)
    _box(ax, 1.6, 2.6, 2.4, 0.9,
         "exit 0\n本轮完成",
         fc=COLOR_EXIT0, fontsize=8.5, fontweight="bold")
    _arrow(ax, 3.4, 2.6, 2.8, 2.6)
    _arrow(ax, cx, 3.55, cx, 3.05, label="nnp-train 成功", label_offset=(1.4, 0),
           label_color=COLOR_OK)

    # failure path
    _box(ax, 8.4, 4.0, 2.6, 0.9,
         "exit 1\nnnp-scaling / nnp-train\n失败或缺输入",
         fc=COLOR_EXIT1, fontsize=8, fontweight="bold")
    _arrow(ax, 7.6, 4.0, 7.1, 4.0, label="失败", label_offset=(0, 0.2),
           label_color=COLOR_ERR)

    # load_env strip
    _box(ax, cx, 1.1, 9.4, 1.0,
         "pei_n2p2_univ_load_env (source, 非 exec):\n"
         "  source module.sh -> load eigen/3.8.8 gsl/2.5 openmpi/2.0.4\n"
         "  export PATH=$N2P2_BIN_DIR:$PATH  (nnp-train / nnp-scaling 入 PATH)\n"
         "  pei_n2p2_check_inputs input.nn ...  (校验输入文件)",
         fc=COLOR_PANEL_BG, ec="#d1d5db", fontsize=7.0)
    _arrow(ax, cx, 2.15, cx, 1.6, lw=1.0, color="#7a8aa0")

    # legend
    legend_items = [
        ("exit 0 = 完成", COLOR_EXIT0),
        ("exit 10 = 已完成跳过", COLOR_EXIT10),
        ("exit 1 = 失败", COLOR_EXIT1),
        ("判定/哨兵", COLOR_DECISION),
    ]
    for i, (txt, col) in enumerate(legend_items):
        lx = 0.3 + i * 2.5
        ax.add_patch(mpatches.Rectangle((lx, 0.15), 0.24, 0.18, fc=col,
                                        ec=COLOR_NODE_BORDER, lw=0.6, zorder=3))
        ax.text(lx + 0.32, 0.24, txt, ha="left", va="center", fontsize=7,
                color=COLOR_TEXT, zorder=4)

    _turn_off(ax)


# ---------------------------------------------------------------------------
# Panel 5: pei_n2p2_univ_clean_train two-stage strategy
# ---------------------------------------------------------------------------
def draw_panel_clean_train(ax):
    """Table of pei_n2p2_univ_clean_train two-stage clean-up strategy."""
    ax.set_title("面板5 · pei_n2p2_univ_clean_train 两阶段清理策略", loc="left")
    ax.axis("off")

    headers = ["阶段", "目标目录", "清理内容", "安全机制", "选项"]
    rows = [
        ("阶段 1\n(train)",
         "y_dir/<job>/\n(nnp-train 工作目录)",
         "裁剪按 epoch 编号的诊断族:\nneuron-stats / trainpoints /\ntestpoints / trainforces /\ntestforces\n每 50 步留 1 个末尾 epoch",
         "逐个明确路径 rm -- <file>\n绝不递归/通配\nweights.<elem>.<epoch>.out\n默认不动",
         "--skip-train 关闭\n--include-weights\n连 weights 也只留末尾"),
        ("阶段 2\n(post)",
         "y_post/*/properties/\ny_epoch_scan/y_dir/\n(逐 epoch LAMMPS 计算)",
         "回收整棵 y_dir 树\n(单棵可达 ~78k 文件 / 22GB)\n势函数与汇总数据已被\npost 复制/聚合出来",
         "epoch 覆盖校验:\np_post_epoch_{cij,gsfe,stretch}.txt\n齐备+非空+epoch 列覆盖全\n形状护栏: */y_post/*/.../y_dir\n非 epoch 子目录永拦",
         "--skip-post 关闭\n--force 跳过覆盖校验\n(形状护栏仍生效)"),
    ]

    col_widths = [0.10, 0.20, 0.30, 0.25, 0.15]
    n_cols = len(headers)
    n_rows = len(rows) + 1
    table = ax.table(
        cellText=rows,
        colLabels=headers,
        colWidths=col_widths,
        cellLoc="left",
        loc="upper center",
        bbox=[0.01, 0.22, 0.98, 0.68],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(7.4)
    table.scale(1.0, 2.6)

    for j in range(n_cols):
        cell = table[0, j]
        cell.set_facecolor(COLOR_HEADER_BG)
        cell.set_text_props(color="white", fontweight="bold", fontsize=8.5)
        cell.set_edgecolor("white")

    stage_colors = {1: "#dbeafe", 2: "#ede9fe"}
    for i in range(1, n_rows):
        cell0 = table[i, 0]
        cell0.set_facecolor(stage_colors[i])
        cell0.set_text_props(ha="center", color=COLOR_TEXT, fontweight="bold", fontsize=8)
        for j in range(1, n_cols):
            cell = table[i, j]
            cell.set_text_props(color=COLOR_TEXT, fontsize=7.0)
            if i % 2 == 0:
                cell.set_facecolor("#f7f9fb")

    # options strip
    ax.text(0.5, 0.14, "命令行选项", ha="center", va="center",
            fontsize=9, fontweight="bold", color=COLOR_TEXT,
            transform=ax.transAxes)
    opt_text = (
        "--apply            真正执行删除 (默认 dry-run 只预览)        --size   统计将释放空间 (慢)\n"
        "--skip-train       跳过阶段 1, 只处理 y_post                 --skip-post  跳过阶段 2, 只处理 y_dir\n"
        "--include-weights  阶段 1 同时裁剪 weights (默认全保留)       --force  阶段 2 跳过 epoch 覆盖校验强删"
    )
    ax.text(0.5, 0.075, opt_text, ha="center", va="center",
            fontsize=7.0, color=COLOR_TEXT, transform=ax.transAxes,
            family=["Microsoft YaHei", "Noto Sans SC", "SimHei", "DejaVu Sans Mono", "DejaVu Sans"],
            bbox=dict(boxstyle="round,pad=0.3", fc=COLOR_CODE_BG, ec="#d1d5db", lw=0.8))

    # exit code note
    ax.text(0.5, 0.015,
            "退出码: 0 = 全部处理完毕无待关注;  1 = 参数/结构错误, 或阶段 2 有 job 校验未过被拒删 (y_dir 原样保留)",
            ha="center", va="bottom", fontsize=7.0, color=COLOR_SKIP,
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
    # 3 rows x 2 cols, last cell (panel 5) spans bottom-full-width
    gs = fig.add_gridspec(
        3, 2, hspace=0.34, wspace=0.12,
        left=0.03, right=0.98, top=0.95, bottom=0.03,
    )
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    ax5 = fig.add_subplot(gs[2, :])

    draw_panel_workflow(ax1)
    draw_panel_active_sf(ax2)
    draw_panel_sf_params(ax3)
    draw_panel_univ_run(ax4)
    draw_panel_clean_train(ax5)

    # suptitle
    fig.suptitle(
        "n2p2_utils 子包功能总览",
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
    path_image = path_output / "n2p2_utils_overview.png"

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
        description="Render the n2p2_utils overview figure (n2p2-free).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("docs/_build/example-n2p2-utils"),
        help="Output directory for the rendered PNG (default: docs/_build/example-n2p2-utils)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary = run_example(args.output)
    # deterministic assertions
    assert summary["size_kb"] > 10.0, "PNG smaller than 10 KB: " + str(summary)
    assert summary["pixel_deviation"] > 0.002, "PNG effectively blank: " + str(summary)
    print("OK: n2p2_utils overview figure verified non-blank.")


if __name__ == "__main__":
    main()
