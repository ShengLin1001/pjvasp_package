#!/usr/bin/env python3
"""Overview figure for the vasp_utils/vasp_workflow_bulk + neb_utils sub-packages.

This tutorial is VASP-free. It renders a single multi-panel figure that summarises:

  * the 10 ``pei_vasp_run_*`` workflow scripts grouped into 5 functional classes,
  * the shared ``y_full_relax`` -> ``y_<workflow>/y_dir/<case>`` lifecycle,
  * the ``pei_vasp_plot_all`` dispatcher mapping (option -> y_* dir -> plot script),
  * the NEB migration-barrier tool chain across neb_utils' 8 perl/python/bash scripts,
  * the ``pei_vasp_run_stretch.py`` strain-type comparison and how it differs from
    ``pei_vasp_run_eos.py``'s V/V0 volume-ratio scan.

No VASP binary, ``sbatch``, ``nebmake.pl``, ``n2p2`` or LAMMPS is invoked.  Every
energy / strain / timing series below is synthetic (generated with numpy) and is
only meant to illustrate what each workflow produces.  The figure is generated
with matplotlib only.
"""

from __future__ import annotations

import argparse
import shutil
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
COLOR_DATA1 = "#2b6cb0"   # blue  – primary data series
COLOR_DATA2 = "#c05621"   # orange – secondary data series
COLOR_FIT = "#6b7280"     # grey   – fitted curve
COLOR_HIGHLIGHT = "#9b2c2c"

# Category colours for Panel 1 (5 functional classes)
CAT_VOLUME = "#dbeafe"        # volume/strain family
CAT_VOLUME_BORDER = "#3b82f6"
CAT_CONVERGENCE = "#dcfce7"   # convergence/benchmark family
CAT_CONVERGENCE_BORDER = "#22c55e"
CAT_SURFACE = "#fef3c7"       # surface/interface family
CAT_SURFACE_BORDER = "#d97706"
CAT_ELECTRON = "#ede9fe"      # electronic-structure family
CAT_ELECTRON_BORDER = "#8b5cf6"
CAT_NEB = "#fce7f3"           # transition-state family
CAT_NEB_BORDER = "#ec4899"

# Language tag colours for the NEB tool chain (Panel 4)
LANG_PERL = "#fef3c7"      # perl  – warm yellow
LANG_PERL_BORDER = "#d97706"
LANG_PY = "#dcfce7"        # python – green
LANG_PY_BORDER = "#16a34a"
LANG_BASH = "#dbeafe"      # bash  – blue
LANG_BASH_BORDER = "#2563eb"


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
            "font.family": ["Microsoft YaHei", "Noto Sans SC", "SimHei", "DejaVu Sans"],
            "axes.unicode_minus": False,
            "text.color": COLOR_TEXT,
            "axes.labelcolor": COLOR_TEXT,
            "axes.edgecolor": COLOR_TEXT,
            "xtick.color": COLOR_TEXT,
            "ytick.color": COLOR_TEXT,
            "axes.titlesize": 11.5,
            "axes.titleweight": "bold",
            "axes.titlepad": 6,
            "font.size": 8.5,
        }
    )


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------
def _box(ax, x, y, w, h, text, fc=COLOR_NODE, ec=COLOR_NODE_BORDER, fontsize=8.0,
         fontweight="normal", text_color=None, style="round,pad=0.02", align="center"):
    """Draw a labelled box centred at (x, y) and return its bounding artist."""
    box = mpatches.FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle=style, fc=fc, ec=ec, lw=1.1, zorder=3,
    )
    ax.add_patch(box)
    ha = {"center": "center", "left": "left", "right": "right"}[align]
    ax.text(
        x, y, text, ha=ha, va="center",
        fontsize=fontsize, fontweight=fontweight,
        color=text_color or COLOR_TEXT, zorder=4,
    )
    return box


def _code(ax, x, y, w, h, text, fontsize=7.6):
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
        family=["DejaVu Sans Mono", "Microsoft YaHei", "Noto Sans SC", "SimHei", "DejaVu Sans"],
    )
    return box


def _arrow(ax, x1, y1, x2, y2, label=None, color=COLOR_ARROW, lw=1.4,
           label_offset=(0.0, 0.0), label_color=None):
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
            mx, my, label, ha="center", va="center", fontsize=7.0,
            color=label_color or color, fontweight="bold", zorder=5,
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85),
        )


def _turn_off(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)


# ---------------------------------------------------------------------------
# Panel 1: workflow classification (10 run_* scripts in 5 functional families)
# ---------------------------------------------------------------------------
def draw_panel_classification(ax):
    """Group the 10 run_* scripts into 5 functional classes, each with its y_dir."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_title("面板1 · 10 个 run_* 脚本按功能分 5 类  (体积/应变 · 收敛/基准 · 表面/界面 · 电子结构 · 过渡态)",
                 loc="left", fontsize=9.5)

    # 5 category headers across the top, each spanning a column of workflow cards.
    categories = [
        (1.05, "体积 / 应变类", CAT_VOLUME, CAT_VOLUME_BORDER, [
            ("pei_vasp_run_eos.py",        "y_eos",            "V/V0 0.94..1.06  (isif 4/2)"),
            ("pei_vasp_run_stretch.py",    "y_stretch",        "strain ±0.4%  (xyz/xy/z)"),
            ("pei_vasp_run_cij_energy",    "y_cij_energy",     "2nd-order Cij  (energy-strain)"),
            ("pei_vasp_run_hoec_energy",   "y_hoec_energy",    "2nd/3rd/4th  (Wang-Li)"),
        ]),
        (3.05, "收敛 / 基准类", CAT_CONVERGENCE, CAT_CONVERGENCE_BORDER, [
            ("pei_vasp_run_convergence",   "y_convergence",    "ENCUT / KPOINTS"),
            ("pei_vasp_run_kpar_ncore",    "y_kpar_ncore",     "KPAR×NCORE  (21 对)"),
        ]),
        (5.05, "表面 / 界面类", CAT_SURFACE, CAT_SURFACE_BORDER, [
            ("pei_vasp_run_surface_energy", "y_surface_energy", "bulk vs slab  (γ)"),
            ("pei_vasp_run_cohesive",       "y_cohesive",       "scale 0.60..4.00  (42 pts)"),
        ]),
        (7.05, "电子结构类", CAT_ELECTRON, CAT_ELECTRON_BORDER, [
            ("pei_vasp_run_dos_band",      "y_dos_band",       "DOS + band  (KPATH.in)"),
        ]),
        (9.05, "过渡态类", CAT_NEB, CAT_NEB_BORDER, [
            ("pei_vasp_run_neb",           "y_neb",            "ini→N img→fin  (需 ini/fin CONTCAR)"),
        ]),
    ]

    # Header band at top
    header_y = 9.15
    header_h = 0.75
    for x, title, fc, ec, _ in categories:
        _box(ax, x, header_y, 1.78, header_h, title,
             fc=fc, ec=ec, fontsize=8.4, fontweight="bold")

    # Workflow cards under each header
    card_h = 1.18
    card_gap = 0.16
    for x, _title, fc, ec, items in categories:
        y = header_y - header_h / 2 - card_h / 2 - card_gap
        for script_name, y_dir, one_liner in items:
            card_text = "%s\n%s\n%s" % (script_name, y_dir, one_liner)
            _box(ax, x, y, 1.78, card_h, card_text,
                 fc=fc, ec=ec, fontsize=6.7,
                 text_color=COLOR_TEXT)
            y -= card_h + card_gap

    # Shared entry-point band at the bottom
    _box(ax, 5, 0.95, 9.6, 1.0,
         "统一入口: y_full_relax/  (CONTCAR / INCAR / KPOINTS / POTCAR / Y_CONSTR_LATT / sub.*)\n"
         "所有 run_* 脚本在此级运行, 仅准备 y_<workflow>/y_dir/<case>/ 目录 (不 sbatch); "
         "由 pei_vasp_univ_sbatch_jobs 提交, pei_vasp_univ_post 取数据, pei_vasp_plot_all 出图",
         fc=COLOR_DECISION, ec=COLOR_DECISION_BORDER, fontsize=7.2)

    # connectors: each category header -> entry point (thin grey)
    for x, _t, _fc, _ec, _items in categories:
        _arrow(ax, x, header_y - header_h / 2, x, 1.46, lw=0.5, color="#cbd5e0")

    _turn_off(ax)


# ---------------------------------------------------------------------------
# Panel 2: shared workflow lifecycle (directory tree + pipeline)
# ---------------------------------------------------------------------------
def draw_panel_lifecycle(ax):
    """Show y_full_relax -> run_* -> y_dir/<case> -> sbatch -> post -> plot."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_title("面板2 · 通用 workflow 生命周期  (y_full_relax → run_* → sbatch → post → plot)",
                 loc="left", fontsize=9.5)

    # Stage 1: entry point
    _box(ax, 5, 9.2, 7.6, 0.85,
         "y_full_relax/   ← 已完成结构弛豫 (CONTCAR / INCAR / KPOINTS / POTCAR / sub.*)",
         fc=COLOR_DECISION, fontsize=7.4)

    # Stage 2: run_* script (prepares dirs only)
    _box(ax, 5, 7.95, 7.6, 0.95,
         "pei_vasp_run_<workflow>   (10 个之一, 见面板1)\n"
         "在 y_full_relax 级运行 → 建 y_<workflow>/y_dir/<case>/  (只准备目录, 不 sbatch)",
         fc=COLOR_RUN, ec="#5f9e6f", fontsize=7.2)
    _arrow(ax, 5, 8.78, 5, 8.43, lw=1.3)

    # Stage 3: directory tree (code block on the left)
    _code(ax, 2.55, 6.05, 4.5, 2.0,
          "y_full_relax/\n"
          "├─ y_<workflow>/\n"
          "│  └─ y_dir/\n"
          "│     ├─ <case_1>/   ← POSCAR INCAR KPOINTS\n"
          "│     │                  POTCAR Y_CONSTR_LATT\n"
          "│     │                  sub.slurm\n"
          "│     ├─ <case_2>/\n"
          "│     └─ .../<case_N>/\n"
          "└─ (源 CONTCAR / INCAR / KPOINTS)",
          fontsize=6.6)
    _arrow(ax, 2.55, 7.48, 2.55, 7.06, lw=1.1)

    # Stage 4: sbatch submission (right of tree)
    _box(ax, 7.45, 6.55, 4.5, 1.4,
         "pei_vasp_univ_sbatch_jobs\n"
         "find ./ -type d -name y_dir\n"
         "→ 逐目录 srun vasp_std\n"
         "(exit 0/1/10: 完成/失败/已跳过)",
         fc="#dbeafe", ec="#5a7fb0", fontsize=7.0)
    _arrow(ax, 4.8, 6.55, 5.2, 6.55, lw=1.1)

    # Stage 5: post-processing
    _box(ax, 7.45, 4.55, 4.5, 1.4,
         "pei_vasp_univ_post\n"
         "扫 y_dir → OUTCAR 取\n"
         "energy / stress / volume / time\n"
         "→ y_post_data.txt / y_post_time.txt",
         fc="#ede9fe", ec="#8b5cf6", fontsize=7.0)
    _arrow(ax, 7.45, 5.85, 7.45, 5.25, lw=1.1)

    # Stage 6: plot dispatcher
    _box(ax, 5, 2.75, 9.4, 1.15,
         "pei_vasp_plot_all <workflow>   (见面板3 调度器映射表)\n"
         "find y_<workflow> → 进每个目录调用对应 pei_vasp_plot_*.py → 出图",
         fc="#f5e6d3", ec="#b9802f", fontsize=7.2)
    _arrow(ax, 7.45, 3.85, 7.0, 3.33, lw=1.0, color="#7a8aa0")
    _arrow(ax, 2.55, 5.05, 3.5, 3.33, lw=1.0, color="#7a8aa0")

    # Shared safety policy footer
    _box(ax, 5, 1.1, 9.4, 1.05,
         "共享安全策略 confirm_prepare_dir():  已有 y_* 目录时交互式确认;  "
         "非 tty / 空白回答 → abort (不删除);  -deleteold / -force 跳过提示;  "
         "全部经 y_full_relax_temp 中转, 完成后清理",
         fc="#fff7ed", ec="#d97706", fontsize=6.9)

    _turn_off(ax)


# ---------------------------------------------------------------------------
# Panel 3: pei_vasp_plot_all dispatcher mapping table
# ---------------------------------------------------------------------------
def draw_panel_plotall_table(ax):
    """Tabulate plot_all option -> y_* dir -> plot script."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_title("面板3 · pei_vasp_plot_all 调度器映射表  (option → y_* 目录 → plot 脚本)",
                 loc="left", fontsize=9.5)

    rows = [
        ("plot_all 选项",        "find 的目录",         "调用的 plot 脚本",                True),
        ("-convergence",         "y_convergence",       "pei_vasp_plot_convergence.py",    False),
        ("-stretch",             "y_stretch",           "pei_vasp_plot_stretch.py",        False),
        ("-neb",                 "y_neb",               "pei_vasp_plot_neb.py",            False),
        ("-hoec_energy",         "y_hoec_energy",       "pei_vasp_plot_hoec_energy.py",    False),
        ("-kpar_nocre (-ncore)", "y_kpar_ncore",        "pei_vasp_plot_kpar_ncore",        False),
        ("-E_in_1_2_bulk",       "y_E_in_1_2_bulk",     "pei_vasp_plot_E_in_1_2_bulk.py*", False),
    ]

    cx = [1.85, 4.55, 8.05]
    cw = [3.0, 2.6, 3.6]
    y_top = 8.7
    row_h = 1.05
    y = y_top
    for opt, wdir, scr, is_header in rows:
        if is_header:
            for x, w, txt in zip(cx, cw, (opt, wdir, scr)):
                _box(ax, x, y, w, row_h, txt, fc="#2c3e50",
                     text_color="white", fontsize=7.4, fontweight="bold")
        else:
            fc = "#f7f9fb"
            for x, w, txt, col in zip(cx, cw, (opt, wdir, scr),
                                      (COLOR_DATA1, COLOR_TEXT, COLOR_TEXT)):
                _box(ax, x, y, w, row_h, txt, fc=fc, ec=COLOR_NODE_BORDER,
                     fontsize=6.9,
                     fontweight="bold" if x == cx[0] else "normal",
                     text_color=col)
            _arrow(ax, cx[0] + cw[0] / 2, y, cx[1] - cw[1] / 2, y, lw=0.8, color="#7a8aa0")
            _arrow(ax, cx[1] + cw[1] / 2, y, cx[2] - cw[2] / 2, y, lw=0.8, color="#7a8aa0")
        y -= row_h + 0.08

    _box(ax, 5, 0.95, 9.6, 1.5,
         "*  -E_in_1_2_bulk 的 plot 脚本位于 vasp_utils/vasp_workflow_others/ (非 workflow_bulk 子包);\n"
         "   其余 5 个 plot 脚本均在 vasp_utils/vasp_workflow_bulk/ 下。\n"
         "   -kpar_nocre 是公开拼写, -kpar_ncore 作为别名也被接受 (见 plot_all case 块)。\n"
         "   每个 plot 脚本在其 y_<workflow> 目录内运行, 会重新跑 pei_vasp_univ_post; 仅读 OUTCAR, 不提交作业。",
         fc=COLOR_CODE_BG, fontsize=6.8)

    _turn_off(ax)


# ---------------------------------------------------------------------------
# Panel 4: NEB tool chain (8 scripts, perl/python/bash tagged)
# ---------------------------------------------------------------------------
def draw_panel_neb_chain(ax):
    """NEB pipeline: ini/fin CONTCAR -> run_neb -> 00..0N -> nebbarrier/nebef/nebresults -> nebmovie -> plot -> select."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_title("面板4 · NEB 工具链流程  (neb_utils 8 脚本: perl / python / bash)",
                 loc="left", fontsize=9.5)

    # left column: run_neb setup + barrier curve; right column: neb_utils 8-script chain
    _code(ax, 2.5, 9.15, 4.5, 1.1,
          "$ pei_vasp_run_neb 6   (bash)\n"
          "  源: y_full_relax_neb/{ini,fin}/CONTCAR\n"
          "  推荐 N_IMAGE = DISTANCE / 0.8",
          fontsize=6.8)
    _box(ax, 2.5, 7.75, 4.5, 0.75,
         "nebmake.pl 生成 N images → nebavoid.pl ≥1Å → y_neb/y_dir/00..0N/",
         fc=COLOR_CODE_BG, fontsize=6.6)
    _arrow(ax, 2.5, 8.6, 2.5, 8.13, lw=1.0)

    # barrier curve inset
    axL = ax.inset_axes([0.03, 0.04, 0.45, 0.50], transform=ax.transAxes)
    n_img = 6
    n_pts = n_img + 2
    idx = np.arange(n_pts)
    x_frac = idx / (n_pts - 1)
    E = 0.80 * np.sin(np.pi * x_frac) ** 2
    xs = np.linspace(0, 1, 200)
    Es = 0.80 * np.sin(np.pi * xs) ** 2
    axL.plot(xs, Es, color=COLOR_FIT, lw=1.4, zorder=2)
    axL.scatter(x_frac, E, color=COLOR_DATA1, s=34, zorder=4,
                edgecolor="white", lw=0.6)
    isad = int(np.argmax(E))
    axL.annotate("鞍点 img %d\nE‡ ≈ %.2f eV" % (isad, E[isad]),
                 xy=(x_frac[isad], E[isad]),
                 xytext=(x_frac[isad] - 0.28, E[isad] + 0.06),
                 fontsize=6.6, color=COLOR_HIGHLIGHT,
                 arrowprops=dict(arrowstyle="->", color=COLOR_HIGHLIGHT, lw=1.0))
    axL.set_xlabel("反应坐标 (ini → 6 img → fin)", fontsize=6.8)
    axL.set_ylabel("ΔE (eV)", fontsize=6.8)
    axL.set_title("NEB 能量曲线 (合成, 8 帧)", fontsize=7.8)
    axL.set_xticks(np.linspace(0, 1, n_pts))
    axL.set_xticklabels(["ini"] + ["%02d" % i for i in range(1, n_img + 1)] + ["fin"],
                        fontsize=5.8, rotation=30)
    axL.tick_params(labelsize=6.5)
    axL.grid(True, lw=0.4, alpha=0.4)

    # right: neb_utils 8-script chain (vertical, language-tagged)
    axR = ax.inset_axes([0.52, 0.03, 0.47, 0.92], transform=ax.transAxes)
    axR.set_xlim(0, 10)
    axR.set_ylim(0, 11)
    axR.set_title("neb_utils 8 脚本分工 (perl / python / bash)", fontsize=8.0, loc="left")
    _turn_off(axR)

    stages = [
        (10.3, "OUTCAR (每帧 00..0N/)",                                    "#dbeafe", "#2563eb", "input"),
        (9.05, "pei_vasp_univ_neb_nebbarrier\n→ neb.dat",                   LANG_PERL,  LANG_PERL_BORDER,  "perl"),
        (7.80, "pei_vasp_univ_neb_nebbarrier.py\n→ p_neb_py.dat",           LANG_PY,    LANG_PY_BORDER,    "python"),
        (6.55, "pei_vasp_univ_neb_nebbarrier_full.py\n→ p_neb_full.dat",    LANG_PY,    LANG_PY_BORDER,    "python"),
        (5.30, "pei_vasp_univ_neb_nebef\n→ force/E 汇总",                   LANG_PERL,  LANG_PERL_BORDER,  "perl"),
        (4.05, "pei_vasp_univ_neb_nebresults\n→ 解压 + s.dat 汇总",         LANG_PERL,  LANG_PERL_BORDER,  "perl"),
        (2.80, "pei_vasp_univ_neb_nebmovie\n→ movie_POSCAR.xyz",            LANG_BASH,  LANG_BASH_BORDER,  "bash"),
        (1.55, "pei_vasp_univ_neb_plot.py  (-start/-end)\n→ p_post_neb_full_selected.pdf", LANG_PY, LANG_PY_BORDER, "python"),
    ]
    box_w = 8.6
    box_h = 0.95
    for y, text, fc, ec, lang in stages:
        _box(axR, 5, y, box_w, box_h, text, fc=fc, ec=ec, fontsize=6.7)
        axR.text(8.9, y, lang, ha="right", va="center", fontsize=6.0,
                 fontweight="bold", color=ec, zorder=5,
                 bbox=dict(boxstyle="round,pad=0.12", fc="white", ec=ec, alpha=0.9))
    for i in range(len(stages) - 1):
        _arrow(axR, 5, stages[i][0] - box_h / 2, 5, stages[i + 1][0] + box_h / 2, lw=1.0)

    # neb_select branch (from nebbarrier_full)
    _box(axR, 1.6, 0.35, 4.2, 0.6,
         "pei_vasp_univ_neb_select.py -frame N\n(python, 从 p_neb_full.dat 提取指定帧)",
         fc=LANG_PY, ec=LANG_PY_BORDER, fontsize=5.9)
    _arrow(axR, 3.5, 0.65, 4.5, 1.08, lw=0.8, color="#7a8aa0")

    axR.text(0.25, 5.6, "数据流:\nOUTCAR\n→ neb.dat\n→ p_neb_full.dat\n→ PDF",
             fontsize=6.0, color=COLOR_SKIP, va="center",
             family=["Microsoft YaHei", "Noto Sans SC", "SimHei", "DejaVu Sans"])

    _turn_off(ax)


# ---------------------------------------------------------------------------
# Panel 5: strain-type comparison (stretch -type table + EOS vs stretch)
# ---------------------------------------------------------------------------
def draw_panel_strain_compare(ax):
    """Compare stretch -type (xyz/xy/z, keepvolume) and EOS V/V0 vs stretch strain."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_title("面板5 · pei_vasp_run_stretch.py 应变类型对比  (-type xyz/xy/z, -keepvolume)  &  EOS V/V0 vs stretch strain",
                 loc="left", fontsize=9.2)

    # left: stretch -type comparison table
    hx = [1.75, 4.75]
    hw = [2.6, 3.2]
    y = 8.7
    rh = 0.75
    for x, w, txt in zip(hx, hw, ("-type", "行为 (cell 缩放方式)")):
        _box(ax, x, y, w, rh, txt, fc="#2c3e50", text_color="white",
             fontsize=7.2, fontweight="bold")
    y -= rh + 0.08
    stretch_rows = [
        ("xyz",        "三轴等应变 (各向同性)\ncell 三方向同乘 (1+ε)"),
        ("xy",         "双轴面内 (a,b 拉, c 不变)\n或 -keepvolume: c 反向缩放"),
        ("z",          "单轴 z 方向 (c 拉, a,b 不变)\n或 -keepvolume: a,b 反向缩放"),
        ("-keepvolume","体积守恒: 未拉伸方向\n按 (1+ε)^(-1/n) 反向补偿"),
    ]
    for i, (typ, behav) in enumerate(stretch_rows):
        fc = "#f7f9fb" if i % 2 == 0 else "white"
        _box(ax, hx[0], y, hw[0], 0.95, typ, fc=fc, ec=COLOR_NODE_BORDER,
             fontsize=7.0, fontweight="bold", text_color=COLOR_DATA1)
        _box(ax, hx[1], y, hw[1], 0.95, behav, fc=fc, ec=COLOR_NODE_BORDER,
             fontsize=6.5)
        y -= 0.95 + 0.08

    _box(ax, 3.25, 3.55, 6.4, 0.85,
         "完整 choices: xyz / xy / xz / yz / x / y / z;  默认 strain ±0.4% step 0.05% (17 pts);\n"
         "-strains=... 显式列表 (负值需用 '=' 避免被读成 flag)",
         fc=COLOR_CODE_BG, fontsize=6.4)

    # right: EOS V/V0 vs stretch strain visual contrast
    axR = ax.inset_axes([0.56, 0.30, 0.42, 0.50], transform=ax.transAxes)
    ratios = np.linspace(0.94, 1.06, 200)
    E_eos = 0.5 * 50.0 * (ratios - 1.0) ** 2
    axR.plot(ratios, E_eos, color=COLOR_DATA1, lw=1.6, label="EOS  V/V0 (体积比)")
    strains = np.linspace(-0.004, 0.004, 200)
    E_st = 0.5 * 1.6e4 * strains ** 2
    strain_mapped = 1.0 + strains
    axR.plot(strain_mapped, E_st, color=COLOR_DATA2, lw=1.6, ls="--",
             label="stretch  ε (工程应变)")
    axR.axvline(1.0, color=COLOR_SKIP, lw=0.7, ls=":", alpha=0.7)
    axR.set_xlabel("EOS: V/V0   |   stretch: 1+ε  (同轴对照)", fontsize=6.8)
    axR.set_ylabel("ΔE (eV, 合成)", fontsize=6.8)
    axR.set_title("EOS 体积比 vs stretch 工程应变", fontsize=7.8)
    axR.tick_params(labelsize=6.5)
    axR.grid(True, lw=0.4, alpha=0.4)
    axR.legend(fontsize=6.0, loc="upper center", framealpha=0.9)

    # bottom: EOS vs stretch distinction note
    _box(ax, 5, 1.05, 9.6, 1.35,
         "区别:  pei_vasp_run_eos.py 用体积比 V/V0, cell 按 ratio^(1/3) 各向同性缩放 (三轴同时变), "
         "拟合 Birch-Murnaghan → E0/B0/B0';\n"
         "        pei_vasp_run_stretch.py 用工程应变 ε, 只拉指定方向 (xyz/xy/z/...), 未拉方向不变 "
         "(或 -keepvolume 反向补偿保体积), 二次拟合 E(ε)=E0+½Bε² → 平衡晶格常数 a0。\n"
         "两脚本都在 y_full_relax 级运行, 仅准备目录 (不 sbatch)。",
         fc=COLOR_PANEL_BG, ec="#d1d5db", fontsize=6.7)

    _turn_off(ax)


# ---------------------------------------------------------------------------
# Figure assembly
# ---------------------------------------------------------------------------
def render_overview(path_image: Path) -> None:
    """Build the 5-panel overview figure and save it to path_image."""
    _set_rcparams()
    fig = plt.figure(figsize=(16, 12), dpi=150)
    gs = fig.add_gridspec(
        3, 2, height_ratios=[1.0, 1.0, 1.0], hspace=0.28, wspace=0.10,
        left=0.03, right=0.985, top=0.945, bottom=0.030,
    )
    ax1 = fig.add_subplot(gs[0, 0])    # panel 1: classification
    ax2 = fig.add_subplot(gs[0, 1])    # panel 2: lifecycle
    ax3 = fig.add_subplot(gs[1, 0])    # panel 3: plot_all table
    ax4 = fig.add_subplot(gs[1, 1])    # panel 4: NEB chain
    ax5 = fig.add_subplot(gs[2, :])    # panel 5: strain-type compare (wide)

    draw_panel_classification(ax1)
    draw_panel_lifecycle(ax2)
    draw_panel_plotall_table(ax3)
    draw_panel_neb_chain(ax4)
    draw_panel_strain_compare(ax5)

    fig.suptitle(
        "vasp_utils / vasp_workflow_bulk + neb_utils 子包功能总览",
        fontsize=16, fontweight="bold", color=COLOR_TEXT, y=0.992,
    )

    fig.savefig(path_image, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # ---- non-blank self-check (mandatory per pjvasp-examples skill) ----
    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError("overview image was not created: " + str(path_image))
    image_rgb = plt.imread(path_image)[..., :3]
    if float(np.mean(np.abs(image_rgb - 1.0))) < 0.002:
        raise AssertionError("overview image is effectively blank: " + str(path_image))


def _static_image_path() -> Path:
    """Resolve docs/source/_static/images/generated/ relative to this script."""
    repo_root = Path(__file__).resolve().parents[2]   # docs/examples/x.py -> repo
    return repo_root / "docs" / "source" / "_static" / "images" / "generated"


def run_example(path_output: Path) -> dict:
    """Resolve output dir, render the PNG, mirror to _static, run assertions."""
    path_output = Path(path_output).expanduser().resolve()
    path_output.mkdir(parents=True, exist_ok=True)
    path_image = path_output / "vasp_workflow_bulk_overview.png"

    render_overview(path_image)

    # Mirror the PNG into the Sphinx static tree (for direct doc references),
    # matching the layout used by docs/scripts/generate_structure_images.py.
    static_dir = _static_image_path()
    try:
        static_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path_image, static_dir / path_image.name)
    except OSError:
        # Non-fatal: the --output copy is the authoritative artifact.
        pass

    size_kb = path_image.stat().st_size / 1024.0
    image_rgb = plt.imread(path_image)[..., :3]
    deviation = float(np.mean(np.abs(image_rgb - 1.0)))

    summary = {
        "image_path": str(path_image),
        "static_path": str(static_dir / path_image.name),
        "size_kb": round(size_kb, 1),
        "pixel_deviation": round(deviation, 4),
    }
    print("wrote: " + str(path_image))
    print("mirror: " + str(static_dir / path_image.name))
    print("size_kb: " + str(round(size_kb, 1)))
    print("pixel_deviation: " + str(round(deviation, 4)))
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render the vasp_utils/vasp_workflow_bulk + neb_utils overview figure (VASP-free).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("docs/_build/example-vasp-workflow"),
        help="Output directory for the rendered PNG (default: docs/_build/example-vasp-workflow)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary = run_example(args.output)
    # deterministic assertions
    assert summary["size_kb"] > 10.0, "PNG smaller than 10 KB: " + str(summary)
    assert summary["pixel_deviation"] > 0.002, "PNG effectively blank: " + str(summary)
    print("OK: vasp_workflow_bulk overview figure verified non-blank.")


if __name__ == "__main__":
    main()
