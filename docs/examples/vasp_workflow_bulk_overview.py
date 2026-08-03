#!/usr/bin/env python3
"""Overview figure for the vasp_utils/vasp_workflow_bulk + neb_utils sub-packages.

This tutorial is VASP-free.  It renders a single multi-panel figure that
summarises the core functionality of every ``pei_vasp_run_*`` workflow script
in ``vasp_utils/vasp_workflow_bulk/`` plus the ``neb_utils`` NEB post-processing
tool chain:

  * Panel 1 — the shared ``y_full_relax`` entry point and the
    ``y_<workflow>/y_dir/<strain_or_param>/`` directory layout (ASCII tree +
    block diagram) for all 10 workflow scripts.
  * Panel 2 — EOS (V/V0 7-point Birch-Murnaghan E-V curve) and stretch
    (17-point strain-energy quadratic fit), with the run-script entry points.
  * Panel 3 — NEB migration-barrier energy curve (ini → 6 images → fin) and
    the 8-script neb_utils tool chain with its data-flow arrows.
  * Panel 4 — a 2×3 sub-grid of data panels: convergence (ENCUT/KPOINTS),
    cij_energy (7-strain Cij fit), cohesive (42 scale-point E-k curve),
    hoec_energy (cubic 11 / hex 20 modes), and kpar_ncore (21-pair timing
    heatmap).
  * Panel 5 — dos_band (DOS + band path), surface_energy (bulk vs slab with
    vacuum), and the ``pei_vasp_plot_all`` dispatcher option→script mapping.

No VASP binary, ``sbatch``, ``nebmake.pl``, ``n2p2`` or LAMMPS is invoked.
Every energy / strain / timing series below is synthetic (generated with
numpy) and is only meant to illustrate what each workflow produces.  The
figure is generated with matplotlib only.
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
COLOR_DATA3 = "#7c3aed"   # purple – tertiary data series
COLOR_FIT = "#6b7280"     # grey   – fitted curve
COLOR_HIGHLIGHT = "#9b2c2c"
COLOR_HEATMAP = "YlOrRd"

# Language tag colours for the NEB tool chain (Panel 3)
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
            "axes.titlesize": 11.0,
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


# ---------------------------------------------------------------------------
# Panel 1: workflow overview + directory tree
# ---------------------------------------------------------------------------
def draw_panel_overview_tree(ax):
    """y_full_relax unified entry → 10 y_* workflows, each y_dir/<param>/ layout."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_title("面板1 · workflow 总览与目录树  (y_full_relax → 10 个 y_* / y_dir/<param>)",
                 loc="left", fontsize=9.5)

    # top: unified entry point
    _box(ax, 5, 9.35, 8.8, 0.8,
         "统一入口: y_full_relax/   (已收敛 CONTCAR / INCAR / KPOINTS / POTCAR / sub.*)",
         fc=COLOR_DECISION, ec=COLOR_DECISION_BORDER, fontsize=7.6, fontweight="bold")

    # 10 y_* workflow chips in two rows (5 per row)
    workflows = [
        ("y_eos",            "pei_vasp_run_eos.py",        "V/V0 0.94..1.06"),
        ("y_stretch",        "pei_vasp_run_stretch.py",    "strain ±0.4%"),
        ("y_neb",            "pei_vasp_run_neb",           "ini→N img→fin"),
        ("y_convergence",    "pei_vasp_run_convergence",   "ENCUT/KPOINTS"),
        ("y_cij_energy",     "pei_vasp_run_cij_energy",    "7 应变 Cij"),
        ("y_cohesive",       "pei_vasp_run_cohesive",      "scale 42 pts"),
        ("y_dos_band",       "pei_vasp_run_dos_band",      "DOS+band"),
        ("y_surface_energy", "pei_vasp_run_surface_energy","bulk vs slab"),
        ("y_hoec_energy",    "pei_vasp_run_hoec_energy",   "2/3/4 阶"),
        ("y_kpar_ncore",     "pei_vasp_run_kpar_ncore",    "KPAR×NCORE 21"),
    ]
    chip_w = 1.74
    chip_h = 0.92
    x0 = 0.95 + chip_w / 2
    gap = 0.02
    row_ys = [8.05, 6.95]
    for i, (ydir, script, one) in enumerate(workflows):
        col = i % 5
        row = i // 5
        cx = x0 + col * (chip_w + gap)
        cy = row_ys[row]
        _box(ax, cx, cy, chip_w, chip_h,
             "%s\n%s\n%s" % (ydir, script, one),
             fc=COLOR_NODE, ec=COLOR_NODE_BORDER, fontsize=6.3)
        # connector to entry
        _arrow(ax, cx, cy + chip_h / 2, 5, 8.95, lw=0.4, color="#cbd5e0")

    # connectors row1->row2 omitted (chips independent of each other)

    # bottom: ASCII directory tree (code block) + run-only note
    _code(ax, 2.7, 4.15, 4.6, 2.7,
          "y_full_relax/\n"
          "├─ y_eos/\n"
          "│  └─ y_dir/\n"
          "│     ├─ 1_0.94/   POSCAR INCAR KPOINTS\n"
          "│     │            POTCAR sub.slurm\n"
          "│     ├─ 2_0.96/\n"
          "│     └─ .../7_1.06/\n"
          "├─ y_stretch/\n"
          "│  └─ y_dir/  (1_-0.004 ... 17_+0.004)\n"
          "├─ y_cij_energy/\n"
          "│  └─ y_dir/  (1_-0.003 ... 7_+0.003)\n"
          "├─ y_cohesive/\n"
          "│  └─ y_dir/  (1_0.60 ... 42_4.00)\n"
          "└─ y_<其余 workflow>/y_dir/<param>/",
          fontsize=6.3)

    # right: pipeline note
    _box(ax, 7.5, 4.6, 4.4, 2.0,
         "每个 run_* 脚本:\n"
         "· 在 y_full_relax 级运行\n"
         "· 仅准备 y_<wf>/y_dir/<case>/\n"
         "· 不 sbatch (cohesive 例外)\n"
         "· confirm_prepare_dir 闸:\n"
         "  非 tty / 空回答 → abort\n"
         "  (不删除已有结果)",
         fc=COLOR_CODE_BG, fontsize=6.8)

    # bottom strip: submission + post + plot
    _box(ax, 2.7, 2.55, 4.6, 0.7,
         "pei_vasp_univ_sbatch_jobs → srun vasp_std\n(exit 0/1/10: 完成/失败/跳过)",
         fc="#dbeafe", ec="#5a7fb0", fontsize=6.6)
    _box(ax, 7.5, 2.55, 4.4, 0.7,
         "pei_vasp_univ_post → y_post_*.txt\n(energy/stress/volume/time)",
         fc="#ede9fe", ec="#8b5cf6", fontsize=6.6)

    _box(ax, 5, 1.35, 9.4, 0.95,
         "pei_vasp_plot_all <option>  →  find y_<wf>  →  调用对应 pei_vasp_plot_*.py 出图\n"
         "(见面板5: -convergence/-stretch/-neb/-hoec_energy/-kpar_ncore/-E_in_1_2_bulk)",
         fc="#f5e6d3", ec="#b9802f", fontsize=7.0)
    _arrow(ax, 2.7, 2.2, 3.5, 1.83, lw=0.8, color="#7a8aa0")
    _arrow(ax, 7.5, 2.2, 6.5, 1.83, lw=0.8, color="#7a8aa0")

    _box(ax, 5, 0.45, 9.4, 0.6,
         "全部 run_* 仅准备目录与输入文件; 由 sbatch_jobs 提交, univ_post 取数据, plot_all 出图",
         fc=COLOR_PANEL_BG, ec="#d1d5db", fontsize=6.6)

    _turn_off(ax)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)


# ---------------------------------------------------------------------------
# Panel 2: EOS + stretch
# ---------------------------------------------------------------------------
def draw_panel_eos_stretch(ax):
    """Left: EOS V/V0 7-point Birch-Murnaghan; Right: stretch 17-point quad fit."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_title("面板2 · EOS (V/V0 7点 Birch-Murnaghan) + stretch (17点 二次拟合)",
                 loc="left", fontsize=9.5)

    # ---- left inset: EOS ----
    axL = ax.inset_axes([0.04, 0.10, 0.42, 0.66], transform=ax.transAxes)
    vv0 = np.array([0.94, 0.96, 0.98, 1.00, 1.02, 1.04, 1.06])
    # synthetic Birch-Murnaghan: E(V) = E0 + 9/16 * B0 * V0 * [(η²-1)^3 * B0' - (η²-1)^2], η=(V/V0)^(1/3)
    E0, B0, B0p, V0 = -10.0, 50.0, 4.0, 1.0
    eta = np.cbrt(vv0)
    xi = eta ** 2 - 1.0
    E_eos = E0 + (9.0 / 16.0) * B0 * V0 * (xi ** 3 * B0p - xi ** 2)
    # dense curve
    vv0_dense = np.linspace(0.93, 1.07, 300)
    eta_d = np.cbrt(vv0_dense)
    xi_d = eta_d ** 2 - 1.0
    E_dense = E0 + (9.0 / 16.0) * B0 * V0 * (xi_d ** 3 * B0p - xi_d ** 2)
    axL.plot(vv0_dense, E_dense, color=COLOR_FIT, lw=1.3, zorder=2, label="Birch-Murnaghan 拟合")
    axL.scatter(vv0, E_eos, color=COLOR_DATA1, s=32, zorder=4, edgecolor="white", lw=0.6,
                label="EOS 7 点 (V/V0)")
    axL.axvline(1.0, color=COLOR_SKIP, lw=0.7, ls=":", alpha=0.7)
    axL.annotate("V0 (V/V0=1.0)\nE0=%.2f eV" % E0, xy=(1.0, E0),
                 xytext=(0.955, E0 - 0.15), fontsize=6.4, color=COLOR_HIGHLIGHT,
                 arrowprops=dict(arrowstyle="->", color=COLOR_HIGHLIGHT, lw=0.9))
    axL.set_xlabel("V/V0  (0.94..1.06, step 0.02, 7 点)", fontsize=7.0)
    axL.set_ylabel("E (eV, 合成)", fontsize=7.0)
    axL.set_title("pei_vasp_run_eos.py  -isif 4", fontsize=8.0, color=COLOR_DATA1)
    axL.tick_params(labelsize=6.6)
    axL.grid(True, lw=0.4, alpha=0.4)
    axL.legend(fontsize=5.8, loc="upper center", framealpha=0.9)

    # ---- right inset: stretch ----
    axR = ax.inset_axes([0.54, 0.10, 0.42, 0.66], transform=ax.transAxes)
    strain = np.linspace(-0.004, 0.004, 17)
    # synthetic quadratic: E = E0 + 0.5 * K * strain^2 + small noise
    K = 1.6e4
    E_st = 0.5 * K * strain ** 2
    # dense fit curve
    strain_d = np.linspace(-0.0045, 0.0045, 300)
    E_fit = 0.5 * K * strain_d ** 2
    axR.plot(strain_d, E_fit, color=COLOR_FIT, lw=1.3, zorder=2, label="二次拟合 E=E0+½Bε²")
    axR.scatter(strain, E_st, color=COLOR_DATA2, s=22, zorder=4, edgecolor="white", lw=0.5,
                label="stretch 17 点 (±0.4%)")
    axR.axvline(0.0, color=COLOR_SKIP, lw=0.7, ls=":", alpha=0.7)
    axR.set_xlabel("strain ε  (-0.004..0.004, step 0.0005, 17 点)", fontsize=7.0)
    axR.set_ylabel("ΔE (eV, 合成)", fontsize=7.0)
    axR.set_title("pei_vasp_run_stretch.py  -type xyz", fontsize=8.0, color=COLOR_DATA2)
    axR.tick_params(labelsize=6.6)
    axR.grid(True, lw=0.4, alpha=0.4)
    axR.legend(fontsize=5.8, loc="upper center", framealpha=0.9)

    # bottom distinction note
    _box(ax, 5, 1.15, 9.4, 1.25,
         "区别:  pei_vasp_run_eos.py 用体积比 V/V0, cell 按 ratio^(1/3) 各向同性缩放, 拟合 Birch-Murnaghan → E0/B0/B0';\n"
         "        pei_vasp_run_stretch.py 用工程应变 ε, 只拉指定方向 (-type xyz/xy/z/...), 二次拟合 E(ε)=E0+½Bε² → 平衡晶格常数 a0。\n"
         "两脚本均在 y_full_relax 级运行, 仅准备 y_eos/y_stretch 目录 (不 sbatch)。",
         fc=COLOR_PANEL_BG, ec="#d1d5db", fontsize=6.8)

    _turn_off(ax)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)


# ---------------------------------------------------------------------------
# Panel 3: NEB (energy curve + 8-script tool chain)
# ---------------------------------------------------------------------------
def draw_panel_neb(ax):
    """NEB: ini→6 img→fin energy curve + neb_utils 8-script chain + data-flow."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_title("面板3 · NEB 迁移能垒 + neb_utils 8 脚本工具链",
                 loc="left", fontsize=9.5)

    # ---- left: run_neb header + energy curve ----
    _code(ax, 2.55, 9.2, 4.5, 0.95,
          "$ pei_vasp_run_neb 6   (bash)\n"
          "  源: y_full_relax_neb/{ini,fin}/CONTCAR\n"
          "  推荐 N_IMAGE = DISTANCE / 0.8",
          fontsize=6.7)

    axL = ax.inset_axes([0.03, 0.34, 0.46, 0.50], transform=ax.transAxes)
    n_img = 6
    n_pts = n_img + 2
    idx = np.arange(n_pts)
    x_frac = idx / (n_pts - 1)
    E = 0.80 * np.sin(np.pi * x_frac) ** 2
    xs = np.linspace(0, 1, 200)
    Es = 0.80 * np.sin(np.pi * xs) ** 2
    axL.plot(xs, Es, color=COLOR_FIT, lw=1.4, zorder=2)
    axL.scatter(x_frac, E, color=COLOR_DATA1, s=30, zorder=4,
                edgecolor="white", lw=0.6)
    isad = int(np.argmax(E))
    axL.annotate("鞍点 img %d\nE‡ ≈ %.2f eV" % (isad, E[isad]),
                 xy=(x_frac[isad], E[isad]),
                 xytext=(x_frac[isad] - 0.30, E[isad] + 0.08),
                 fontsize=6.4, color=COLOR_HIGHLIGHT,
                 arrowprops=dict(arrowstyle="->", color=COLOR_HIGHLIGHT, lw=1.0))
    axL.set_xlabel("反应坐标 (ini → 6 img → fin)", fontsize=6.8)
    axL.set_ylabel("ΔE (eV)", fontsize=6.8)
    axL.set_title("NEB 能量曲线 (合成, 8 帧, 能垒 ~0.8 eV)", fontsize=7.4)
    axL.set_xticks(np.linspace(0, 1, n_pts))
    axL.set_xticklabels(["ini"] + ["%02d" % i for i in range(1, n_img + 1)] + ["fin"],
                        fontsize=5.6, rotation=30)
    axL.tick_params(labelsize=6.4)
    axL.grid(True, lw=0.4, alpha=0.4)

    # data-flow arrow strip under the curve
    _box(ax, 2.55, 1.55, 4.5, 1.0,
         "数据流:  OUTCAR → nebbarrier → neb.dat →\n"
         "neb_plot.py → p_post_neb_full_selected.pdf",
         fc=COLOR_CODE_BG, fontsize=6.6)

    # ---- right: neb_utils 8-script chain (vertical, language-tagged) ----
    axR = ax.inset_axes([0.53, 0.04, 0.45, 0.88], transform=ax.transAxes)
    axR.set_xlim(0, 10)
    axR.set_ylim(0, 11)
    axR.set_title("neb_utils 8 脚本分工 (perl / python / bash)", fontsize=7.8, loc="left")
    axR.set_xticks([])
    axR.set_yticks([])
    for s in axR.spines.values():
        s.set_visible(False)

    stages = [
        (10.4, "OUTCAR (每帧 00..0N/)",                                    "#dbeafe", "#2563eb", "input"),
        (9.15, "pei_vasp_univ_neb_nebbarrier\n→ neb.dat  (perl vtst)",     LANG_PERL,  LANG_PERL_BORDER,  "perl"),
        (7.90, "pei_vasp_univ_neb_nebbarrier.py\n→ p_neb_py.dat",           LANG_PY,    LANG_PY_BORDER,    "python"),
        (6.65, "pei_vasp_univ_neb_nebbarrier_full.py\n→ p_neb_full.dat (每帧)", LANG_PY, LANG_PY_BORDER, "python"),
        (5.40, "pei_vasp_univ_neb_nebef\n→ force/E 汇总",                   LANG_PERL,  LANG_PERL_BORDER,  "perl"),
        (4.15, "pei_vasp_univ_neb_nebresults\n→ 解压 + s.dat 汇总",         LANG_PERL,  LANG_PERL_BORDER,  "perl"),
        (2.90, "pei_vasp_univ_neb_nebmovie\n→ movie_POSCAR.xyz",            LANG_BASH,  LANG_BASH_BORDER,  "bash"),
        (1.55, "pei_vasp_univ_neb_plot.py  (-start/-end)\n→ p_post_neb_full_selected.pdf", LANG_PY, LANG_PY_BORDER, "python"),
    ]
    box_w = 8.8
    box_h = 0.95
    for y, text, fc, ec, lang in stages:
        _box(axR, 5, y, box_w, box_h, text, fc=fc, ec=ec, fontsize=6.5)
        axR.text(9.0, y, lang, ha="right", va="center", fontsize=5.8,
                 fontweight="bold", color=ec, zorder=5,
                 bbox=dict(boxstyle="round,pad=0.12", fc="white", ec=ec, alpha=0.9))
    for i in range(len(stages) - 1):
        _arrow(axR, 5, stages[i][0] - box_h / 2, 5, stages[i + 1][0] + box_h / 2, lw=1.0)

    # neb_select branch (from nebbarrier_full)
    _box(axR, 1.6, 0.35, 4.4, 0.6,
         "pei_vasp_univ_neb_select.py -frame N\n(python, 从 p_neb_full.dat 提取指定帧)",
         fc=LANG_PY, ec=LANG_PY_BORDER, fontsize=5.8)
    _arrow(axR, 3.5, 0.65, 4.5, 1.08, lw=0.8, color="#7a8aa0")

    _turn_off(ax)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)


# ---------------------------------------------------------------------------
# Panel 4: convergence + cij + cohesive + hoec + kpar_ncore (2×3 sub-grid)
# ---------------------------------------------------------------------------
def draw_panel_multi(ax):
    """2×3 sub-grid of data panels: convergence, cij, cohesive, hoec, kpar_ncore."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_title("面板4 · convergence + cij_energy + cohesive + hoec_energy + kpar_ncore",
                 loc="left", fontsize=9.5)
    ax.axis("off")

    # 6 inset cells: 5 data plots + 1 legend cell
    positions = [
        (0.04, 0.52, 0.30, 0.42),   # (a) convergence
        (0.36, 0.52, 0.30, 0.42),   # (b) cij_energy
        (0.68, 0.52, 0.30, 0.42),   # (c) cohesive
        (0.04, 0.04, 0.30, 0.42),   # (d) hoec_energy
        (0.36, 0.04, 0.30, 0.42),   # (e) kpar_ncore heatmap
        (0.68, 0.04, 0.30, 0.42),   # (f) legend / summary
    ]

    # ---- (a) convergence: ENCUT(400,500) / KPOINTS(3x3x3,4x4x4) ----
    ax_a = ax.inset_axes(positions[0], transform=ax.transAxes)
    encut = np.array([300, 350, 400, 450, 500, 550, 600])
    E_400 = -10.0 + 50.0 / (encut - 280.0) ** 1.2 + np.array([0, 0, 0, 0, 0, 0, 0]) * 0.0
    E_500 = -10.02 + 45.0 / (encut - 280.0) ** 1.2
    ax_a.plot(encut, E_400, "-o", color=COLOR_DATA1, lw=1.2, ms=3.5, label="ENCUT scan")
    ax_a.axhline(E_400[-1], color=COLOR_DATA1, lw=0.7, ls=":", alpha=0.6)
    ax_a.set_xlabel("ENCUT (eV)", fontsize=6.6)
    ax_a.set_ylabel("E (eV, 合成)", fontsize=6.6)
    ax_a.set_title("(a) convergence: ENCUT 400/500", fontsize=7.2)
    # KPOINTS inset text
    ax_a.text(0.97, 0.30, "KPOINTS:\n3×3×3 → 4×4×4\nΔE < 1 meV",
              transform=ax_a.transAxes, fontsize=5.8, ha="right", va="center",
              bbox=dict(boxstyle="round,pad=0.2", fc="#fff7ed", ec="#d97706", alpha=0.9))
    ax_a.tick_params(labelsize=6.0)
    ax_a.grid(True, lw=0.4, alpha=0.4)
    ax_a.legend(fontsize=5.6, loc="lower right")

    # ---- (b) cij_energy: 7 strain points (-0.003..0.003) quad fit, mark Cij ----
    ax_b = ax.inset_axes(positions[1], transform=ax.transAxes)
    strain_cij = np.linspace(-0.003, 0.003, 7)
    C11 = 200.0  # GPa synthetic
    E_cij = 0.5 * C11 * strain_cij ** 2 * 1e3  # scale for visibility
    strain_d = np.linspace(-0.0035, 0.0035, 200)
    E_fit = 0.5 * C11 * strain_d ** 2 * 1e3
    ax_b.plot(strain_d, E_fit, color=COLOR_FIT, lw=1.2, zorder=2)
    ax_b.scatter(strain_cij, E_cij, color=COLOR_DATA2, s=26, zorder=4, edgecolor="white", lw=0.5)
    ax_b.axvline(0, color=COLOR_SKIP, lw=0.6, ls=":", alpha=0.6)
    ax_b.annotate("C11 ≈ %d GPa\n(½ ∂²E/∂ε²)" % C11, xy=(0.0025, 0.5 * C11 * 0.0025 ** 2 * 1e3),
                  xytext=(-0.0015, 0.55), fontsize=6.0, color=COLOR_HIGHLIGHT,
                  arrowprops=dict(arrowstyle="->", color=COLOR_HIGHLIGHT, lw=0.8))
    ax_b.set_xlabel("strain ε  (-0.003..0.003, 7 点)", fontsize=6.6)
    ax_b.set_ylabel("ΔE (meV, 合成)", fontsize=6.6)
    ax_b.set_title("(b) cij_energy: 7 应变二次拟合", fontsize=7.2)
    ax_b.tick_params(labelsize=6.0)
    ax_b.grid(True, lw=0.4, alpha=0.4)

    # ---- (c) cohesive: 42 scale points (0.60..4.00) E-k curve, mark k=4.00 ref ----
    ax_c = ax.inset_axes(positions[2], transform=ax.transAxes)
    scale = np.linspace(0.60, 4.00, 42)
    # synthetic cohesive: Rose-Vinet-like, min near k=1.0, plateau by k=4.0
    k = scale
    E_coh = -3.5 * np.exp(-2.0 * (k - 1.0)) + 0.8 * (k - 1.0) ** 2 - 3.0
    ax_c.plot(k, E_coh, "-", color=COLOR_DATA3, lw=1.3, zorder=2)
    ax_c.scatter(k, E_coh, color=COLOR_DATA3, s=10, zorder=4)
    ax_c.axhline(E_coh[-1], color=COLOR_HIGHLIGHT, lw=0.8, ls="--", alpha=0.8,
                 label="k=4.00 参考平台")
    ax_c.axvline(4.00, color=COLOR_SKIP, lw=0.7, ls=":", alpha=0.6)
    ax_c.annotate("k=4.00\n参考平台\nE=%.2f eV" % E_coh[-1], xy=(4.00, E_coh[-1]),
                  xytext=(3.0, E_coh[-1] - 1.2), fontsize=5.8, color=COLOR_HIGHLIGHT,
                  arrowprops=dict(arrowstyle="->", color=COLOR_HIGHLIGHT, lw=0.8))
    ax_c.set_xlabel("scale k  (0.60..4.00, 42 点)", fontsize=6.6)
    ax_c.set_ylabel("E (eV, 合成)", fontsize=6.6)
    ax_c.set_title("(c) cohesive: E-k 曲线", fontsize=7.2)
    ax_c.tick_params(labelsize=6.0)
    ax_c.grid(True, lw=0.4, alpha=0.4)
    ax_c.legend(fontsize=5.6, loc="lower right")

    # ---- (d) hoec_energy: cubic 11 modes / hex 20 modes示意 ----
    ax_d = ax.inset_axes(positions[3], transform=ax.transAxes)
    modes_cubic = np.arange(1, 12)
    modes_hex = np.arange(1, 21)
    E_cubic = 0.5 * np.array([3.0, 2.5, 2.0, 1.8, 1.5, 1.3, 1.1, 0.9, 0.8, 0.7, 0.6]) * np.abs(
        np.sin(modes_cubic * 0.7)) + 0.3 * modes_cubic ** 2 * 1e-3
    ax_d.bar(modes_cubic - 0.18, E_cubic, width=0.32, color=COLOR_DATA1, label="cubic 11 modes")
    E_hex = 0.5 * np.abs(np.sin(modes_hex * 0.5)) + 0.02 * modes_hex
    ax_d.bar(modes_hex + 0.18, E_hex, width=0.32, color=COLOR_DATA2, alpha=0.7, label="hex 20 modes")
    ax_d.set_xlabel("mode index", fontsize=6.6)
    ax_d.set_ylabel("ΔE (meV, 合成)", fontsize=6.6)
    ax_d.set_title("(d) hoec_energy: cubic 11 / hex 20 modes", fontsize=7.2)
    ax_d.tick_params(labelsize=6.0)
    ax_d.grid(True, lw=0.4, alpha=0.4, axis="y")
    ax_d.legend(fontsize=5.4, loc="upper left")

    # ---- (e) kpar_ncore: 21-pair timing heatmap (KPAR×NCORE) ----
    ax_e = ax.inset_axes(positions[4], transform=ax.transAxes)
    # synthetic: KPAR in [1,2,4,8], NCORE in [1,2,4,8,16] -> 20 cells, +1 = 21 pairs
    kpar_vals = np.array([1, 2, 4, 8])
    ncore_vals = np.array([1, 2, 4, 8, 16])
    # timing: lower is better; optimum around KPAR=4,NCORE=4
    timing = np.zeros((len(ncore_vals), len(kpar_vals)))
    for i, nc in enumerate(ncore_vals):
        for j, kp in enumerate(kpar_vals):
            # synthetic: time decreases with product but penalty for imbalance
            prod = kp * nc
            timing[i, j] = 1000.0 / (1.0 + 0.15 * prod) + 20.0 * abs(np.log2(kp) - 2.0)
    im = ax_e.imshow(timing, cmap=COLOR_HEATMAP, aspect="auto", origin="lower")
    ax_e.set_xticks(np.arange(len(kpar_vals)))
    ax_e.set_xticklabels(["KPAR=%d" % k for k in kpar_vals], fontsize=5.6, rotation=20)
    ax_e.set_yticks(np.arange(len(ncore_vals)))
    ax_e.set_yticklabels(["NC=%d" % n for n in ncore_vals], fontsize=5.6)
    # annotate fastest cell
    iy, ix = np.unravel_index(np.argmin(timing), timing.shape)
    ax_e.add_patch(mpatches.Rectangle((ix - 0.5, iy - 0.5), 1, 1, fill=False,
                                      edgecolor=COLOR_HIGHLIGHT, lw=2.0))
    ax_e.text(ix, iy, "fastest\n%.0f s" % timing[iy, ix], ha="center", va="center",
              fontsize=5.4, color=COLOR_HIGHLIGHT, fontweight="bold")
    ax_e.set_title("(e) kpar_ncore: 21-pair timing 热力图", fontsize=7.2)
    cbar = plt.colorbar(im, ax=ax_e, fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=5.4)
    cbar.set_label("time (s, 合成)", fontsize=6.0)

    # ---- (f) legend / summary cell ----
    ax_f = ax.inset_axes(positions[5], transform=ax.transAxes)
    ax_f.axis("off")
    ax_f.text(0.02, 0.96,
              "(f) 各 run_* 数据维度汇总\n"
              "────────────────────────\n"
              "convergence:\n"
              "  ENCUT {400,500} × KPOINTS\n"
              "  {3×3×3, 4×4×4}\n"
              "────────────────────────\n"
              "cij_energy:\n"
              "  7 应变 (-0.003..0.003)\n"
              "  二次拟合 → Cij\n"
              "────────────────────────\n"
              "cohesive:\n"
              "  42 scale (0.60..4.00)\n"
              "  k=4.00 作参考平台\n"
              "────────────────────────\n"
              "hoec_energy:\n"
              "  cubic 11 modes / hex 20\n"
              "  → SOEC/TOEC/FOEC\n"
              "────────────────────────\n"
              "kpar_ncore:\n"
              "  21 对 (KPAR×NCORE)\n"
              "  timing 热力图 → 最快组合",
              transform=ax_f.transAxes, fontsize=6.2, va="top", ha="left",
              family=["Microsoft YaHei", "Noto Sans SC", "SimHei", "DejaVu Sans"],
              bbox=dict(boxstyle="round,pad=0.3", fc=COLOR_CODE_BG, ec="#d1d5db", lw=0.8))


# ---------------------------------------------------------------------------
# Panel 5: dos_band + surface_energy + plot_all dispatcher
# ---------------------------------------------------------------------------
def draw_panel_dos_surface_plotall(ax):
    """dos_band (DOS+band), surface_energy (bulk vs slab), plot_all dispatcher."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_title("面板5 · dos_band + surface_energy + pei_vasp_plot_all 调度器",
                 loc="left", fontsize=9.5)

    # ---- left: dos_band (DOS + band path) ----
    axL = ax.inset_axes([0.03, 0.10, 0.34, 0.76], transform=ax.transAxes)
    # DOS curve (filled)
    e_range = np.linspace(-6, 6, 300)
    dos = (3.0 * np.exp(-((e_range + 2.0) / 1.2) ** 2)
           + 1.8 * np.exp(-((e_range - 0.5) / 1.8) ** 2)
           + 0.6 * np.exp(-((e_range - 3.0) / 1.0) ** 2))
    axL.fill_betweenx(e_range, 0, dos, color=COLOR_DATA1, alpha=0.35, zorder=2)
    axL.plot(dos, e_range, color=COLOR_DATA1, lw=1.2, zorder=3)
    axL.axhline(0, color=COLOR_SKIP, lw=0.7, ls=":", alpha=0.7)
    axL.set_xlabel("DOS (a.u., 合成)", fontsize=6.6)
    axL.set_ylabel("E - Ef (eV)", fontsize=6.6)
    axL.set_title("dos_band: DOS + band (KPATH.in)", fontsize=7.2)
    axL.tick_params(labelsize=6.0)
    axL.grid(True, lw=0.4, alpha=0.4)
    # band path overlay (second axis)
    axL2 = axL.twinx()
    kpath = np.linspace(0, 1, 100)
    band = 2.0 * np.sin(3 * np.pi * kpath) + 1.0 * np.cos(5 * np.pi * kpath) - 1.0
    axL2.plot(kpath, band, color=COLOR_DATA2, lw=1.0, alpha=0.8)
    axL2.set_yticks([])
    # KPATH labels
    kpoints_lbl = ["Γ", "K", "M", "Γ", "A", "L"]
    for i, lbl in enumerate(kpoints_lbl):
        xpos = i / (len(kpoints_lbl) - 1)
        axL.axvline(xpos * dos.max() * 0.0, color=COLOR_SKIP, lw=0.0)  # noop
    axL2.set_xticks(np.linspace(0, 1, len(kpoints_lbl)))
    axL2.set_xticklabels(kpoints_lbl, fontsize=5.6)
    axL2.tick_params(labelsize=5.8)

    # ---- middle: surface_energy (bulk vs slab, mark vacuum) ----
    axM = ax.inset_axes([0.40, 0.10, 0.30, 0.76], transform=ax.transAxes)
    # bar chart: bulk energy vs slab energy (with vacuum layer)
    # use numeric x positions to avoid StrCategoryConverter issues with "\n"
    bar_x = np.array([0, 1])
    bar_labels = ["bulk\n(无真空)", "slab\n(真空层 ~15Å)"]
    energies = np.array([-8.50, -8.32])
    colors = [COLOR_DATA1, COLOR_DATA2]
    bars = axM.bar(bar_x, energies, color=colors, width=0.55, edgecolor="white", lw=1.0)
    axM.set_xticks(bar_x)
    axM.set_xticklabels(bar_labels, fontsize=6.0)
    axM.set_ylabel("E (eV, 合成)", fontsize=6.6)
    axM.set_title("surface_energy: bulk vs slab", fontsize=7.2)
    axM.tick_params(labelsize=6.0)
    axM.grid(True, lw=0.4, alpha=0.4, axis="y")
    # annotate vacuum
    axM.annotate("真空层\n~15 Å\n(加 z 方向)", xy=(1, energies[1]),
                 xytext=(0.4, -7.6), fontsize=5.8, color=COLOR_HIGHLIGHT,
                 ha="center",
                 arrowprops=dict(arrowstyle="->", color=COLOR_HIGHLIGHT, lw=0.8))
    # surface energy formula
    axM.text(0.5, 0.02, "γ = (E_slab − N·E_bulk) / (2·A)",
             transform=axM.transAxes, ha="center", fontsize=5.8, color=COLOR_TEXT,
             bbox=dict(boxstyle="round,pad=0.2", fc="#fff7ed", ec="#d97706", alpha=0.9))

    # ---- right: plot_all dispatcher table ----
    axR = ax.inset_axes([0.73, 0.10, 0.25, 0.76], transform=ax.transAxes)
    axR.axis("off")
    axR.set_xlim(0, 1)
    axR.set_ylim(0, 1)
    axR.text(0.5, 0.97, "pei_vasp_plot_all 调度器", ha="center", va="top",
             fontsize=7.4, fontweight="bold", color=COLOR_TEXT,
             transform=axR.transAxes)
    rows = [
        ("选项", "y_* 目录", "plot 脚本"),
        ("-convergence", "y_convergence", "plot_convergence.py"),
        ("-stretch", "y_stretch", "plot_stretch.py"),
        ("-neb", "y_neb", "plot_neb.py"),
        ("-hoec_energy", "y_hoec_energy", "plot_hoec_energy.py"),
        ("-kpar_ncore", "y_kpar_ncore", "plot_kpar_ncore"),
        ("-E_in_1_2_bulk", "y_E_in_1_2_bulk", "plot_E_in_1_2_bulk.py*"),
    ]
    y_top = 0.88
    rh = 0.115
    y = y_top
    for ri, (opt, wdir, scr) in enumerate(rows):
        is_header = (ri == 0)
        fc = "#2c3e50" if is_header else ("#f7f9fb" if ri % 2 == 0 else "white")
        tc = "white" if is_header else COLOR_TEXT
        fw = "bold" if is_header or (ri > 0) else "normal"
        axR.text(0.08, y, opt, ha="left", va="center", fontsize=5.8,
                 color=("white" if is_header else COLOR_DATA1), fontweight="bold",
                 transform=axR.transAxes)
        axR.text(0.46, y, wdir, ha="left", va="center", fontsize=5.6,
                 color=tc, fontweight=fw, transform=axR.transAxes)
        axR.text(0.74, y, scr, ha="left", va="center", fontsize=5.4,
                 color=tc, fontweight=fw, transform=axR.transAxes,
                 family=["DejaVu Sans Mono", "Microsoft YaHei", "DejaVu Sans"])
        # row background
        if not is_header:
            axR.add_patch(mpatches.Rectangle((0.02, y - rh / 2), 0.96, rh,
                                             fc=fc, ec="none", zorder=0))
            # redraw text above bg
            axR.text(0.08, y, opt, ha="left", va="center", fontsize=5.8,
                     color=COLOR_DATA1, fontweight="bold", zorder=2,
                     transform=axR.transAxes)
            axR.text(0.46, y, wdir, ha="left", va="center", fontsize=5.6,
                     color=COLOR_TEXT, zorder=2, transform=axR.transAxes)
            axR.text(0.74, y, scr, ha="left", va="center", fontsize=5.4,
                     color=COLOR_TEXT, zorder=2, transform=axR.transAxes,
                     family=["DejaVu Sans Mono", "Microsoft YaHei", "DejaVu Sans"])
        else:
            axR.add_patch(mpatches.Rectangle((0.02, y - rh / 2), 0.96, rh,
                                             fc=fc, ec="none", zorder=0))
            axR.text(0.08, y, opt, ha="left", va="center", fontsize=5.8,
                     color="white", fontweight="bold", zorder=2, transform=axR.transAxes)
            axR.text(0.46, y, wdir, ha="left", va="center", fontsize=5.6,
                     color="white", fontweight="bold", zorder=2, transform=axR.transAxes)
            axR.text(0.74, y, scr, ha="left", va="center", fontsize=5.4,
                     color="white", fontweight="bold", zorder=2, transform=axR.transAxes,
                     family=["DejaVu Sans Mono", "Microsoft YaHei", "DejaVu Sans"])
        y -= rh + 0.005
    # footer note
    axR.text(0.5, 0.02,
             "* E_in_1_2_bulk 的 plot 脚本\n  位于 vasp_workflow_others/\n  其余 5 个在 workflow_bulk/",
             ha="center", va="bottom", fontsize=5.2, color=COLOR_SKIP,
             transform=axR.transAxes)

    _turn_off(ax)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)


# ---------------------------------------------------------------------------
# Figure assembly
# ---------------------------------------------------------------------------
def render_overview(path_image: Path) -> None:
    """Build the 5-panel overview figure and save it to path_image."""
    _set_rcparams()
    fig = plt.figure(figsize=(16, 14), dpi=150)
    # 3 rows x 2 cols; panel 1 (top-left) + panel 2 (top-right) + panel 3
    # (mid-left) + panel 4 (mid-right, wide-ish) + panel 5 (bottom, full width)
    gs = fig.add_gridspec(
        3, 2, height_ratios=[1.0, 1.05, 1.0], hspace=0.22, wspace=0.08,
        left=0.03, right=0.985, top=0.962, bottom=0.022,
    )
    ax1 = fig.add_subplot(gs[0, 0])    # panel 1: overview + tree
    ax2 = fig.add_subplot(gs[0, 1])    # panel 2: EOS + stretch
    ax3 = fig.add_subplot(gs[1, 0])    # panel 3: NEB
    ax4 = fig.add_subplot(gs[1, 1])    # panel 4: 2×3 multi data grid
    ax5 = fig.add_subplot(gs[2, :])    # panel 5: dos_band + surface + plot_all

    draw_panel_overview_tree(ax1)
    draw_panel_eos_stretch(ax2)
    draw_panel_neb(ax3)
    draw_panel_multi(ax4)
    draw_panel_dos_surface_plotall(ax5)

    fig.suptitle(
        "vasp_utils / vasp_workflow_bulk + neb_utils 子包功能总览",
        fontsize=16, fontweight="bold", color=COLOR_TEXT, y=0.988,
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
        default=Path("docs/_build/example-vasp-workflow-bulk"),
        help="Output directory for the rendered PNG (default: docs/_build/example-vasp-workflow-bulk)",
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
