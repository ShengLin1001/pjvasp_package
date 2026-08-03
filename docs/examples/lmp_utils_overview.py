#!/usr/bin/env python3
"""Overview figure for the lmp_utils sub-package.

This tutorial is LAMMPS-free. It demonstrates the core functionality of the
lmp_utils sub-package (template/ + post/ + lmp_universal/) by rendering a single
multi-panel figure that summarises:

  * the pei_lmp_run_properties runner flow: stretch → Cij → gsfe with post
    processing and a final Summary,
  * the template-structure comparison of the three calculations (stretch /
    Cij_energy / gsfe) covering relaxation mode, deformation method, and output,
  * how each ``.in`` file ``include`` s its ``.mod`` module files,
  * the GSFE slip-system registry (gsfe_type → bp1/bp2 shear offsets), and
  * how pei_lmp_run_properties uses ``sed`` to substitute placeholders inside
    the template files at runtime.

No LAMMPS binary, no ``srun``, no real calculation is invoked. The figure is
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
COLOR_STRETCH = "#dbeafe"
COLOR_STRETCH_BORDER = "#2563eb"
COLOR_CIJ = "#dcfce7"
COLOR_CIJ_BORDER = "#16a34a"
COLOR_GSFE = "#ede9fe"
COLOR_GSFE_BORDER = "#7c3aed"
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
            # CJK + Latin per-glyph fallback: Chinese glyphs resolve to
            # Microsoft YaHei, while ✓ / – / Latin fall back to DejaVu Sans.
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
# Panel 1: pei_lmp_run_properties runner flow
# ---------------------------------------------------------------------------
def draw_panel_runner_flow(ax):
    """Flowchart of pei_lmp_run_properties: stretch → Cij → gsfe → Summary."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 11)
    ax.set_title("面板1 · pei_lmp_run_properties 总流程 (stretch → Cij → gsfe)", loc="left")

    # ---- stretch stage ----
    _box(ax, 2.0, 10.3, 3.0, 0.55, "[1] stretch 阶段\nfor lat in fcc bcc hcp",
         fc=COLOR_STRETCH, ec=COLOR_STRETCH_BORDER, fontsize=8, fontweight="bold")
    # three lattice boxes
    for i, (lat, aa) in enumerate([("fcc", "4.08"), ("bcc", "3.20"), ("hcp", "2.88")]):
        lx = 0.7 + i * 1.3
        _box(ax, lx, 9.3, 1.15, 0.7,
             f"{lat}\nini_aa={aa}\nlatnum", fc=COLOR_NODE, fontsize=6.6)
        _arrow(ax, 2.0, 10.0, lx, 9.68, lw=0.9, color="#7a8aa0")
    # lmp run + dir label
    _box(ax, 2.0, 8.15, 3.4, 0.6,
         "srun lmp -in stretch_template.in\n→ y_stretch/<lat>/dump/",
         fc=COLOR_RUN, fontsize=6.8)
    _arrow(ax, 2.0, 8.95, 2.0, 8.47, lw=1.0, color="#7a8aa0")
    # post
    _box(ax, 2.0, 7.15, 3.0, 0.6,
         "post/stretch.py\npost_lammps_stretch\n→ p_post_stretch.txt",
         fc=COLOR_CODE_BG, fontsize=6.8)
    _arrow(ax, 2.0, 7.85, 2.0, 7.47, lw=1.0, color="#7a8aa0")

    # ---- Cij stage ----
    _box(ax, 5.0, 6.15, 3.0, 0.55, "[2] Cij_energy 阶段\nfor lat in fcc bcc hcp",
         fc=COLOR_CIJ, ec=COLOR_CIJ_BORDER, fontsize=8, fontweight="bold")
    _box(ax, 5.0, 5.15, 3.4, 0.7,
         "读 a0/c0 ← p_post_stretch.txt\nsrun lmp -in Cij_energy_template.in\n→ y_Cij_energy/<lat>/dump/",
         fc=COLOR_RUN, fontsize=6.6)
    _arrow(ax, 5.0, 5.87, 5.0, 5.52, lw=1.0, color="#7a8aa0")
    _box(ax, 5.0, 4.15, 3.0, 0.6,
         "post/Cij_energy.py\npost_lammps_Cij_energy\n→ y_post_cij_energy.txt",
         fc=COLOR_CODE_BG, fontsize=6.8)
    _arrow(ax, 5.0, 4.80, 5.0, 4.47, lw=1.0, color="#7a8aa0")

    # ---- gsfe stage ----
    _box(ax, 8.0, 3.15, 3.0, 0.55, "[3] gsfe 阶段\nfor lat in fcc hcp",
         fc=COLOR_GSFE, ec=COLOR_GSFE_BORDER, fontsize=8, fontweight="bold")
    _box(ax, 8.0, 2.15, 3.6, 0.7,
         "fcc: FCC_111, FCC_100\nhcp: HCP_basal/prism1w/pyr1w/pyr2\nsrun lmp -in gsfe_template.in\n→ y_gsfe/<lat>/<type>/dump/",
         fc=COLOR_RUN, fontsize=6.2)
    _arrow(ax, 8.0, 2.87, 8.0, 2.52, lw=1.0, color="#7a8aa0")
    _box(ax, 8.0, 1.15, 3.0, 0.6,
         "post/gsfe.py\npost_gsfe\n→ y_post_gsfe.txt",
         fc=COLOR_CODE_BG, fontsize=6.8)
    _arrow(ax, 8.0, 1.80, 8.0, 1.47, lw=1.0, color="#7a8aa0")

    # ---- summary ----
    _box(ax, 5.0, 0.35, 9.6, 0.75,
         "Summary: failed_steps[] 收集 → post/stretch.py | Cij_energy.py | gsfe.py 逐项 [OK]/[FAIL]\n"
         "全部成功 → exit 0 ; 任一失败 → exit 1",
         fc=COLOR_DECISION, ec=COLOR_DECISION_BORDER, fontsize=7.4)
    _arrow(ax, 2.0, 6.85, 3.8, 0.85, lw=1.0, color="#7a8aa0")
    _arrow(ax, 5.0, 3.85, 5.0, 0.78, lw=1.0, color="#7a8aa0")
    _arrow(ax, 8.0, 0.85, 6.2, 0.55, lw=1.0, color="#7a8aa0")

    # stage-to-stage arrows (stretch→Cij→gsfe)
    _arrow(ax, 3.5, 7.15, 3.8, 6.5, lw=1.2, color=COLOR_ARROW,
           label="a0/c0 传递", label_offset=(0.0, 0.15))
    _arrow(ax, 6.5, 4.15, 6.8, 3.5, lw=1.2, color=COLOR_ARROW,
           label="a0/c0 传递", label_offset=(0.0, 0.15))

    # legend strip
    legend_items = [
        ("stretch (平衡晶格)", COLOR_STRETCH),
        ("Cij (弹性常数)", COLOR_CIJ),
        ("gsfe (层错能)", COLOR_GSFE),
        ("LAMMPS 运行", COLOR_RUN),
        ("post 后处理", COLOR_CODE_BG),
    ]
    for i, (txt, col) in enumerate(legend_items):
        lx = 0.2 + i * 1.95
        ax.add_patch(mpatches.Rectangle((lx, 10.78), 0.24, 0.16, fc=col,
                                        ec=COLOR_NODE_BORDER, lw=0.6, zorder=3))
        ax.text(lx + 0.30, 10.86, txt, ha="left", va="center", fontsize=6.4,
                color=COLOR_TEXT, zorder=4)
    _turn_off(ax)


# ---------------------------------------------------------------------------
# Panel 2: template structure comparison table
# ---------------------------------------------------------------------------
def draw_panel_template_table(ax):
    """Table comparing stretch / Cij_energy / gsfe template structure."""
    ax.set_title("面板2 · 三种计算的模板结构对比", loc="left")
    ax.axis("off")

    headers = ["计算类型", "模板文件", "弛豫方式", "变形方式", "输出文件"]
    rows = [
        ("stretch",
         "stretch_template.in",
         "full_relax\n(box/relax iso / couple xy)",
         "面内+面外拉伸\nchange_box x/y/z delta\n±0.003 strain × 101 步",
         "p_post_stretch.txt\n(a0, c0, E, stress)"),
        ("Cij_energy",
         "Cij_energy_template.in",
         "constrained_relax\n(fcc/bcc: 无; hcp: box/relax z)",
         "5 个 Voigt 方向\nC11/C12/C13/C33/C44\ndirn=1..5 × 101 步",
         "y_post_cij_energy.txt\n(Cij 矩阵)"),
        ("gsfe",
         "gsfe_template.in",
         "full_relax → change_box\n(box/relax z, setforce 0 0 NULL)",
         "tilted-cell (bp1/bp2)\nΔxz=bp1·lx·factor\nΔyz=bp2·ly·factor × 21 步",
         "y_post_gsfe.txt\n(GSFE 曲线)"),
    ]

    col_widths = [0.12, 0.20, 0.22, 0.27, 0.19]
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
    table.scale(1.0, 3.2)

    # header
    for j in range(n_cols):
        cell = table[0, j]
        cell.set_facecolor(COLOR_HEADER_BG)
        cell.set_text_props(color="white", fontweight="bold", fontsize=8.2)
        cell.set_edgecolor("white")

    # row tinting by calculation type
    row_colors = {
        1: (COLOR_STRETCH, COLOR_STRETCH_BORDER),
        2: (COLOR_CIJ, COLOR_CIJ_BORDER),
        3: (COLOR_GSFE, COLOR_GSFE_BORDER),
    }
    for i in range(1, n_rows):
        type_cell = table[i, 0]
        fc, ec = row_colors[i]
        type_cell.set_facecolor(fc)
        type_cell.set_text_props(fontweight="bold", fontsize=8, color=COLOR_TEXT)
        for j in range(1, n_cols):
            cell = table[i, j]
            cell.set_text_props(color=COLOR_TEXT, fontsize=6.8)
            if i % 2 == 0:
                cell.set_facecolor("#f7f9fb")
        if i % 2 == 0:
            type_cell.set_facecolor("#e8f0f8")

    ax.text(0.5, 0.02,
            "注: 三种计算共用 general_init.mod / general_potential.mod / general_mass.mod; "
            "Cij 与 gsfe 的 a0/c0 来自 stretch 的 p_post_stretch.txt 第 21 行。",
            ha="center", va="bottom", fontsize=6.8, color=COLOR_SKIP,
            transform=ax.transAxes)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


# ---------------------------------------------------------------------------
# Panel 3: .mod module dependency graph
# ---------------------------------------------------------------------------
def draw_panel_mod_deps(ax):
    """Show how .in files include .mod module files."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 11)
    ax.set_title("面板3 · .mod 模块文件依赖关系 (include 链)", loc="left")

    # three columns: stretch / gsfe / Cij
    col_x = [1.8, 5.0, 8.2]
    col_titles = ["stretch_template.in", "gsfe_template.in", "Cij_energy_template.in"]
    col_colors = [COLOR_STRETCH, COLOR_GSFE, COLOR_CIJ]
    col_borders = [COLOR_STRETCH_BORDER, COLOR_GSFE_BORDER, COLOR_CIJ_BORDER]

    for x, title, fc, ec in zip(col_x, col_titles, col_colors, col_borders):
        _box(ax, x, 10.3, 2.8, 0.55, title,
             fc=fc, ec=ec, fontsize=8, fontweight="bold")

    # ---- stretch deps ----
    stretch_deps = [
        ("general_init.mod\nunits metal / boundary", COLOR_NODE),
        ("stretch_model.mod\nlattice (aa, lat)", COLOR_NODE),
        ("general_potential.mod\npair_style/coeff + mass", COLOR_NODE),
        ("stretch_full_relax.mod\nbox/relax iso / couple xy", COLOR_NODE),
        ("stretch.mod\nchange_box loop (101 步)", COLOR_RUN),
    ]
    for i, (txt, fc) in enumerate(stretch_deps):
        y = 9.1 - i * 1.55
        _box(ax, col_x[0], y, 2.8, 1.0, txt, fc=fc, fontsize=6.8)
        _arrow(ax, col_x[0], y + 0.55, col_x[0], y + 0.50, lw=0.9, color="#7a8aa0")
    # first arrow from .in to first dep
    _arrow(ax, col_x[0], 10.0, col_x[0], 9.62, lw=0.9, color="#7a8aa0")

    # ---- gsfe deps ----
    gsfe_deps = [
        ("general_init.mod\nunits metal / boundary", COLOR_NODE),
        ("gsfe_model.py\ncreate_gsfe_model → data.ini\n(python_path_template)", COLOR_NODE),
        ("general_potential.mod\npair_style/coeff + mass", COLOR_NODE),
        ("minimize (full relax)\nwrite_restart restart.equil", COLOR_NODE),
        ("gsfe.mod\nchange_box tilt loop (21 步)\nbp1/bp2 shear", COLOR_RUN),
    ]
    for i, (txt, fc) in enumerate(gsfe_deps):
        y = 9.1 - i * 1.55
        _box(ax, col_x[1], y, 2.8, 1.0, txt, fc=fc, fontsize=6.6)
        _arrow(ax, col_x[1], y + 0.55, col_x[1], y + 0.50, lw=0.9, color="#7a8aa0")
    _arrow(ax, col_x[1], 10.0, col_x[1], 9.62, lw=0.9, color="#7a8aa0")

    # ---- Cij deps ----
    cij_deps = [
        ("general_init.mod\nunits metal / boundary", COLOR_NODE),
        ("stretch_model.mod\nlattice (aa, lat)", COLOR_NODE),
        ("general_potential.mod\npair_style/coeff + mass", COLOR_NODE),
        ("stretch_constrained_relax.mod\nfcc/bcc: 无; hcp: box/relax z", COLOR_NODE),
        ("Cij_energy.mod\ndirn=1..5 × 101 步\n+ Cij_energy_output.mod", COLOR_RUN),
    ]
    for i, (txt, fc) in enumerate(cij_deps):
        y = 9.1 - i * 1.55
        _box(ax, col_x[2], y, 2.8, 1.0, txt, fc=fc, fontsize=6.6)
        _arrow(ax, col_x[2], y + 0.55, col_x[2], y + 0.50, lw=0.9, color="#7a8aa0")
    _arrow(ax, col_x[2], 10.0, col_x[2], 9.62, lw=0.9, color="#7a8aa0")

    # shared-modules note at bottom
    _box(ax, 5.0, 0.35, 9.6, 0.8,
         "共用模块: general_init.mod (units/boundary/atom_style) · general_mass.mod (运行时生成) ·\n"
         "general_potential.mod (含 general_mass.mod) · general_output.mod / general_structural_info.mod",
         fc=COLOR_CODE_BG, fontsize=6.8)

    _turn_off(ax)


# ---------------------------------------------------------------------------
# Panel 4: GSFE slip-system registry table
# ---------------------------------------------------------------------------
def draw_panel_gsfe_table(ax):
    """Table of gsfe_type → bp1/bp2 shear offsets."""
    ax.set_title("面板4 · GSFE 滑移系类型 (gsfe_type → bp1/bp2 剪切偏移)", loc="left")
    ax.axis("off")

    headers = ["gsfe_type", "晶格", "bp1 (Δxz)", "bp2 (Δyz)", "LAMMPS 表达式", "说明"]
    rows = [
        ("FCC_111", "fcc", "-0.5", "1/3", "1.0/3.0", "{111}⟨112⟩ 主滑移系"),
        ("FCC_100", "fcc", "0.5", "0.5", "0.5", "{100}⟨110⟩"),
        ("HCP_basal", "hcp", "-0.5", "1/3", "1.0/3.0", "(0001) 基面"),
        ("HCP_prism1w", "hcp", "0.5", "0", "0.0", "(10-10) 柱面 I"),
        ("HCP_pyr1w", "hcp", "0", "0.5", "0.5", "(10-11) 锥面 I"),
        ("HCP_pyr2", "hcp", "0", "0.5", "0.5", "(10-11) 锥面 II"),
    ]

    col_widths = [0.16, 0.08, 0.11, 0.11, 0.16, 0.38]
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
    table.scale(1.0, 1.85)

    # header
    for j in range(n_cols):
        cell = table[0, j]
        cell.set_facecolor(COLOR_HEADER_BG)
        cell.set_text_props(color="white", fontweight="bold", fontsize=8)
        cell.set_edgecolor("white")

    # body: color fcc vs hcp rows
    for i in range(1, n_rows):
        lat = rows[i - 1][1]
        type_cell = table[i, 0]
        lat_cell = table[i, 1]
        if lat == "fcc":
            band = COLOR_STRETCH
        else:
            band = COLOR_GSFE
        type_cell.set_facecolor(band)
        type_cell.set_text_props(fontweight="bold", fontsize=7.6)
        lat_cell.set_text_props(fontweight="bold", color=COLOR_TEXT)
        for j in range(2, n_cols):
            cell = table[i, j]
            cell.set_text_props(color=COLOR_TEXT, fontsize=7.2)
            if j == 4:  # LAMMPS expression column → monospace feel
                cell.set_text_props(color=COLOR_WARN, fontweight="bold", fontsize=7.0)
            if i % 2 == 0:
                cell.set_facecolor("#f7f9fb")
        if i % 2 == 0:
            type_cell.set_facecolor("#e8f0f8" if lat == "fcc" else "#f3eefb")
            lat_cell.set_facecolor("#f7f9fb")

    ax.text(0.5, 0.02,
            "注: bp1/bp2 与 DFT vasp_utils/.../pei_vasp_run_gsfe 完全一致; "
            "Δxz = bp1·lx·factor, Δyz = bp2·ly·factor, factor: 0 → 1 (21 步)。",
            ha="center", va="bottom", fontsize=6.8, color=COLOR_SKIP,
            transform=ax.transAxes)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


# ---------------------------------------------------------------------------
# Panel 5: sed template substitution flow
# ---------------------------------------------------------------------------
def draw_panel_sed_flow(ax):
    """Show how pei_lmp_run_properties uses sed to substitute placeholders."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 11)
    ax.set_title("面板5 · sed 模板替换流程 (pei_lmp_run_properties)", loc="left")

    # left: runner invokes sed
    _box(ax, 2.0, 10.2, 3.6, 0.6,
         "pei_lmp_run_properties\nfor lat in fcc bcc hcp",
         fc=COLOR_HEADER_BG, text_color="white", fontweight="bold", fontsize=8)

    # placeholder table on the left
    _box(ax, 2.0, 8.9, 3.6, 1.7,
         "占位符 → 实际值 (按 lat / gsfe_type):\n"
         "aa_template  → ini_aa / a0\n"
         "  fcc:4.08  bcc:3.20  hcp:2.88\n"
         "lat_template → latnum\n"
         "  fcc:2  bcc:3  hcp:1",
         fc=COLOR_CODE_BG, fontsize=6.8)
    _arrow(ax, 2.0, 9.9, 2.0, 9.78, lw=1.0, color="#7a8aa0")

    _box(ax, 2.0, 7.1, 3.6, 1.5,
         "势函数占位符:\n"
         "pair_style_template → eam / ...\n"
         "pair_coeff_template → Au_u3.eam / ...\n"
         "(由 $1 / $2 参数传入)",
         fc=COLOR_CODE_BG, fontsize=6.8)
    _arrow(ax, 2.0, 8.03, 2.0, 7.87, lw=1.0, color="#7a8aa0")

    _box(ax, 2.0, 5.3, 3.6, 1.5,
         "gsfe 占位符 (仅 gsfe 阶段):\n"
         "gsfe_type_template → FCC_111 / ...\n"
         "gsfe_bp1_template  → -0.5 / 0.5 / ...\n"
         "gsfe_bp2_template  → 1.0/3.0 / 0.5 / ...\n"
         "cc_template → c0 (hcp)",
         fc=COLOR_CODE_BG, fontsize=6.6)
    _arrow(ax, 2.0, 6.33, 2.0, 6.07, lw=1.0, color="#7a8aa0")

    # middle: sed command
    _code(ax, 5.5, 8.9, 2.6, 1.0,
          "sed -i\n  s/aa_template/${aa}/g\n  s/lat_template/${latnum}/g\n  s|pair_style_template|...|g",
          fontsize=6.6)
    _arrow(ax, 3.85, 8.9, 4.2, 8.9, lw=1.2, color=COLOR_ARROW,
           label="sed -i", label_offset=(0.0, 0.18))

    _code(ax, 5.5, 7.1, 2.6, 1.0,
          "sed -i\n  s|pair_coeff_template|...|g\n  printf '%b' $mass_content\n  > general_mass.mod",
          fontsize=6.6)
    _arrow(ax, 3.85, 7.1, 4.2, 7.1, lw=1.2, color=COLOR_ARROW)

    _code(ax, 5.5, 5.3, 2.6, 1.0,
          "sed -i\n  s|gsfe_type_template|${gsfe_type}|g\n  s|gsfe_bp1_template|${bp1}|g\n  s|gsfe_bp2_template|${bp2}|g",
          fontsize=6.4)
    _arrow(ax, 3.85, 5.3, 4.2, 5.3, lw=1.2, color=COLOR_ARROW)

    # right: target files
    targets = [
        (8.5, 8.9, "stretch_model.mod\n→ variable aa / lat"),
        (8.5, 7.1, "general_potential.mod\n→ pair_style / pair_coeff"),
        (8.5, 5.3, "gsfe_model.py + gsfe.mod\n→ gsfe_type / bp1 / bp2"),
    ]
    for x, y, txt in targets:
        _box(ax, x, y, 2.6, 1.0, txt,
             fc=COLOR_NODE, ec=COLOR_OK, fontsize=6.8)
    _arrow(ax, 6.85, 8.9, 7.2, 8.9, lw=1.0, color=COLOR_OK)
    _arrow(ax, 6.85, 7.1, 7.2, 7.1, lw=1.0, color=COLOR_OK)
    _arrow(ax, 6.85, 5.3, 7.2, 5.3, lw=1.0, color=COLOR_OK)

    # bottom: mass_content generation
    _box(ax, 5.0, 3.5, 9.6, 1.0,
         "general_mass.mod 运行时生成 (非 sed):\n"
         "printf '%b\\n' \"$mass_content\" > general_mass.mod\n"
         "mass_content 由 runner 第 4 参数传入, 默认 'mass 1 196.97'",
         fc=COLOR_DECISION, ec=COLOR_DECISION_BORDER, fontsize=7.0)

    # final: srun lmp
    _box(ax, 5.0, 1.8, 9.6, 1.0,
         "替换完成后:\n"
         "run_lmp_retry -in <xxx>_template.in -screen none\n"
         "(srun step 创建失败时重试, 默认 LMP_MAX_TRY=99, 间隔 5s)",
         fc=COLOR_RUN, fontsize=7.2)
    _arrow(ax, 5.0, 3.0, 5.0, 2.32, lw=1.0, color="#7a8aa0")

    # summary note
    _box(ax, 5.0, 0.45, 9.6, 0.75,
         "failed_steps[] 收集所有 lmp/post 失败 → Summary 逐项打印 → exit 0 (全成功) / exit 1 (有失败)",
         fc=COLOR_PANEL_BG, ec="#d1d5db", fontsize=7.0)
    _arrow(ax, 5.0, 1.3, 5.0, 0.83, lw=1.0, color="#7a8aa0")

    _turn_off(ax)


# ---------------------------------------------------------------------------
# Figure assembly
# ---------------------------------------------------------------------------
def render_overview(path_image: Path) -> None:
    """Build the 5-panel overview figure and save it to path_image."""
    _set_rcparams()
    fig = plt.figure(figsize=(16, 12), dpi=150)
    # 3 rows x 2 cols, last cell (panel 5) spans bottom-right
    gs = fig.add_gridspec(
        3, 2, hspace=0.34, wspace=0.12,
        left=0.03, right=0.98, top=0.95, bottom=0.03,
    )
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    ax5 = fig.add_subplot(gs[2, :])

    draw_panel_runner_flow(ax1)
    draw_panel_template_table(ax2)
    draw_panel_mod_deps(ax3)
    draw_panel_gsfe_table(ax4)
    draw_panel_sed_flow(ax5)

    # suptitle
    fig.suptitle(
        "lmp_utils 子包功能总览",
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
    path_image = path_output / "lmp_utils_overview.png"

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
        description="Render the lmp_utils overview figure (LAMMPS-free).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("docs/_build/example-lmp-utils"),
        help="Output directory for the rendered PNG (default: docs/_build/example-lmp-utils)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary = run_example(args.output)
    # deterministic assertions
    assert summary["size_kb"] > 10.0, "PNG smaller than 10 KB: " + str(summary)
    assert summary["pixel_deviation"] > 0.002, "PNG effectively blank: " + str(summary)
    print("OK: lmp_utils overview figure verified non-blank.")


if __name__ == "__main__":
    main()
