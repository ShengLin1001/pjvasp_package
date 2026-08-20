#!/usr/bin/env python3
"""Plot gallery demo for ``mymetal.universal.plot``.

This tutorial is VASP-free. It demonstrates the public functions of every
module under :mod:`mymetal.universal.plot` using synthetic data and the
package's own styling entry points (:func:`my_plot` /
:func:`general_set_all_rcParams`). Each section builds a small, deterministic
figure and writes a PNG to the output directory so the Sphinx gallery pages
can embed it.

Covered modules:

* :mod:`mymetal.universal.plot.general` — style / annotation helpers
* :mod:`mymetal.universal.plot.plot` — figure / canvas creators
* :mod:`mymetal.universal.plot.workflow` — workflow post-processing curves
* :mod:`mymetal.universal.plot.energy` — energy component decomposition
* :mod:`mymetal.universal.plot.atominfo` — layer distance / z-position / RDF
* :mod:`mymetal.universal.plot.n2p2` — n2p2 training diagnostics
* :mod:`mymetal.universal.plot.plotting` — pymatgen-style pretty plots
* :mod:`mymetal.universal.plot.render` — OVITO pipeline rendering (illustrated)
* :mod:`mymetal.universal.plot.ppt` — PowerPoint slide export (illustrated)
* :mod:`mymetal.universal.plot.oldplotdos` — DOS / IDOS plotting helpers
"""

from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from mymetal.universal.plot.general import (  # noqa: E402
    add_arrow,
    add_circle_number,
    add_color_band,
    general_add_vlines_hlines,
    general_adjust_text,
    general_modify_legend,
    general_modify_line,
    generate_gradient_colors,
    general_margin_bin,
)
from mymetal.universal.plot.plot import my_plot, my_plot_colorbar  # noqa: E402


# Hardcoded synthetic seeds for determinism.
SEED = 20260820


# ---------------------------------------------------------------------------
# Output directory helpers
# ---------------------------------------------------------------------------
def _assert_non_blank(path_png: Path) -> None:
    """Re-read a saved PNG and assert it is neither tiny nor blank-white."""
    if not path_png.is_file() or path_png.stat().st_size < 1_000:
        raise AssertionError("image was not created or too small: " + str(path_png))
    rgb = plt.imread(path_png)[..., :3]
    if float(np.mean(np.abs(rgb - 1.0))) < 0.002:
        raise AssertionError("image is effectively blank: " + str(path_png))


# ---------------------------------------------------------------------------
# general.py
# ---------------------------------------------------------------------------
def render_general(path_out: Path) -> Path:
    """Demonstrate style/annotation helpers from ``general.py``."""
    np.random.seed(SEED)

    fig, axes = my_plot(fig_subp=[2, 2], fig_sharex=False)
    # axes is a 2x2 ndarray; flatten for convenience.
    axa = axes.flatten()

    # (0) add_color_band + generate_gradient_colors
    x = np.linspace(0, 1, 50)
    y = np.sin(2 * np.pi * x)
    axa[0].plot(x, y, marker="o", label="sin")
    gradient = generate_gradient_colors(
        if_cmap_color=True, if_reshape=True,
        reshape_M_N_L=[1, 1000, 4], if_reverse=True,
    )
    add_color_band(ax=axa[0], extent=[0.25, 0.45, -1.2, 1.2], gradient=gradient, alpha=0.4)
    axa[0].set_xlabel("x")
    axa[0].set_ylabel("signal")
    general_modify_legend(axa[0].legend(loc="upper right"))

    # (1) add_circle_number + add_arrow
    x2 = np.linspace(0, 5, 20)
    y2 = 0.3 * x2 ** 2
    axa[1].plot(x2, y2, "-o", label="growth")
    add_circle_number(ax=axa[1], positions=[3.0, 0.3 * 9.0], number=1, radiusx_ratio=0.05)
    add_arrow(ax=axa[1], text="peak", start=[3.2, 3.5], end=[4.2, 7.5],
              arrowprops=dict(arrowstyle="->", color="red"))
    axa[1].set_xlabel("step")
    axa[1].set_ylabel("value")
    general_modify_legend(axa[1].legend(loc="upper left"))

    # (2) general_add_vlines_hlines + general_modify_line
    x3 = np.linspace(0, 6, 30)
    y3a = np.sin(x3)
    y3b = 0.5 + 0.3 * np.cos(x3)
    axa[2].plot(x3, y3a, "--", label="raw sin")
    axa[2].plot(x3, y3b, ":", label="shift cos")
    general_add_vlines_hlines(axa[2], vlines=[2.0, 4.0], hlines=[0.5])
    general_modify_line(axa[2], if_change_color=True)
    axa[2].set_xlabel("x")
    axa[2].set_ylabel("amplitude")
    general_modify_legend(axa[2].legend(loc="upper right"))

    # (3) general_adjust_text (avoid overlapping text)
    xs = np.random.rand(8)
    ys = np.random.rand(8)
    axa[3].scatter(xs, ys, c="steelblue", s=80)
    texts = [axa[3].text(xs[i], ys[i], "P" + str(i), fontsize=14) for i in range(len(xs))]
    general_adjust_text(texts=texts, ax=axa[3], x=list(xs), y=list(ys),
                        ensure_inside_axes=True, iter_lim=200)
    general_margin_bin(axa[3], x_margin=0.15, y_margin=0.15)
    axa[3].set_xlabel("u")
    axa[3].set_ylabel("v")

    fig.savefig(path_out, bbox_inches="tight")
    plt.close(fig)
    _assert_non_blank(path_out)
    return path_out


# ---------------------------------------------------------------------------
# plot.py
# ---------------------------------------------------------------------------
def render_plot(path_out: Path) -> Path:
    """Demonstrate the figure creators in ``plot.py``."""
    fig, axes = my_plot(fig_subp=[1, 2], fig_sharex=False)
    x = np.linspace(0, 2 * np.pi, 100)
    axes[0].plot(x, np.sin(x), "-o", label="sin(x)")
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("y")
    general_modify_legend(axes[0].legend(loc="upper right"))
    axes[1].plot(x, np.cos(x), "-s", label="cos(x)", color="darkred")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("y")
    general_modify_legend(axes[1].legend(loc="upper right"))

    fig.savefig(path_out, bbox_inches="tight")
    plt.close(fig)
    _assert_non_blank(path_out)
    return path_out


def render_colorbar(path_out: Path) -> Path:
    """Demonstrate ``my_plot_colorbar``.

    ``my_plot_colorbar`` returns ``(fig, (ax_main, ax_cbar))``: the second
    element is a tuple of the main axis and the dedicated colorbar axis.
    """
    fig, axes = my_plot_colorbar(grid=False)
    ax_main, ax_cbar = axes[0], axes[1]
    data = np.random.RandomState(SEED).rand(20, 20)
    im = ax_main.imshow(data, origin="lower", cmap="viridis")
    fig.colorbar(im, cax=ax_cbar)
    fig.savefig(path_out, bbox_inches="tight")
    plt.close(fig)
    _assert_non_blank(path_out)
    return path_out


# ---------------------------------------------------------------------------
# workflow.py
# ---------------------------------------------------------------------------
def render_convergence(path_out: Path) -> Path:
    """Demonstrate ``my_plot_convergence`` with synthetic energy-vs-cutoff data.

    The function computes ``y - y[-1]`` internally when ``if_difference=True``;
    pass ``y`` as a numpy array (not a list) so the subtraction works. The
    function saves the figure itself when ``if_save=True``.
    """
    from mymetal.universal.plot.workflow import my_plot_convergence

    encuts = [300, 400, 500, 600, 700, 800]
    energies = np.array([-3.214, -3.217, -3.218, -3.2185, -3.2186, -3.2187])
    my_plot_convergence(
        x=encuts, y=energies, encuts=True, if_save=True,
        if_difference=True, savefile=str(path_out),
    )
    plt.close("all")
    _assert_non_blank(path_out)
    return path_out


def render_relax_convergence(path_out: Path) -> Path:
    """Demonstrate ``my_plot_relax_convergence``.

    The function saves and optionally closes the figure itself, so we pass
    ``if_save=True`` and ``savefile=path_out`` rather than saving afterwards.
    """
    from mymetal.universal.plot.workflow import my_plot_relax_convergence

    n = 20
    lframe = list(range(1, n + 1))
    base = -10.0
    rng = np.random.RandomState(SEED)
    lenergy = [base + 0.5 * np.exp(-i / 5) + 0.01 * rng.rand() for i in range(n)]
    lforce = [2.0 * np.exp(-i / 4) + 0.05 * rng.rand() for i in range(n)]
    my_plot_relax_convergence(
        lframe=lframe, lenergy=lenergy, lforce=lforce,
        natoms=4, ediffg=-0.01, isif=2, yscale="log",
        if_save=True, savefile=str(path_out), if_close=False,
    )
    plt.close("all")
    _assert_non_blank(path_out)
    return path_out


def render_kpar_ncore(path_out: Path) -> Path:
    """Demonstrate ``my_plot_kpar_ncore``.

    The function expects ``dict_time`` and ``dict_delta_energy`` keyed by
    ``(kpar, ncore)`` tuples, and saves the figure itself when
    ``if_save=True``.
    """
    from mymetal.universal.plot.workflow import my_plot_kpar_ncore

    lkpar = [4, 8, 16, 32, 64]
    lncore = [16, 8, 4, 2, 1]
    # Elapsed time (min) and delta-energy (meV/atom) keyed by (kpar, ncore).
    times = [1200, 620, 340, 200, 180, 800, 520, 380, 320, 350,
             600, 430, 350, 320, 360, 520, 400, 350, 340, 380,
             480, 390, 350, 340, 390]
    deltas = [0.0, 0.2, 0.5, 1.1, 2.4, 0.0, 0.1, 0.3, 0.6, 1.0,
              0.0, 0.1, 0.2, 0.4, 0.7, 0.0, 0.05, 0.2, 0.3, 0.5,
              0.0, 0.03, 0.1, 0.2, 0.3]
    dict_time = {}
    dict_delta_energy = {}
    idx = 0
    for kpar in lkpar:
        for ncore in lncore:
            dict_time[(kpar, ncore)] = times[idx]
            dict_delta_energy[(kpar, ncore)] = deltas[idx]
            idx += 1
    my_plot_kpar_ncore(
        dict_time=dict_time, dict_delta_energy=dict_delta_energy,
        lkpar=lkpar, lncore=lncore,
        if_save=True, savefile=str(path_out), if_close=False,
    )
    plt.close("all")
    _assert_non_blank(path_out)
    return path_out


def render_stretch(path_out: Path) -> Path:
    """Demonstrate ``my_plot_stretch``.

    The function multiplies ``rvectors_ref`` by the fitted equilibrium factor,
    so a non-None reference row vector must be supplied.
    """
    from mymetal.universal.plot.workflow import my_plot_stretch

    jobn = ["0.97", "0.98", "0.99", "1.00", "1.01", "1.02", "1.03"]
    factors = [float(j) for j in jobn]
    e0 = -12.5
    Etot = [e0 + 10.0 * (f - 1.0) ** 2 for f in factors]
    lca = [1.621 * f for f in factors]
    rvectors_ref = np.array([2.556, 2.556, 4.073])  # synthetic HCP-like row vector
    my_plot_stretch(
        jobn=jobn, Etot=Etot, natoms=4, stretch_type="c",
        rvectors_ref=rvectors_ref, lca=lca,
        if_save=True, savefile=str(path_out),
    )
    plt.close("all")
    _assert_non_blank(path_out)
    return path_out


# ---------------------------------------------------------------------------
# energy.py
# ---------------------------------------------------------------------------
def render_energy_components(path_out: Path) -> Path:
    """Demonstrate ``my_plot_energy_components`` with synthetic OUTCAR-like dicts.

    The function iterates over the full fixed tag list (TOTEN, PSCENC, TEWEN,
    DENC, EXHF, XCENC, Double1, Double2, Double, EENTRO, EBANDS, EATOM,
    Ediel_sol), so every dict must contain all of these keys (lower-cased).
    """
    from mymetal.universal.plot.energy import my_plot_energy_components

    dirs = [0.0, 0.1, 0.2, 0.3, 0.4]
    # Full tag set expected by my_plot_energy_components (lower-cased from tagsref).
    tags = ["toten", "pscenc", "tewen", "denc", "exhf", "xcenc",
            "double1", "double2", "double", "eentro", "ebands",
            "eatom", "ediel_sol"]
    base = {"toten": -10.0, "pscenc": 5.0, "tewen": 12.0, "denc": -18.0,
            "exhf": -1.0, "xcenc": -4.0, "double1": 0.5, "double2": -0.5,
            "double": 0.0, "eentro": -0.1, "ebands": 6.0, "eatom": -2.0,
            "ediel_sol": 0.0}

    def build_dict(scale: float) -> dict:
        d = {"directory": dirs}
        for t in tags:
            d[t] = np.array([base[t] * (1.0 + 0.05 * x * scale) for x in dirs])
        return d

    dic_list = [build_dict(0.5), build_dict(1.0)]
    label_list = ["set A", "set B"]
    my_plot_energy_components(
        dic_list=dic_list, label_list=label_list, atoms_number=4,
        if_tight_layout=False, save_path=str(path_out), poly_fit_e=2,
    )
    plt.close("all")
    _assert_non_blank(path_out)
    return path_out


# ---------------------------------------------------------------------------
# atominfo.py
# ---------------------------------------------------------------------------
def _build_synthetic_slab() -> "Atoms":  # type: ignore[name-defined]
    """Build a thin synthetic slab (VASP-free) for atominfo demos."""
    from ase import Atoms
    from ase.build import bulk

    atoms = bulk("Cu", "fcc", a=3.61, cubic=True)
    atoms = atoms.repeat((2, 2, 6))
    atoms.cell[2, 2] += 6.0  # add vacuum
    atoms.wrap()
    return atoms


def render_interlayer_distance(path_out: Path) -> Path:
    """Demonstrate ``my_plot_interlayer_distance``."""
    from mymetal.universal.plot.atominfo import my_plot_interlayer_distance

    slab = _build_synthetic_slab()
    fig, ax = my_plot_interlayer_distance(
        atoms=slab, if_plot=True, if_save=True, if_save_txt=False,
        save_plot_path=str(path_out),
    )
    fig.savefig(path_out, bbox_inches="tight")
    plt.close(fig)
    _assert_non_blank(path_out)
    return path_out


def render_zpositions(path_out: Path) -> Path:
    """Demonstrate ``my_plot_zpositions``."""
    from mymetal.universal.plot.atominfo import my_plot_zpositions

    slab = _build_synthetic_slab()
    fig, axes = my_plot_zpositions(
        atoms=slab, if_save=True, save_path=str(path_out),
    )
    fig.savefig(path_out, bbox_inches="tight")
    plt.close(fig)
    _assert_non_blank(path_out)
    return path_out


def render_rdf(path_out: Path) -> Path:
    """Demonstrate ``my_plot_rdf``."""
    from mymetal.universal.plot.atominfo import my_plot_rdf

    atoms = _build_synthetic_slab()
    fig, axes = my_plot_rdf(atoms=atoms, bins=30, cutoff=6.0, title="Synthetic Cu RDF")
    fig.savefig(path_out, bbox_inches="tight")
    plt.close(fig)
    _assert_non_blank(path_out)
    return path_out


# ---------------------------------------------------------------------------
# n2p2.py
# ---------------------------------------------------------------------------
def render_learning_curve(path_out: Path) -> Path:
    """Demonstrate ``my_plot_learning_curve``."""
    from mymetal.universal.plot.n2p2 import my_plot_learning_curve

    epochs = np.arange(1, 31)
    e_rmse = 10.0 * np.exp(-0.15 * epochs) + 0.5
    f_rmse = 120.0 * np.exp(-0.12 * epochs) + 5.0
    my_plot_learning_curve(
        epochs=epochs, e_rmse_mev=e_rmse, f_rmse_mev=f_rmse,
        file_path=str(path_out), label="Train",
    )
    _assert_non_blank(path_out)
    return path_out


def render_compare(path_out: Path) -> Path:
    """Demonstrate ``my_plot_compare``."""
    from mymetal.universal.plot.n2p2 import my_plot_compare

    rng = np.random.RandomState(SEED)
    n = 200
    e_ref = rng.normal(-3.5, 0.3, n)
    e_nnp = e_ref + rng.normal(0, 0.05, n)
    tag_e = np.array(["bulk"] * 100 + ["surface"] * 100)
    f_ref = rng.normal(0, 1.5, n * 3)
    f_nnp = f_ref + rng.normal(0, 0.2, n * 3)
    tag_f = np.array(["bulk"] * 300 + ["surface"] * 300)
    my_plot_compare(
        e_ref=e_ref, e_nnp=e_nnp, tag_e=tag_e,
        f_ref=f_ref, f_nnp=f_nnp, tag_f=tag_f,
        file_path=str(path_out),
        text_e="RMSE_E = 5.0 meV/atom", text_f="RMSE_F = 20.0 meV/A",
    )
    _assert_non_blank(path_out)
    return path_out


def render_rmse_by_tag(path_out: Path) -> Path:
    """Demonstrate ``my_plot_rmse_by_tag``."""
    from mymetal.universal.plot.n2p2 import my_plot_rmse_by_tag
    import pandas as pd

    df = pd.DataFrame({
        "tag": ["bulk", "surface", "vacancy", "TOTAL"],
        "E_RMSE_meV/at": [3.2, 5.1, 8.7, 5.6],
        "F_RMSE_meV/A": [45.0, 80.0, 130.0, 85.0],
    })
    my_plot_rmse_by_tag(df=df, file_path=str(path_out))
    _assert_non_blank(path_out)
    return path_out


def render_epoch_rmse(path_out: Path) -> Path:
    """Demonstrate ``my_plot_epoch_rmse``."""
    from mymetal.universal.plot.n2p2 import my_plot_epoch_rmse
    import pandas as pd

    epochs = np.arange(1, 21)
    df = pd.DataFrame({
        "epoch": epochs,
        "E_RMSE_meV/at": 10.0 * np.exp(-0.15 * epochs) + 0.5,
        "F_RMSE_meV/A": 120.0 * np.exp(-0.12 * epochs) + 5.0,
    })
    my_plot_epoch_rmse(df=df, file_path=str(path_out))
    _assert_non_blank(path_out)
    return path_out


# ---------------------------------------------------------------------------
# plotting.py
# ---------------------------------------------------------------------------
def render_pretty_plot(path_out: Path) -> Path:
    """Demonstrate ``pretty_plot`` and ``pretty_polyfit_plot``."""
    from mymetal.universal.plot.plotting import pretty_plot, pretty_polyfit_plot

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    x = np.linspace(0, 10, 30)
    y = 0.5 * x + 0.2 + 0.3 * np.random.RandomState(SEED).randn(30)

    plt.sca(axes[0])
    pretty_plot(width=6, ax=axes[0])
    axes[0].plot(x, y, "o", label="data")
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("y")
    axes[0].legend()

    plt.sca(axes[1])
    pretty_polyfit_plot(x, y, deg=1, xlabel="x", ylabel="y")
    axes[1].set_title("polyfit deg=1")

    fig.savefig(path_out, bbox_inches="tight")
    plt.close(fig)
    _assert_non_blank(path_out)
    return path_out


def render_periodic_heatmap(path_out: Path) -> Path:
    """Demonstrate ``periodic_table_heatmap`` (matplotlib backend)."""
    from mymetal.universal.plot.plotting import periodic_table_heatmap

    data = {
        "Fe": 211, "Cu": 140, "Na": 9.5, "K": 3.5,
        "W": 310, "Au": 180, "Mg": 45, "Al": 76,
    }
    periodic_table_heatmap(
        elemental_data=data, cbar_label="Bulk modulus (GPa)",
        pymatviz=False, show_plot=False,
    )
    plt.savefig(path_out, bbox_inches="tight")
    plt.close()
    _assert_non_blank(path_out)
    return path_out


def render_van_arkel(path_out: Path) -> Path:
    """Demonstrate ``van_arkel_triangle``."""
    from mymetal.universal.plot.plotting import van_arkel_triangle

    fig, ax = plt.subplots(figsize=(8, 7))
    materials = [["Na", "Cl"], ["Mg", "O"], ["Al", "O"], ["Ga", "As"], ["Cu", "Br"]]
    van_arkel_triangle(materials, annotate=True)
    fig.savefig(path_out, bbox_inches="tight")
    plt.close(fig)
    _assert_non_blank(path_out)
    return path_out


# ---------------------------------------------------------------------------
# render.py, ppt.py, oldplotdos.py — illustrated (the real functions need
# OVITO / PowerPoint / a vasprun.xml file that we deliberately avoid here).
# ---------------------------------------------------------------------------
def render_render_info(path_out: Path) -> Path:
    """Illustrate ``my_render`` (OVITO) with a schematic panel.

    The real ``my_render`` requires an OVITO pipeline and a structure file,
    which is outside the VASP-free scope of this gallery. The panel below
    summarises the call signature and the typical render workflow.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.axis("off")
    text = (
        "my_render(pipeline, imagefile, size, renderer, camera_dir, viewtype)\n\n"
        "Typical workflow:\n"
        "  pipeline = import_file('POSCAR')\n"
        "  pipeline.modifiers.append(CommonNeighborAnalysisModifier())\n"
        "  my_render(pipeline=pipeline, imagefile='render.png',\n"
        "            size=(3200, 2400),\n"
        "            renderer=TachyonRenderer(shadows=True),\n"
        "            camera_dir=(1, 1, 1))\n\n"
        "Requires the optional ovito package."
    )
    ax.text(0.05, 0.5, text, fontsize=14, family="monospace",
            verticalalignment="center", transform=ax.transAxes,
            bbox=dict(facecolor="#f8f9fa", edgecolor="black", boxstyle="Round,pad=0.5"))
    fig.savefig(path_out, bbox_inches="tight")
    plt.close(fig)
    _assert_non_blank(path_out)
    return path_out


def render_ppt_info(path_out: Path) -> Path:
    """Illustrate ``ppt2picture`` with a schematic panel.

    The real ``ppt2picture`` drives the Windows PowerPoint COM interface,
    which is outside the VASP-free scope of this gallery. The panel summarises
    the call and the export pipeline.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.axis("off")
    text = (
        "ppt2picture(path_ppt, path_output, style='png', dpi=600,\n"
        "            padding_mm=1.0, crop_threshold=8, overwrite=False)\n\n"
        "Exports every slide in a .pptx file as a cropped image.\n"
        "Requires Microsoft PowerPoint and pywin32 on Windows.\n"
        "Supported styles: png, jpg, jpeg, bmp, gif, tif, tiff."
    )
    ax.text(0.05, 0.5, text, fontsize=14, family="monospace",
            verticalalignment="center", transform=ax.transAxes,
            bbox=dict(facecolor="#f8f9fa", edgecolor="black", boxstyle="Round,pad=0.5"))
    fig.savefig(path_out, bbox_inches="tight")
    plt.close(fig)
    _assert_non_blank(path_out)
    return path_out


def render_dos(path_out: Path) -> Path:
    """Demonstrate DOS-style plotting with synthetic data.

    The real ``my_plot_complete_dos`` reads a ``vasprun.xml`` file, which is
    outside the VASP-free scope of this gallery. Here we build a synthetic
    Gaussian DOS and plot it with the package's styling entry points.
    """
    from scipy.special import wofz

    energies = np.linspace(-6, 6, 400)
    ef = 0.0
    e = energies - ef
    sigma = 0.4
    dos_total = np.exp(-(e ** 2) / (2 * sigma ** 2)) / (sigma * np.sqrt(2 * np.pi))
    dos_s = 0.6 * np.exp(-((e - 0.5) ** 2) / (2 * sigma ** 2)) / (sigma * np.sqrt(2 * np.pi))
    dos_p = 0.4 * np.exp(-((e + 0.5) ** 2) / (2 * sigma ** 2)) / (sigma * np.sqrt(2 * np.pi))

    fig, axes = my_plot(fig_subp=[1, 2], fig_sharex=False)
    axes[0].plot(energies, dos_total, "-", label="Total")
    axes[0].fill_between(energies, 0, dos_total, alpha=0.3)
    axes[0].axvline(ef, color="gray", linestyle="--")
    axes[0].set_xlabel(r"$E - E_f$ (eV)")
    axes[0].set_ylabel("DOS (a.u.)")
    axes[0].set_title("Total DOS")
    general_modify_legend(axes[0].legend(loc="upper right"))

    axes[1].plot(energies, dos_s, "-", label="s")
    axes[1].plot(energies, dos_p, "-", label="p")
    axes[1].axvline(ef, color="gray", linestyle="--")
    axes[1].set_xlabel(r"$E - E_f$ (eV)")
    axes[1].set_ylabel("DOS (a.u.)")
    axes[1].set_title("Projected DOS")
    general_modify_legend(axes[1].legend(loc="upper right"))

    fig.savefig(path_out, bbox_inches="tight")
    plt.close(fig)
    _assert_non_blank(path_out)
    return path_out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
SECTION_PLAN = [
    ("general", render_general, "plot_gallery_general.png"),
    ("plot", render_plot, "plot_gallery_plot.png"),
    ("plot_colorbar", render_colorbar, "plot_gallery_colorbar.png"),
    ("workflow_convergence", render_convergence, "plot_gallery_convergence.png"),
    ("workflow_relax", render_relax_convergence, "plot_gallery_relax.png"),
    ("workflow_kpar_ncore", render_kpar_ncore, "plot_gallery_kpar_ncore.png"),
    ("workflow_stretch", render_stretch, "plot_gallery_stretch.png"),
    ("energy", render_energy_components, "plot_gallery_energy.png"),
    ("atominfo_interlayer", render_interlayer_distance, "plot_gallery_interlayer.png"),
    ("atominfo_zpositions", render_zpositions, "plot_gallery_zpositions.png"),
    ("atominfo_rdf", render_rdf, "plot_gallery_rdf.png"),
    ("n2p2_learning", render_learning_curve, "plot_gallery_learning.png"),
    ("n2p2_compare", render_compare, "plot_gallery_compare.png"),
    ("n2p2_rmse_tag", render_rmse_by_tag, "plot_gallery_rmse_tag.png"),
    ("n2p2_epoch_rmse", render_epoch_rmse, "plot_gallery_epoch_rmse.png"),
    ("plotting_pretty", render_pretty_plot, "plot_gallery_pretty.png"),
    ("plotting_periodic", render_periodic_heatmap, "plot_gallery_periodic.png"),
    ("plotting_arkel", render_van_arkel, "plot_gallery_arkel.png"),
    ("render_info", render_render_info, "plot_gallery_render.png"),
    ("ppt_info", render_ppt_info, "plot_gallery_ppt.png"),
    ("dos", render_dos, "plot_gallery_dos.png"),
]


def run_example(path_output: Path) -> list[dict]:
    """Run every section, writing one PNG per entry into ``path_output``."""
    path_output.mkdir(parents=True, exist_ok=True)
    results = []
    for slug, fn, name in SECTION_PLAN:
        path_png = path_output / name
        try:
            fn(path_png)
            ok = True
            err = ""
        except Exception as exc:  # pragma: no cover - report per-section failure
            ok = False
            err = repr(exc)
        results.append({
            "section": slug,
            "image": str(path_png),
            "size_kb": round(path_png.stat().st_size / 1024, 1) if path_png.is_file() else 0.0,
            "ok": ok,
            "error": err,
        })
        print(("wrote: " if ok else "FAILED: ") + str(path_png) + (" " + err if not ok else ""))
    return results


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=Path("docs/_build/example-plot-gallery"),
        help="Output directory for the gallery PNGs.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = run_example(args.output)
    print("\n=== plot gallery summary ===")
    for r in rows:
        print(f"{r['section']:>22}  ok={r['ok']!s:>5}  size_kb={r['size_kb']:>7.1f}  {r['error']}")
    n_ok = sum(1 for r in rows if r["ok"])
    assert n_ok == len(rows), f"only {n_ok}/{len(rows)} sections succeeded"
    print("OK: all assertions passed.")


if __name__ == "__main__":
    main()
