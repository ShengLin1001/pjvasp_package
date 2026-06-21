"""Plots for n2p2 training post-processing.

Built on :func:`mymetal.universal.plot.plot.my_plot` so the figures share the
package's house style (size, fonts, ticks, grid). Three plots:

    - plot_learning_curve : energy/force RMSE vs epoch (log-y)
    - plot_compare        : DFT-vs-NNP scatter (full range + zoom), per tag
    - plot_rmse_by_tag    : per-tag energy/force RMSE bar chart (log-y)

Data come from :mod:`mymetal.ml.n2p2.calculate.post`.
"""

import matplotlib
matplotlib.use('Agg')  # 集群无显示，强制非交互后端
import matplotlib.pyplot as plt
import numpy as np

from mymetal.universal.plot.plot import my_plot
from mymetal.universal.plot.general import general_modify_legend

__all__ = ['plot_learning_curve', 'plot_compare', 'plot_rmse_by_tag']


def _tag_colors(ltag: list) -> dict:
    """Stable tag -> color map (tab20), shared by scatter and bar plots."""
    cmap = plt.get_cmap('tab20')
    return {t: cmap(i % 20) for i, t in enumerate(ltag)}


def _format_plain(value: float, digits: int = 2) -> str:
    """Short decimal label without scientific notation."""
    if not np.isfinite(value):
        return ''
    if value == 0:
        return '0'
    decimals = max(0, min(4, digits - 1 - int(np.floor(np.log10(abs(value))))))
    return f'{value:.{decimals}f}'.rstrip('0').rstrip('.')


def my_plot_learning_curve(epochs: np.ndarray = None, e_rmse_mev: np.ndarray = None,
                        f_rmse_mev: np.ndarray = None, file_path: str = None,
                        label: str = 'Train') -> tuple:
    """Plot energy/force training RMSE versus epoch (log-y).

    Args:
        epochs: Epoch numbers.
        e_rmse_mev: Energy RMSE in meV/atom.
        f_rmse_mev: Force RMSE in meV/Angstrom.
        file_path: Output figure path (e.g. learning_curve.pdf).
        label: Legend label for the curve (e.g. 'Train' or 'Validation').
    """
    fig, ax = my_plot(fig_subp=[1, 2], fig_sharex=True)
    ax[0].semilogy(epochs, e_rmse_mev, '-', label=label.capitalize())
    ax[0].set_xlabel('Epoch (-)')
    ax[0].set_ylabel('Energy RMSE (meV/atom)')
    ax[1].semilogy(epochs, f_rmse_mev, '-', label=label.capitalize())
    ax[1].set_xlabel('Epoch (-)')
    ax[1].set_ylabel(r'Force RMSE (meV/$\mathrm{\AA}$)')
    for a in ax:
        general_modify_legend(a.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95)))
    fig.savefig(file_path)

    return fig, ax


def my_plot_compare(e_ref: np.ndarray = None, e_nnp: np.ndarray = None, tag_e: np.ndarray = None,
                 f_ref: np.ndarray = None, f_nnp: np.ndarray = None, tag_f: np.ndarray = None,
                 file_path: str = None, text_e: str = '', text_f: str = '') -> tuple:
    """DFT-vs-NNP scatter for energies and forces, colored by tag.

    Top row spans the full range; bottom row is clipped to DFT-value
    percentiles so the bulk of the data is not squashed by extreme structures
    (e.g. compressed dimers).

    Args:
        e_ref: DFT energies per atom (eV/atom).
        e_nnp: NNP energies per atom (eV/atom).
        tag_e: Tag label of each energy point.
        f_ref: DFT force components (eV/Angstrom).
        f_nnp: NNP force components (eV/Angstrom).
        tag_f: Tag label of each force component.
        file_path: Output figure path (e.g. dncompare.pdf).
        text_e: Annotation text for the energy panel (e.g. RMSE summary).
        text_f: Annotation text for the force panel.
    """
    tag_e, tag_f = np.asarray(tag_e), np.asarray(tag_f)
    ltag = sorted(set(tag_e))
    color = _tag_colors(ltag)

    fig, ax = my_plot(fig_subp=[2, 2], fig_sharex=False)
    panels = [(ax[0, 0], e_ref, e_nnp, tag_e, None,         'Energies',           r'DFT energy (eV/atom)',          r'NNP energy (eV/atom)'),
              (ax[0, 1], f_ref, f_nnp, tag_f, None,         'Force components',   r'DFT force (eV/$\mathrm{\AA}$)', r'NNP force (eV/$\mathrm{\AA}$)'),
              (ax[1, 0], e_ref, e_nnp, tag_e, (0.0, 98.0),  'Energies (zoom)',    r'DFT energy (eV/atom)',          r'NNP energy (eV/atom)'),
              (ax[1, 1], f_ref, f_nnp, tag_f, (0.5, 99.5),  'Force (zoom)',       r'DFT force (eV/$\mathrm{\AA}$)', r'NNP force (eV/$\mathrm{\AA}$)')]
    for a, ref, nnp, tags, q, title, xl, yl in panels:
        ref, nnp = np.asarray(ref), np.asarray(nnp)
        if q is None:
            lo, hi = min(ref.min(), nnp.min()), max(ref.max(), nnp.max())
        else:
            # 表示裁掉最小 0.5% 和最大 0.5% 的 DFT force 值
            lo, hi = np.percentile(ref, q)
        pad = 0.05 * (hi - lo)
        a.plot([lo - pad, hi + pad], [lo - pad, hi + pad], 'k-', lw=1.0, zorder=1)
        for t in ltag:
            mask = tags == t
            if mask.any():
                # rasterized: 几万个点的 pdf 不会爆体积, 会把散点栅格化，但保留坐标轴、文字、线条为矢量
                a.scatter(ref[mask], nnp[mask], s=8, alpha=0.5, color=color[t],
                          label=t, edgecolors='none', rasterized=True, zorder=2)
        a.set_xlim(lo - pad, hi + pad)
        a.set_ylim(lo - pad, hi + pad)
        a.set_aspect('equal')
        a.set_title(title)
        a.set_xlabel(xl)
        a.set_ylabel(yl)

    leg = ax[0, 0].legend(
        fontsize=16,
        ncol=2,
        loc='center left',
        bbox_to_anchor=(1.05, 0.5),
        borderpad=0.0,
        columnspacing=0.0,
        handletextpad=0.0
    )

    for h in leg.legend_handles:
        h.set_sizes([80])   # 只改变 legend 中 scatter marker 的面积

    general_modify_legend( leg )

    if text_e:
        ax[0, 0].text(0.04, 0.96, text_e, transform=ax[0, 0].transAxes, va='top', fontsize=22)
    if text_f:
        ax[0, 1].text(0.04, 0.96, text_f, transform=ax[0, 1].transAxes, va='top', fontsize=22)
    fig.savefig(file_path)
    return fig, ax


def my_plot_rmse_by_tag(df=None, file_path: str = None) -> tuple:
    """Bar chart of per-tag energy/force RMSE (log-y), value labels on bars.

    Args:
        df: DataFrame from ``build_rmse_by_tag_df`` (the TOTAL row is shown
            separately as a reference line, not as a bar).
        file_path: Output figure path (e.g. rmse_by_tag.pdf).
    """
    d = df[df['tag'] != 'TOTAL'].reset_index(drop=True)
    total = df[df['tag'] == 'TOTAL']
    tags = list(d['tag'])
    x = np.arange(len(tags))
    color = [_tag_colors(tags)[t] for t in tags]

    fig, ax = my_plot(fig_subp=[2, 1], fig_sharex=True, if_keep_wspace_hspace=True,
                      axes_width = 7.31 * 1.5, grid = False, hspace=2.315 * 1.5)
    # my_plot 依据 left/top 边距推出默认 hspace；这里只在返回的 fig 上压缩上下两个共享 x 子图的间距，
    # 不改动 my_plot 本身。值越小两栏越紧凑（默认约 0.39），按需调整。
    fig.subplots_adjust(hspace=0.2)
    specs = [(ax[0], 'E_RMSE_meV/at', 'RMSE of energy (meV/atom)'),
             (ax[1], 'F_RMSE_meV/A', r'RMSE of force (meV/$\mathrm{\AA}$)')]
    for a, col, yl in specs:
        vals = d[col].values
        a.bar(x, vals, color=color, edgecolor=None, zorder=2)
        a.set_yscale('log')
        a.set_ylabel(yl)
        if len(total) and np.isfinite(total[col].values[0]):
            tot = total[col].values[0]
            a.axhline(tot, ls='--', color='gray', zorder=1,
                      label=f'Total = {_format_plain(tot, 3)}')
            general_modify_legend(a.legend(loc='lower left', bbox_to_anchor=(0.00, 1.01)))
        for xi, v in zip(x, vals):
            if np.isfinite(v):
                a.text(xi, v, _format_plain(v, 2), ha='center', va='bottom', rotation=90)
        y_low, _ = a.get_ylim()
        a.margins(y=0.2)
        _, y_high = a.get_ylim()
        a.set_ylim(bottom=y_low, top=y_high)
    ax[0].set_xlabel('')  # my_plot/general_axes 默认会写 'X-axis Label'，共享 x 的上排清掉
    ax[1].set_xticks(x)
    ax[1].set_xticklabels(tags, rotation=45)
    ax[1].set_xlabel('Tag')
    fig.savefig(file_path)
    return fig, ax
