"""Plots for n2p2 training post-processing.

Built on :func:`mymetal.universal.plot.plot.my_plot` so the figures share the
package's house style (size, fonts, ticks, grid).

Per-training-run RMSE / scatter plots (data from
:mod:`mymetal.ml.n2p2.calculate.post`):

    - my_plot_learning_curve : energy/force RMSE vs epoch (log-y)
    - my_plot_compare        : DFT-vs-NNP scatter (full range + zoom), per tag
    - my_plot_rmse_by_tag    : per-tag energy/force RMSE bar chart (log-y)

LAMMPS (pair_style hdnnp) physical properties vs training epoch (data
aggregated from the per-epoch mymetal post readers ``my_read_stretch`` /
``read_cij_energy`` / ``read_output``):

    - my_plot_epoch_stretch  : a / c (c/a for HCP) / E vs epoch, 3 rows x phase cols
    - my_plot_epoch_cij      : C11/C12/C13/C33/C44 vs epoch, 5 rows x phase cols
    - my_plot_epoch_gsfe     : stable/unstable SFE vs epoch, 2 rows x slip-system cols

All three epoch grids share the x-axis (epoch); each panel is one
default-coloured (C0) line with a C1 circle marker every 50 epochs, labelled by a
column title (phase / slip system) and its own y-axis label instead of a
legend. y-limits track the converged tail (last 75% of the data, 5% margin);
grid is off; axes keep the my_plot preset size with a tightened inter-row gap
(if_keep_wspace_hspace), and figures are saved with bbox_inches='tight'.
"""

import matplotlib
matplotlib.use('Agg')  # 集群无显示，强制非交互后端
import matplotlib.pyplot as plt
import numpy as np

from mymetal.universal.plot.plot import my_plot
from mymetal.universal.plot.general import general_modify_legend

__all__ = ['my_plot_learning_curve', 'my_plot_compare', 'my_plot_rmse_by_tag',
           'my_plot_epoch_stretch', 'my_plot_epoch_cij', 'my_plot_epoch_gsfe']


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
    plt.close(fig)

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
    plt.close(fig)
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
    plt.close(fig)
    return fig, ax


# 共享 x 的 epoch 网格统一参数：
#   _EPOCH_HSPACE   行间留白（英寸）。配合 my_plot(if_keep_wspace_hspace=True) 让每个 axes 保持
#                   预设大小（7.31 x 5.89），只缩小竖直间距而不是把 axes 撑大；wspace 不传、保持默认。
#   _EPOCH_MARK_EVERY  每隔多少 epoch 打一个 marker。
_EPOCH_HSPACE = 1.0
_EPOCH_MARK_EVERY = 50


def _clear_default_labels(ax) -> None:
    """Drop the placeholder 'X-/Y-axis Label' that general_axes writes on every axis."""
    for a in np.atleast_2d(ax).ravel():
        a.set_xlabel('')
        a.set_ylabel('')


def _dft_ref_val(dft: dict, key: str):
    """Return the finite DFT reference for ``key``, or ``None`` when unavailable.

    ``dft`` is one of the flat sub-dicts returned by
    :func:`mymetal.ml.n2p2.dataset.read_dft_reference` (keys identical to the
    epoch-scan columns, e.g. ``a_fcc`` / ``C11_hcp`` / ``usf_FCC_111``). Returns
    ``None`` when ``dft`` is falsy, lacks the key, or the value is NaN (that
    phase / slip system is simply absent from the DFT archive). Shared by
    :func:`_epoch_panel` (to widen the y-limit) and :func:`_epoch_dft_ref` (to
    draw the reference line), so both agree on what counts as a valid reference.
    """
    if not dft:
        return None
    val = dft.get(key)
    if val is None or not np.isfinite(val):
        return None
    return float(val)


def _epoch_panel(a, x, y, mark_every: int = _EPOCH_MARK_EVERY,
                 dft: dict = None, key: str = None) -> None:
    """Draw one epoch series with house defaults and a data-driven y-limit.

    Full line in the default style (C0, default linewidth) with a circle marker
    every ``mark_every`` epochs placed via ``markevery`` on the same artist. The
    markers are coloured C1 so they stand out against the C0 line;
    marker size keeps its default. The y-limits span the last 75% of the (finite)
    data with a 5% margin top and bottom, so the early, still-converging epochs
    don't squash the converged tail. When a finite DFT reference (``dft[key]``)
    lies outside that data range it is folded in as an extra bound, so the grey
    dashed reference line drawn by :func:`_epoch_dft_ref` stays inside the axes
    instead of falling off it. An all-NaN panel (e.g. a slip system without a
    stable SFE) is left empty.
    """
    x = np.asarray(x)
    y = np.asarray(y, dtype=float)
    if not np.isfinite(y).any():
        return
    mark_idx = np.nonzero(x % mark_every == 0)[0]          # marker 落在 epoch 的整 mark_every 倍处
    a.plot(x, y, marker='o', markevery=mark_idx.tolist(),
           markeredgecolor='C1')     # 线默认 C0，marker 统一用 C1

    tail = y[int(round(0.25 * len(y))):]                   # 后 75% 数据定 ylim
    tail = tail[np.isfinite(tail)]
    if tail.size:
        lo, hi = float(tail.min()), float(tail.max())
        ref = _dft_ref_val(dft, key)                       # DFT 参考也算作数据边界，避免虚线跑出坐标轴
        if ref is not None:
            lo, hi = min(lo, ref), max(hi, ref)
        span = hi - lo
        pad = 0.05 * span if span > 0 else (0.05 * abs(hi) if hi else 0.05)
        a.set_ylim(lo - pad, hi + pad)


def _epoch_dft_ref(a, dft: dict, key: str) -> None:
    """Overlay the DFT reference for ``key`` as a grey dashed line on panel ``a``.

    No-op when :func:`_dft_ref_val` reports no finite reference for ``key``.
    Linewidth stays at the house default so the line reads as a flat reference,
    not another data series.
    """
    val = _dft_ref_val(dft, key)
    if val is None:
        return
    a.axhline(val, ls='--', color='gray', zorder=1)


def _pretty_type(t: str) -> str:
    """'HCP_basal' -> 'HCP Basal' (capitalise each token, keep all-caps tokens like FCC/HCP)."""
    return ' '.join(p if p.isupper() else p.capitalize() for p in t.split('_'))


def my_plot_epoch_stretch(df=None, file_path: str = None,
                          phases=('hcp', 'fcc', 'bcc'), dft=None) -> tuple:
    """Equilibrium stretch quantities versus training epoch (3 rows x phase cols).

    One column per phase (titles ``HCP|FCC|BCC``), shared x (epoch). Rows top to
    bottom: lattice constant ``a``, ``c`` (the HCP column shows ``c/a`` instead),
    and energy per atom ``E``. One default-coloured (C0) line per panel (circle
    marker every 50 epochs); phases are told apart by the column + its own y
    label, not by colour, so there is no legend.

    Args:
        df: DataFrame with column ``epoch`` and (per phase) ``a_<p>``, ``c_<p>``,
            ``E_<p>`` plus ``ca_hcp`` (built in ``PeiN2p2.post_epoch_scan``).
        file_path: Output figure path (e.g. p_post_epoch_stretch.pdf).
        phases: Column order (default hcp/fcc/bcc).
        dft: Optional ``read_dft_reference(...)['stretch']`` sub-dict; each
            matching quantity is drawn as a grey dashed horizontal reference line.
    """
    d = df.copy()
    ep = d['epoch'].to_numpy()

    fig, ax = my_plot(fig_subp=[3, len(phases)], fig_sharex=True, grid=False,
                      if_keep_wspace_hspace=True, hspace=_EPOCH_HSPACE)
    ax = np.atleast_2d(ax)
    _clear_default_labels(ax)

    for j, p in enumerate(phases):
        if f'a_{p}' in d:                                  # row 0: a
            _epoch_panel(ax[0, j], ep, d[f'a_{p}'].to_numpy(), dft=dft, key=f'a_{p}')
        _epoch_dft_ref(ax[0, j], dft, f'a_{p}')            # 叠 DFT 灰色虚线参考
        ax[0, j].set_ylabel(r'$a$ ($\mathrm{\AA}$)')
        if p == 'hcp':                                     # row 1: HCP 列画 c/a
            if 'ca_hcp' in d:
                _epoch_panel(ax[1, j], ep, d['ca_hcp'].to_numpy(), dft=dft, key='ca_hcp')
            _epoch_dft_ref(ax[1, j], dft, 'ca_hcp')
            ax[1, j].set_ylabel(r'$c/a$ (-)')
        else:                                              # row 1: 立方相画 c
            if f'c_{p}' in d:
                _epoch_panel(ax[1, j], ep, d[f'c_{p}'].to_numpy(), dft=dft, key=f'c_{p}')
            _epoch_dft_ref(ax[1, j], dft, f'c_{p}')
            ax[1, j].set_ylabel(r'$c$ ($\mathrm{\AA}$)')
        if f'E_{p}' in d:                                  # row 2: E
            _epoch_panel(ax[2, j], ep, d[f'E_{p}'].to_numpy(), dft=dft, key=f'E_{p}')
        _epoch_dft_ref(ax[2, j], dft, f'E_{p}')
        ax[2, j].set_ylabel('Energy (eV/atom)')
        ax[0, j].set_title(p.upper())                      # 列标题：相
        ax[-1, j].set_xlabel('Epoch (-)')                  # 仅底排写 x 轴

    fig.savefig(file_path, bbox_inches='tight')
    plt.close(fig)
    return fig, ax


def my_plot_epoch_cij(df=None, file_path: str = None,
                      phases=('hcp', 'fcc', 'bcc'),
                      comps=('11', '12', '13', '33', '44'), dft=None) -> tuple:
    """Elastic constants Cij versus training epoch (5 rows x phase cols).

    One column per phase (titles ``HCP|FCC|BCC``), shared x (epoch); rows top to
    bottom are C11, C12, C13, C33, C44. One default-coloured (C0) line per panel
    (circle marker every 50 epochs); no legend (the component is the row + y
    label, the phase is the column).

    Args:
        df: DataFrame with column ``epoch`` and ``C<ij>_<p>`` columns
            (built in ``PeiN2p2.post_epoch_scan``).
        file_path: Output figure path (e.g. p_post_epoch_cij.pdf).
        phases: Column order (default hcp/fcc/bcc).
        comps: Cij components, one row each (default 11/12/13/33/44).
        dft: Optional ``read_dft_reference(...)['cij']`` sub-dict; each matching
            Cij is drawn as a grey dashed horizontal reference line.
    """
    d = df.copy()
    ep = d['epoch'].to_numpy()

    fig, ax = my_plot(fig_subp=[len(comps), len(phases)], fig_sharex=True, grid=False,
                      if_keep_wspace_hspace=True, hspace=_EPOCH_HSPACE)
    ax = np.atleast_2d(ax)
    _clear_default_labels(ax)

    for i, c in enumerate(comps):
        for j, p in enumerate(phases):
            col = f'C{c}_{p}'
            if col in d:
                _epoch_panel(ax[i, j], ep, d[col].to_numpy(), dft=dft, key=col)
            _epoch_dft_ref(ax[i, j], dft, col)             # 叠 DFT 灰色虚线参考
            ax[i, j].set_ylabel(rf'$C_{{{c}}}$ (GPa)')     # 每个面板都标 y
            if i == 0:
                ax[i, j].set_title(p.upper())              # 列标题：相
    for j in range(ax.shape[1]):
        ax[-1, j].set_xlabel('Epoch (-)')                  # 仅底排写 x 轴

    fig.savefig(file_path, bbox_inches='tight')
    plt.close(fig)
    return fig, ax


def my_plot_epoch_gsfe(df=None, file_path: str = None, types=None, dft=None) -> tuple:
    """Stacking-fault energies versus training epoch (2 rows x slip-system cols).

    One column per slip system (title = system name), shared x (epoch); top row
    is the stable SFE ``gamma_sf`` (last gamma), bottom row the unstable SFE
    ``gamma_usf`` (max gamma). One default-coloured (C0) line per panel (circle
    marker every 50 epochs); no legend (the system is the column + y label). A
    system with a valid GSFE table has both panels populated.

    Args:
        df: DataFrame with column ``epoch`` and ``usf_<type>`` / ``sf_<type>``
            columns (built in ``PeiN2p2.post_epoch_scan``).
        file_path: Output figure path (e.g. p_post_epoch_gsfe.pdf).
        types: Column order of slip systems. Defaults to the Au FCC/HCP set
            (FCC_111, FCC_100, then the four HCP systems).
        dft: Optional ``read_dft_reference(...)['gsfe']`` sub-dict; each matching
            usf/sf is drawn as a grey dashed horizontal reference line.
    """
    if types is None:
        types = ['FCC_111', 'FCC_100', 'HCP_basal', 'HCP_prism1w', 'HCP_pyr1w', 'HCP_pyr2']

    d = df.copy()
    ep = d['epoch'].to_numpy()

    fig, ax = my_plot(fig_subp=[2, len(types)], fig_sharex=True, grid=False,
                      if_keep_wspace_hspace=True, hspace=_EPOCH_HSPACE)
    ax = np.atleast_2d(ax)
    _clear_default_labels(ax)

    # 行：0 = 稳定层错能 sf，1 = 不稳定层错能 usf（与 post_epoch_scan 的列名一致）
    ylabels = {'sf': r'$\gamma_{\mathrm{sf}}$ (mJ/m$^2$)',
               'usf': r'$\gamma_{\mathrm{usf}}$ (mJ/m$^2$)'}
    for j, t in enumerate(types):
        for i, pref in enumerate(('sf', 'usf')):
            col = f'{pref}_{t}'
            if col in d:
                _epoch_panel(ax[i, j], ep, d[col].to_numpy(), dft=dft, key=col)
            _epoch_dft_ref(ax[i, j], dft, col)             # 叠 DFT 灰色虚线参考
            ax[i, j].set_ylabel(ylabels[pref])             # 每个面板都标 y
        ax[0, j].set_title(_pretty_type(t))                # 列标题：滑移系（首字母大写）
        ax[-1, j].set_xlabel('Epoch (-)')                  # 仅底排写 x 轴

    fig.savefig(file_path, bbox_inches='tight')
    plt.close(fig)
    return fig, ax
