import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mymetal.universal.plot.general import general_font, general_set_all_rcParams
from mymetal.universal.plot.plot import my_plot


def test_fontsize_scales_unspecified_axis_and_tick_labels():
    with plt.rc_context():
        general_set_all_rcParams(backend="Agg", fontsize=12)

        assert plt.rcParams["axes.labelsize"] == 12
        assert plt.rcParams["axes.linewidth"] < 1.3
        assert plt.rcParams["grid.linewidth"] < 0.9
        assert plt.rcParams["lines.linewidth"] < 1.3
        assert plt.rcParams["xtick.labelsize"] == 12
        assert plt.rcParams["ytick.labelsize"] == 12


def test_general_font_resets_axis_and_tick_labels():
    with plt.rc_context():
        general_set_all_rcParams(backend="Agg", fontsize=10)
        general_font(fontsize=28)

        assert plt.rcParams["axes.labelsize"] == 28
        assert plt.rcParams["xtick.labelsize"] == 28
        assert plt.rcParams["ytick.labelsize"] == 28


def test_my_plot_accepts_compact_style_without_changing_defaults():
    with plt.rc_context():
        fig, ax = my_plot(
            fontsize=18,
            legend_fontsize=16,
            markersize=9,
            linewidth=1.5,
            markeredgewidth=1.0,
            tick_width=1.2,
        )

        assert plt.rcParams["font.size"] == 18
        assert plt.rcParams["lines.linewidth"] == 1.5
        assert ax.spines["left"].get_linewidth() == 1.5
        assert ax.xaxis.majorTicks[0].tick1line.get_markeredgewidth() == 1.2
        plt.close(fig)
