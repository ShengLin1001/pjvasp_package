import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from typing import List, Tuple, Union
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
from brokenaxes import brokenaxes
from matplotlib.ticker import MaxNLocator

def my_plot(
    one_fig_wh: List[float] = [10.72, 8.205],    
    fig_subp: List[int] = [1, 1], 
    fig_sharex: bool = True, 
    grid: bool = True, 
    labelpad: int = 15, 
    tick_pad: int = 10, 
    left: float = 1.918, 
    top: float = 0.9517,
    axes_height: float = 5.89,
    axes_width: float = 7.31
) -> Tuple[Figure, List[Axes]]:
    """Creates a customized matplotlib figure with specific layout adjustments.

    This function sets global matplotlib settings, creates a figure with 
    subplots, and adjusts their layout according to the given parameters. 
    It also enables minor ticks and customizes the axes' appearance, margins, 
    and tick settings.

    Args:
        one_fig_wh (List[float]): Width and height of the figure in inches 
            (default: [10.72, 8.205]).
        fig_subp (List[int]): Number of rows and columns of subplots 
            (default: [1, 1]).
        fig_sharex (bool): Whether the subplots share the same x-axis 
            (default: True).
        grid (bool): Whether to display a grid on the plots (default: True).
        labelpad (int): Padding between axis labels and the axis (default: 15).
        tick_pad (int): Padding between tick labels and the axes (default: 10).
        left (float): Left margin in inches (default: 1.918).
        top (float): Top margin in inches (default: 0.9517).
        axes_height (float): Height of each subplot in inches (default: 5.89).
        axes_width (float): Width of each subplot in inches (default: 7.31).

    Returns:
        Tuple[Figure, List[Axes]]: A tuple containing:
            - fig (Figure): The created figure object.
            - axes (List[Axes]): A list of axes objects for the subplots. 
              If there's only one subplot, it returns a list with one item.
    """
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 28
    plt.rcParams['axes.linewidth'] = 3
    plt.rcParams['axes.grid'] = grid
    plt.rcParams['grid.linestyle'] = '--'
    plt.rcParams['grid.linewidth'] = 0.2
    plt.rcParams["savefig.transparent"] = 'True'
    plt.rcParams['lines.linewidth'] = 3
    plt.rcParams['lines.markersize'] = 20
    plt.rcParams['lines.markeredgewidth'] = 3
    plt.rcParams['lines.markerfacecolor'] = 'white'

    # 图例相关全局参数
    plt.rcParams['legend.loc'] = 'upper right'
    plt.rcParams['legend.fontsize'] = 24
    plt.rcParams['legend.frameon'] = True
    plt.rcParams['legend.borderpad'] = 0.0
    plt.rcParams['legend.labelspacing'] = 0.5
    plt.rcParams['legend.columnspacing'] = 0.5

    fig_wh = [a * b for a, b in zip(one_fig_wh, fig_subp)]
    fig, ax = plt.subplots(nrows=fig_subp[0], ncols=fig_subp[1], 
                           sharex=fig_sharex, figsize=(fig_wh[0], fig_wh[1]))
    # inch 10.72* 8.205 is one figsize
    left = left
    right = one_fig_wh[0]-left-axes_width
    top = top
    bottom = one_fig_wh[1]-top-axes_height
    wspace = (left+right)/axes_width  # axes_width 7.31 inch
    hspace = (top+bottom)/axes_height  # axes_height 5.89 inch

    fig.subplots_adjust(left=left/fig_wh[0], right=(fig_wh[0]-right)/fig_wh[0], 
                        top=(fig_wh[1]-top)/fig_wh[1], bottom=(bottom)/fig_wh[1],
                        hspace= (hspace), wspace=(wspace))


    if isinstance(ax, np.ndarray):  # 多子图情况
        axes = ax.flatten()
    else:  # 单子图情况
        axes = [ax]

    for axis in axes:
        # 启用副刻度
        axis.minorticks_on()
        axis.xaxis.set_minor_locator(AutoMinorLocator(2))
        axis.yaxis.set_minor_locator(AutoMinorLocator(2))
        axis.tick_params(which='major', direction='in', length=8, width=3.0, pad = tick_pad)
        axis.tick_params(which='minor', direction='in', length=4, width=3.0, pad = tick_pad)

        # 设置轴标签，并调整与刻度的距离
        axis.set_xlabel('X-axis Label', labelpad=labelpad)
        axis.set_ylabel('Y-axis Label', labelpad=labelpad)

    def general_modify_legend(legend):
        legend.get_frame().set_boxstyle("Square", pad=0.5)
        legend.get_frame().set_linewidth(2.5)
        legend.get_frame().set_edgecolor('black')
        legend.get_frame().set_facecolor('white')
        legend.get_frame().set_alpha(1)
        
    def general_margin_bin(axis,
                            y_margin = 0.1, 
                            y_nbins = 3,
                            x_margin = 0.1, 
                            x_nbins = 4, 
                            prune = 'both',
                            if_autoscale = True):
        axis.margins(x=x_margin, y=y_margin)
        axis.yaxis.set_major_locator(MaxNLocator(nbins=y_nbins, prune=prune))
        axis.xaxis.set_major_locator(MaxNLocator(nbins=x_nbins, prune=prune))
        if if_autoscale:
            # Force axis limits to update based on the new margins
            axis.autoscale(enable=True, axis='both', tight=False)
    
    # 将修改图例函数绑定到 fig 对象，方便调用
    fig.general_modify_legend = general_modify_legend
    fig.general_margin_bin = general_margin_bin

    return fig, ax

def my_plot_brokenaxed(
    one_fig_wh: List[float] = [10.72, 8.205],    
    fig_subp: List[int] = [1, 1], 
    fig_sharex: bool = True, 
    grid: bool = True, 
    tick_pad: int = 10, 
    left: float = 1.918, 
    top: float = 0.9517, 
    ylims: List[Tuple[float, float]] = [(0, 1), (2, 3), (4, 5)],
    axes_height: float = 5.89,
    axes_width: float = 7.31,
) -> Tuple[Figure, brokenaxes]:
    """Creates a broken axes plot with customized layout and legend settings.

    This function creates a figure with subplots using brokenaxes, adjusts 
    the layout, and sets up the appearance of ticks, grid, and legend styles.

    Args:
        one_fig_wh (List[float]): Width and height of the figure in inches (one figure)
            (default: [10.72, 8.205]).
        fig_subp (List[int]): Number of rows and columns of subplots 
            (default: [1, 1]).
        fig_sharex (bool): Whether the subplots share the same x-axis 
            (default: True).
        grid (bool): Whether to display a grid on the plots (default: True).
        tick_pad (int): Padding between tick labels and the axes (default: 10).
        left (float): Left margin in inches (default: 1.918).
        top (float): Top margin in inches (default: 0.9517).
        ylims (List[Tuple[float, float]]): List of y-axis limits for each segment 
            of the broken axes (default: [(0, 1), (2, 3), (4, 5)]).

    Returns:
        Tuple[Figure, brokenaxes]: A tuple containing:
            - fig (Figure): The created figure object.
            - ax (brokenaxes): The brokenaxes object containing the subplots.
    """

    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 28
    plt.rcParams['axes.linewidth'] = 3
    plt.rcParams['axes.grid'] = grid
    plt.rcParams['grid.linestyle'] = '--'
    plt.rcParams['grid.linewidth'] = 0.2
    plt.rcParams["savefig.transparent"] = 'True'
    plt.rcParams['lines.linewidth'] = 3
    plt.rcParams['lines.markersize'] = 20
    plt.rcParams['lines.markeredgewidth'] = 3
    plt.rcParams['lines.markerfacecolor'] = 'white'

    # 图例相关全局参数
    plt.rcParams['legend.loc'] = 'upper right'
    plt.rcParams['legend.fontsize'] = 24
    plt.rcParams['legend.frameon'] = True
    plt.rcParams['legend.borderpad'] = 0.0
    plt.rcParams['legend.labelspacing'] = 0.5
    plt.rcParams['legend.columnspacing'] = 0.5

    fig_wh = [a * b for a, b in zip(one_fig_wh, fig_subp)]
    # inch 10.72* 8.205 is one figsize
    left = left
    right = one_fig_wh[0]-left-axes_width
    top = top
    bottom = one_fig_wh[1]-top-axes_height
    wspace = (left+right)/axes_width  # axes_width 7.31 inch
    hspace = (top+bottom)/axes_height  # axes_height 5.89 inch

    fig = plt.figure(figsize=(fig_wh[0], fig_wh[1]))
    fig.subplots_adjust(left=left/fig_wh[0], right=(fig_wh[0]-right)/fig_wh[0], 
                    top=(fig_wh[1]-top)/fig_wh[1], bottom=(bottom)/fig_wh[1],
                    hspace= (hspace), wspace=(wspace))
    ax = brokenaxes(ylims=ylims, despine=False, hspace=0.05, d=0.01)

    for axis in ax.axs:
        axis.minorticks_on()
        axis.xaxis.set_minor_locator(AutoMinorLocator(2))
        axis.yaxis.set_minor_locator(AutoMinorLocator(2))
        axis.tick_params(which='major', direction='in', length=8, width=3.0, pad = tick_pad)
        axis.tick_params(which='minor', direction='in', length=4, width=3.0, pad = tick_pad)


    def general_modify_legend(legend):
        legend.get_frame().set_boxstyle("Square", pad=0.5)
        legend.get_frame().set_linewidth(2.5)
        legend.get_frame().set_edgecolor('black')
        legend.get_frame().set_facecolor('white')
        legend.get_frame().set_alpha(1)

    fig.general_modify_legend = general_modify_legend

    return fig, ax
