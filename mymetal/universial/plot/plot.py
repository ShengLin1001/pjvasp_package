"""
plot submodule

This submodule contains functions for creating customized plots.

Functions:
    - my_plot: Creates a customized matplotlib figure with specific layout adjustments.
    - my_plot_brokenaxed: Creates a broken axes plot with customized layout and legend settings.
    - my_plot_energy_components: Generates plots for energy components with polynomial fits and comparisons.
    - my_plot_interlayer_distance: Calculates and plots the interlayer distances from the given atomic positions.
    - my_plot_zpositions: Plot z-positions.
"""

from ase import Atoms
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.colors import to_rgb
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, MaxNLocator
from brokenaxes import brokenaxes
from typing import List, Tuple, Union
import matplotlib.colors as mcolors
from itertools import cycle

import matplotlib.pyplot as plt

def check_font_size(ax: plt.Axes):
    x_label_fontsize = ax.xaxis.get_label().get_fontsize()
    y_label_fontsize = ax.yaxis.get_label().get_fontsize()
    print(f"Font size of x-axis label: {x_label_fontsize}")
    print(f"Font size of y-axis label: {y_label_fontsize}")

    x_tick_labels = ax.get_xticklabels()
    y_tick_labels = ax.get_yticklabels()

    if x_tick_labels:
        x_tick_font = x_tick_labels[0].get_fontproperties()
        print(f"Font name of the first major tick label on the x-axis: {x_tick_font.get_name()}")
        print(f"Font size of the first major tick label on the x-axis: {x_tick_font.get_size()}")

    if y_tick_labels:
        y_tick_font = y_tick_labels[0].get_fontproperties()
        print(f"Font name of the first major tick label on the y-axis: {y_tick_font.get_name()}")
        print(f"Font size of the first major tick label on the y-axis: {y_tick_font.get_size()}")


def get_ploted_figure():
    fig = plt.gcf()
    ax = plt.gca()
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    return fig, ax, xlim, ylim

def general_modify_band_plot(ax: plt.Axes) -> plt.Axes:
    xticks = ax.get_xticks()
    xticklabels = ax.get_xticklabels()
    xlabels = ax.get_xlabel()
    ax.tick_params(axis='x', which='both', bottom=False)
    for xtick in xticks:
        ax.axvline(x=xtick, color='gray', linestyle='--')
    return ax

def generate_gradient_colors(start_color, end_color, num_colors):
    cmap = mcolors.LinearSegmentedColormap.from_list("gradient", [start_color, end_color])
    gradient_colors = [cmap(i / (num_colors - 1)) for i in range(num_colors)]
    return gradient_colors

def general_modify_line(ax: plt.Axes,
                    remove_dashed_line: bool = True,
                    remove_line_style: list = ['--', ':'],
                    color_list: list = ['black', 'red', 'blue', 'green', 'orange', 'purple', 
                                        'brown', 'pink', 'gray', 'olive', 'cyan', 'lime', 
                                        'teal', 'lavender', 'tan', 'salmon', 'gold', 'lightcoral', 
                                        'skyblue', 'darkblue', 'darkred', 'darkgreen', 'darkorange', 
                                        ],
                    if_gradient_color: bool = False,
                    gradient_colors: list=['red', 'blue'],
                    if_cmap_color: bool = False,
                    cmap_color: str = 'coolwarm',
                    if_change_color: bool = False,
                    ):
    all_lines = ax.get_lines()
    all_index = [index for index, line in enumerate(all_lines)]
    removed_index = []
    after_removed_lines = []
    if remove_dashed_line:
        for index, line in enumerate(all_lines):
            if line.get_linestyle() in remove_line_style:
                line.set_visible(False)
                removed_index.append(index)
            else:
                after_removed_lines.append(line)
    if if_gradient_color:
        color_list = generate_gradient_colors(gradient_colors[0], gradient_colors[1], len(after_removed_lines))
    elif if_cmap_color:
        cmap = plt.get_cmap(cmap_color)
        values = np.linspace(0, 1, len(after_removed_lines))
        color_list = cmap(values)
    color_cycle = cycle(color_list)
    if if_change_color:
        for index, line in enumerate(after_removed_lines):
            line.set_color(next(color_cycle))
    return ax


def general_modify_ploted_figure(axes: plt.Axes,
                                 grid: bool = True,
                                 grid_linewidth: float = 0.5,
                                ########################### Figure size
                                one_fig_wh: list = [10.72, 8.205],
                                fig_subp: List[int] = [1, 1],
                                axes_height: float = 5.89,
                                axes_width: float = 7.31,
                                left: float = 1.918,
                                top: float = 0.9517,
                                ########################### axes
                                labelpad: int = 15, 
                                tick_pad: int = 10, 
                                xlabel: str = 'Energies (eV)', 
                                ylabel: str = 'Density of States (a.u.)',
                                xlim: list = [-15, 10],
                                ylim: list = [0, 8],
                                if_show_right_top_tick: bool = False,
                                ########################### legend
                                loc: str = 'upper right',
                                bbox_to_anchor: tuple = (0.95, 0.95),
                                ncol: int = 1,
                                if_show_legend: bool = True,
                                ########################### tick
                                y_margin = 0.1, 
                                y_nbins = 3,
                                x_margin = 0.1, 
                                x_nbins = 4, 
                                prune = 'both',
                                if_autoscale = True,
                                change_margin: bool = False,
                                ########################### line, color
                                remove_dashed_line: bool = True,
                                remove_line_style: list = ['--', ':'],
                                color_list: list = ['black', 'red', 'blue', 'green', 'orange', 'purple', 
                                                    'brown', 'pink', 'gray', 'olive', 'cyan', 'lime', 
                                                    'teal', 'lavender', 'tan', 'salmon', 'gold', 'lightcoral', 
                                                    'skyblue', 'darkblue', 'darkred', 'darkgreen', 'darkorange', 
                                                    ],
                                if_gradient_color: bool = False,
                                gradient_colors: list=['red', 'blue'],
                                if_cmap_color: bool = False,
                                cmap_color: str = 'coolwarm',
                                if_change_color: bool = False,
                                ########################### plot type
                                if_band_plot: bool = False,
                                ########################### save
                                if_save: bool = True,
                                save_path: str = './dos.jpg'
                                )  -> tuple:
    """
    Modifies the appearance of a plotted figure, including axes, legend, grid, margins, and saves the plot.

    Notes:
        Please before you draw the figure, you should use the function `general_font(grid, gridwidth)` to predefine some value.
        Because the change after drawing can't influence the figure.
    
    Args:
        axes (plt.Axes): The axes to modify. Now it don't used, so you don't need to assign it.
        grid (bool): Whether to display grid. Default is True.
        grid_linewidth (float): Line width for the grid. Default is 0.5.
        one_fig_wh (list): The figure's width and height in inches. Default is [10.72, 8.205].
        fig_subp (List[int]): Number of subplots in rows and columns. Default is [1, 1].
        axes_height (float): Height of the axes in inches. Default is 5.89.
        axes_width (float): Width of the axes in inches. Default is 7.31.
        left (float): Left margin in inches. Default is 1.918.
        top (float): Top margin in inches. Default is 0.9517.
        labelpad (int): Padding for axis labels. Default is 15.
        tick_pad (int): Padding for axis ticks. Default is 10.
        xlabel (str): Label for the x-axis. Default is 'Energies (eV)'.
        ylabel (str): Label for the y-axis. Default is 'Density of States (a.u.)'.
        xlim (list): x-axis limits. Default is [-15, 10].
        ylim (list): y-axis limits. Default is [0, 8].
        loc (str): Location for the legend. Default is 'upper right'.
        bbox_to_anchor (tuple): Bounding box for legend. Default is (0.95, 0.95).
        y_margin (float): y-axis margin. Default is 0.1.
        y_nbins (int): Number of bins for y-axis. Default is 3.
        x_margin (float): x-axis margin. Default is 0.1.
        x_nbins (int): Number of bins for x-axis. Default is 4.
        prune (str): Whether to prune axis ticks. Default is 'both'.
        if_autoscale (bool): Whether to autoscale the axes. Default is True.
        remove_dashed_line (bool): Whether to remove dashed lines. Default is True.
        change_margin (bool): Whether to change the axis margins. Default is False.
        if_save (bool): Whether to save the figure. Default is True.
        save_path (str): Path to save the figure. Default is './dos.jpg'.

    Returns:
        tuple: The modified figure and axes.
    """

    general_font(grid, grid_linewidth, if_ploted=True)
    ax = axes
    # pyplot
    fig, ax, xlim0, ylim0 = get_ploted_figure()
    print('The type(ax) is: ',type(ax))
    fig_wh = fig_subp[1]*one_fig_wh[0], fig_subp[0]*one_fig_wh[1]
    fig.set_size_inches(fig_wh)
    general_subplots_adjust(fig, one_fig_wh, fig_subp, axes_height, axes_width, left, top)
    if xlim == None or ylim == None:
        xlim = xlim0
        ylim = ylim0
    general_axes(ax, labelpad, tick_pad,xlabel, ylabel, xlim, ylim, True, if_show_right_top_tick)
    if change_margin:
        general_margin_bin(ax, y_margin, y_nbins, x_margin, x_nbins, prune, if_autoscale)
    general_modify_line(ax, remove_dashed_line, remove_line_style, color_list, if_gradient_color, gradient_colors, if_cmap_color, cmap_color, if_change_color)
    if if_band_plot:
        general_modify_band_plot(ax)
    ax.legend().remove()
    if if_show_legend:
        general_modify_legend(ax.legend(loc=loc, bbox_to_anchor=bbox_to_anchor, ncol=ncol))

    if if_save:
        plt.savefig(save_path, dpi=300)
    return fig, ax

def general_modify_legend(legend):
    """
    Customizes the appearance of a legend.

    Args:
        legend: The legend object to modify.

    Returns:
        None
    """
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
    """
    Adjusts the margins and binning for the x and y axes.

    Args:
        axis: The axis object to modify.
        y_margin (float): The y-axis margin. Default is 0.1.
        y_nbins (int): Number of bins for the y-axis. Default is 3.
        x_margin (float): The x-axis margin. Default is 0.1.
        x_nbins (int): Number of bins for the x-axis. Default is 4.
        prune (str): Prune setting for axis ticks. Default is 'both'.
        if_autoscale (bool): Whether to autoscale the axis. Default is True.

    Returns:
        None
    """
    axis.margins(x=x_margin, y=y_margin)
    axis.yaxis.set_major_locator(MaxNLocator(nbins=y_nbins, prune=prune))
    axis.xaxis.set_major_locator(MaxNLocator(nbins=x_nbins, prune=prune))
    if if_autoscale:
        # Force axis limits to update based on the new margins
        axis.autoscale(enable=True, axis='both', tight=False)

def general_font( grid: bool = True, grid_linewidth: float = 0.5, if_ploted: bool = False,
                  fontsize: int = 28, markersize: int = 20, linewidth: int = 3, 
                  markeredgewidth: int = 3, legend_fontsize: int = 24, 
                  markerfacecolor: str = 'white'):
    """
    Sets global font and plot settings for consistency across figures.

    Args:
        grid (bool): Whether to show grid lines. Default is True.
        grid_linewidth (float): Line width for the grid. Default is 0.5.

    Returns:
        None
    """
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['axes.linewidth'] = linewidth
    plt.rcParams['axes.grid'] = grid
    plt.rcParams['grid.linestyle'] = '--'
    plt.rcParams['grid.linewidth'] = grid_linewidth
    plt.rcParams["savefig.transparent"] = 'True'
    plt.rcParams['lines.linewidth'] = linewidth
    plt.rcParams['lines.markersize'] = markersize
    plt.rcParams['lines.markeredgewidth'] = markeredgewidth
    plt.rcParams['lines.markerfacecolor'] = markerfacecolor

    # 图例相关全局参数
    plt.rcParams['legend.loc'] = 'upper right'
    plt.rcParams['legend.fontsize'] = legend_fontsize
    plt.rcParams['legend.frameon'] = True
    plt.rcParams['legend.borderpad'] = 0.0
    plt.rcParams['legend.labelspacing'] = 0.5
    plt.rcParams['legend.columnspacing'] = 0.5

    # 设置数学字体
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = 'Arial'  # 设置数学字体为 Arial

    if if_ploted:
        fig = plt.gcf()
        ax = plt.gca()
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        ax.set_title(ax.get_title(), size=fontsize)
        ax.set_xlabel(ax.get_xlabel(), size=fontsize)
        ax.set_ylabel(ax.get_ylabel(), size=fontsize)



def general_axes(ax,
                labelpad: int = 15, 
                tick_pad: int = 10, 
                xlabel: str = 'X-axis Label',
                ylabel: str = 'Y-axis Label',
                xlim: Union[List[float], None] = None,
                ylim: Union[List[float], None] = None,
                if_set_lim: bool = True,
                if_show_right_top_tick: bool = True):
    """
    Modifies axis properties like ticks, labels, and limits.

    Args:
        ax: The axis object to modify.
        labelpad (int): Padding for axis labels. Default is 15.
        tick_pad (int): Padding for axis ticks. Default is 10.
        xlabel (str): Label for the x-axis. Default is 'X-axis Label'.
        ylabel (str): Label for the y-axis. Default is 'Y-axis Label'.
        xlim (Union[List[float], None]): x-axis limits. Default is None.
        ylim (Union[List[float], None]): y-axis limits. Default is None.

    Returns:
        None
    """
    if isinstance(ax, np.ndarray):
        axes = ax.flatten()
    else:  
        axes = [ax]

    for axis in axes:
        axis.minorticks_on()
        axis.xaxis.set_minor_locator(AutoMinorLocator(2))
        axis.yaxis.set_minor_locator(AutoMinorLocator(2))
        axis.tick_params(which='major', direction='in', length=8, width=3.0, pad = tick_pad)
        axis.tick_params(which='minor', direction='in', length=4, width=3.0, pad = tick_pad)

        axis.set_xlabel(xlabel, labelpad=labelpad)
        axis.set_ylabel(ylabel, labelpad=labelpad)

        if if_show_right_top_tick == False:
            axis.tick_params(axis='y', which='both', right=False)
            axis.tick_params(axis='x', which='both', top=False)

        if if_set_lim:
            axis.set_xlim(xlim)
            axis.set_ylim(ylim)
        

def general_subplots_adjust(fig: Figure,
                            one_fig_wh: List[float] = [10.72, 8.205], 
                            fig_subp: List[int] = [1, 1],   
                            axes_height: float = 5.89,
                            axes_width: float = 7.31,
                            left: float = 1.918, 
                            top: float = 0.9517,
                            ):
    """
    Adjusts the figure's subplots, size, and spacing.

    Args:
        fig: The figure object to adjust.
        one_fig_wh (List[float]): Width and height of one figure. Default is [10.72, 8.205].
        fig_subp (List[int]): Number of subplots. Default is [1, 1].
        axes_height (float): Height of the axes. Default is 5.89.
        axes_width (float): Width of the axes. Default is 7.31.
        left (float): Left margin. Default is 1.918.
        top (float): Top margin. Default is 0.9517.

    Returns:
        None
    """
    fig_wh = fig_subp[1]*one_fig_wh[0], fig_subp[0]*one_fig_wh[1]
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
    axes_width: float = 7.31,
    grid_linewidth: float = 0.5,
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

    general_font(grid, grid_linewidth)

    fig_wh = fig_subp[1]*one_fig_wh[0], fig_subp[0]*one_fig_wh[1]
    fig, ax = plt.subplots(nrows=fig_subp[0], ncols=fig_subp[1], 
                           sharex=fig_sharex, figsize=(fig_wh[0], fig_wh[1]))
    general_subplots_adjust(fig, one_fig_wh, fig_subp, axes_height, axes_width, left, top)
    # can't set xlim/ylim before drawing the figure
    general_axes(ax, labelpad, tick_pad, if_set_lim=False)

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

    # 设置数学字体
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = 'Arial'

    fig_wh = fig_subp[1]*one_fig_wh[0], fig_subp[0]*one_fig_wh[1]
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

def my_plot_energy_components(fig_subp: List[int]=[4,5],
                                one_fig_wh: List[float] = [11.72, 8.205],    
                                fig_sharex: bool = False, 
                                grid: bool = True, 
                                labelpad: int = 15, 
                                tick_pad: int = 10, 
                                left: float = 1.918, 
                                top: float = 0.9517,
                                axes_height: float = 5.89,
                                axes_width: float = 7.31,
                                dic_list: list=None, 
                                label_list: list=None,
                                color_list: list=[to_rgb('#fc9d9a'),to_rgb('#8dc6ff')],
                                xlabel: str="Directory",
                                poly_fit_e: int = 3,
                                fit_count: int = 100,
                                color: str = 'tab:red',
                                atoms_number: int = 1,
                                if_tight_layout: bool= False,
                                vline_list: list= [],
                                vline_colors: list= [],
                                save_path: str='./p_post_energy_components.jpg',
                                bbox_to_anchor: tuple=(0.5, 0.95)) -> Tuple[Figure, Axes]:
    """
    Generates plots for energy components with polynomial fits and comparisons.

    Args:
        fig_subp (List[int]): Subplot grid dimensions [rows, columns]. Defaults to [4, 4].
        one_fig_wh (List[float]): Figure width and height in inches. Defaults to [11.72, 8.205].
        fig_sharex (bool): Whether subplots share the x-axis. Defaults to False.
        grid (bool): Whether to include a grid in the plots. Defaults to True.
        labelpad (int): Padding for axis labels. Defaults to 15.
        tick_pad (int): Padding for tick marks. Defaults to 10.
        left (float): Left margin of the figure. Defaults to 1.918.
        top (float): Top margin of the figure. Defaults to 0.9517.
        axes_height (float): Height of individual axes. Defaults to 5.89.
        axes_width (float): Width of individual axes. Defaults to 7.31.
        dic_list (list): List of dictionaries containing data for plotting.
        label_list (list): List of labels for each dataset in `dic_list`.
        color_list (list): List of colors for plots. Defaults to a predefined color list.
        xlabel (str): Label for the x-axis. Defaults to "Directory".
        poly_fit_e (int): Degree of the polynomial for fitting. Defaults to 3.
        fit_count (int): Number of points for the fitted curve. Defaults to 100.
        color (str): Color for the twin y-axis plots. Defaults to 'tab:red'.
        atoms_number (int): Number of atoms to normalize the data. Defaults to 1.
        if_tight_layout (bool): Whether to apply `tight_layout` for better spacing. Defaults to False.
        vline_list (list): List of x-coordinates for vertical lines. Defaults to an empty list.
        vline_colors (list): List of colors for the vertical lines. Defaults to an empty list.
        save_path (str): Path to save the plot. Defaults to './p_post_energy_components.jpg'.

    Returns:
        Tuple[Figure, Axes]: The figure and axes objects for the generated plots.

    """
    fig_subp = fig_subp
    fig, axes = my_plot(one_fig_wh, fig_subp, fig_sharex, grid, labelpad, tick_pad, 
                        left, top, axes_height, axes_width)

    tagsref=[
        "Directory", "TOTEN", "PSCENC", "TEWEN", "DENC", "EXHF",
        "XCENC", "Double1", "Double2", "Double", "EENTRO", "EBANDS", "EATOM", "Ediel_sol"
    ]
    tags = [tag.lower() for tag in tagsref]
    tags = tags[1:]

    for i, dic in enumerate(dic_list):
        valid = 0
        valide0 = 0
        for tag in tags:
            if tag == 'double':
                continue
            valid = valid + dic[tag]
        valide0 = valid - dic['eentro']
        dic['valid-toten'] = valid
        dic['valid-e0'] = valide0
    tags.insert(1, "valid-toten")
    tags.insert(2, "valid-e0")

    i = 0
    j = 0
    count = 0
    for tag in tags:
        count = count+1
        ax = axes[i, j]
        #print(count, i,j, tag)

        coeff_list=[]
        xmin = []
        xmax = []
        for i, dic in enumerate(dic_list):
            x = dic['directory']
            xmin.append(min(x))
            xmax.append(max(x))
            y = dic[tag]/atoms_number
            coeff_list.append(np.polyfit(x, y, poly_fit_e))
            ax.plot(x, y, marker='o', linestyle='', color = color_list[i], label=label_list[i])
        ax.set_xlabel(f'{xlabel}')
        ax.set_ylabel(f'{tag.upper()}'+' (eV/atom)')
        #ax.set_xlim(-0.05, 0.4)
        ax.grid(True, linestyle='--', linewidth=2)
        yfit_list=[]
        xfit = np.linspace(min(xmin), max(xmax), fit_count)
        for i, coeff in enumerate(coeff_list):
            poly = np.poly1d(coeff)
            yfit = poly(xfit)
            yfit_list.append(yfit)
            ax.plot(xfit, yfit, linestyle='--', color = color_list[i])
        #ax.axvline(x=2.79338179, color=color, linestyle='--', linewidth=3)
        if len(dic_list) == 2:
            diff = yfit_list[1] - yfit_list[0]
            diff_label = rf"$E_{{\mathrm{{{label_list[1]}}}}}-E_{{{label_list[0]}}}$"
            ax2 = ax.twinx() 
            ax2.set_ylabel(diff_label+' (eV/atom)', color=color, labelpad = 15)
            ax2.plot(xfit, diff, 'o', color = color, 
                    label=diff_label, markeredgecolor=color) 
            ax2.axhline(y=0, color=color, linestyle='--', linewidth=3)
            ax2.tick_params(axis='y', labelcolor=color, direction='in', color = color)
            ax2.spines['right'].set_color(color) 
            ax2.tick_params(which='major', direction='in', color = color, length=8, width=3.0, pad = 10)
            ax2.tick_params(which='minor', direction='in', color = color, length=4, width=3.0, pad = 10)
            ax2.set_zorder(ax.get_zorder() - 1)
            ax.spines['right'].set_visible(False)
        fig.general_modify_legend(ax.legend(loc='upper center', bbox_to_anchor=bbox_to_anchor))
        # format title
        title = tag.upper()
        ax.set_title(f'{title}', fontsize=28)
        #ax.text(0.4, 0.9, f'{tag.upper()}', transform=ax.transAxes, fontsize=28)
        for i, vline in enumerate(vline_list):
            ax.axvline(x=vline, color=vline_colors[i], linestyle='--', linewidth=3)
        i=count//fig_subp[1]
        j=count%fig_subp[1]
        if count==(fig_subp[1]*fig_subp[0]):
            break
    i = 3
    j = 0
    for temp in ['ebands', 'pscenc', 'denc', 'xcenc', 'tewen']:
        count = count + 1
        ax = axes[i, j]
        #print(count, i,j, tag)

        coeff_list=[]
        xmin = []
        xmax = []
        for i, dic in enumerate(dic_list):
            x = dic['directory']
            xmin.append(min(x))
            xmax.append(max(x))
            y = dic[temp]/atoms_number
            coeff_list.append(np.polyfit(x, y, poly_fit_e))
            ax.plot(x, y, marker='o', linestyle='', color = color_list[i], label=label_list[i])
        ax.set_xlabel(f'{xlabel}')
        ax.set_ylabel(f'{tag.upper()}'+' (eV/atom)')
        #ax.set_xlim(-0.05, 0.4)
        ax.grid(True, linestyle='--', linewidth=2)
        yfit_list=[]
        xfit = np.linspace(min(xmin), max(xmax), fit_count)
        for i, coeff in enumerate(coeff_list):
            poly = np.poly1d(coeff)
            yfit = poly(xfit)
            yfit_list.append(yfit)
            ax.plot(xfit, yfit, linestyle='--', color = color_list[i])
        #ax.axvline(x=2.79338179, color=color, linestyle='--', linewidth=3)
        if len(dic_list) == 2:
            diff = yfit_list[1] - yfit_list[0]
            diff_label = rf"$E_{{\mathrm{{{label_list[1]}}}}}-E_{{{label_list[0]}}}$"
            ax2 = ax.twinx() 
            ax2.set_ylabel(diff_label+' (eV/atom)', color=color, labelpad = 15)
            ax2.plot(xfit, diff, 'o', color = color, 
                    label=diff_label, markeredgecolor=color) 
            ax2.axhline(y=0, color=color, linestyle='--', linewidth=3)
            ax2.tick_params(axis='y', labelcolor=color, direction='in', color = color)
            ax2.spines['right'].set_color(color) 
            ax2.tick_params(which='major', direction='in', color = color, length=8, width=3.0, pad = 10)
            ax2.tick_params(which='minor', direction='in', color = color, length=4, width=3.0, pad = 10)
            ax2.set_zorder(ax.get_zorder() - 1)
            ax.spines['right'].set_visible(False)
        general_modify_legend(ax.legend(loc='upper center', bbox_to_anchor=bbox_to_anchor))
        # format title
        title = tag.upper()
        if temp == 'ebands':
            title2 = r'$E_{k}: E_{BANDS}$' + '\n' + 'Kinetic energy part of the electronic wavefunction'
        elif temp == 'pscenc':
            title2 = r'$E_{e-i}: PSCENC$' + '\n' + 'Coulomb interaction of electrons in the ionic potential field'
        elif temp == 'denc':
            title2 = r'$E_{e-e}: DENC$' + '\n' + 'Coulomb interaction between electrons'
        elif temp == 'xcenc':
            title2 = r'$E_{xc}: XCENC$' + '\n' + 'Exchange and correlation correction energy'
        elif temp == 'tewen':
            title2 = r'$E_{ii}: TEWEN$' + '\n' + 'Electrostatic Coulomb interaction between ion cores'
        ax.set_title(f'{title2}', fontsize=28)
        #ax.text(0.4, 0.9, f'{tag.upper()}', transform=ax.transAxes, fontsize=28)
        for i, vline in enumerate(vline_list):
            ax.axvline(x=vline, color=vline_colors[i], linestyle='--', linewidth=3)
        i=count//fig_subp[1]
        j=count%fig_subp[1]
        if count==(fig_subp[1]*fig_subp[0]):
            break
    if if_tight_layout:
        fig.tight_layout()
    plt.savefig(save_path, dpi=300)
    return fig, axes

def my_plot_interlayer_distance(atoms: Atoms= None, if_plot: bool = True, if_save: bool = True, save_plot_path: str = './p_post_inter_distance.jpg', 
                                save_txt_path: str = './p_post_inter_distance.txt') -> np.ndarray:
    """
    Calculates and plots the interlayer distances from the given atomic positions.

    Args:
        atoms (Atoms): Atomic structure object containing atom positions.
        if_plot (bool): If True, generates and displays the plot of interlayer distances.
        if_save (bool): If True, saves the plot and data to specified paths.
        save_plot_path (str): Path to save the generated plot image.
        save_txt_path (str): Path to save the interlayer distance data to a text file.

    Returns:
        fig, ax: for refined plot
    
    Notes:
        - The plot visualizes the interlayer distances with respect to the index of layers.
        - Data is saved in three sections: original, absolute, and normalized (relative).
    """
    atom = atoms.copy()
    positions = atom.get_positions()
    zo = positions[:,2]
    z = np.array(sorted(zo))
    z = z[1:] - z[:-1]
    zabs = z.copy()
    zref = z[len(z)//2]
    z = z/zref -1
    if if_plot:
        fig, ax = my_plot(left=3.0)
        ax.set_ylabel(rf'$d_i$/$d_{{{len(z)//2+1}}}$-1')
        ax.set_xlabel(r'index $i$ of interlayer distance $d_i$')
        ax.plot(np.arange(1, len(z) + 1), z, 'o', linestyle='-')
        ax.margins(x=0.1, y=0.1)
        if if_save:
            plt.savefig(save_plot_path, transparent=False)
    if if_save:
        with open(save_txt_path, "w") as f:
            f.write("zoriginal\n")
            for j, value in enumerate(zo, start=1):  
                f.write(f"{j:<4}  {value:>12.8f}\n") 
            f.write("\n")
            f.write("zabsolute\n")
            for j, value in enumerate(zabs, start=1):  
                f.write(f"{j:<4}  {value:>12.8f}\n") 
            f.write("\n")
            f.write("zplot\n")
            for j, value in enumerate(z, start=1):  
                f.write(f"{j:<4}  {value:>12.8f}\n") 
            f.write("\n")
    return fig, ax

def my_plot_zpositions(atoms: Atoms=None, xmargins: float = 0.1, ymargins: float = 0.1, if_save: bool=True, save_path: str = './p_post_zpositions.jpg')-> tuple:
    """
    Plots the z-positions of atoms and displays the thickness of the material.

    Args:
        atoms (Atoms, optional): The Atoms object containing the atomic positions. Defaults to None.
        xmargins (float, optional): The margin space for the x-axis. Defaults to 0.1.
        ymargins (float, optional): The margin space for the y-axis. Defaults to 0.1.
        if_save:
        save_path:

    Returns:
        tuple: A tuple containing the figure and axes objects of the plot.

    """
    positions = atoms.get_positions()
    position = positions[:, 2]

    sorted_position = np.sort(position)
    thick = float(max(position) - min(position))
    row_numbers = np.arange(1, len(sorted_position)+1)

    fig, axes = my_plot(left=3.0)
    ax = axes
    ax.plot(row_numbers, sorted_position, '-o')
    ax.set_xlabel(r'index $i$ of z-position $z_i$')
    ax.set_ylabel(r'$z_i$')
    ax.margins(x=xmargins, y=ymargins)
    ax.text(0.05, 0.95, f'Thickness: {thick:.2f} Å', transform=ax.transAxes,
        verticalalignment='top', horizontalalignment='left',
        bbox=dict(facecolor='white', edgecolor='black', boxstyle='Square,pad=0.3', linewidth=2.5, alpha=1))
    if if_save:
        plt.savefig(save_path, transparent=False)
    return fig, axes

