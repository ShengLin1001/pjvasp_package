"""
general submodule

This module contains general functions for customizing Matplotlib plots.

Functions:
    - check_font_size: Check the font size of labels and ticks on a plot.
    - check_axes_size: Check the size of axes in a plot.
    - check_all_rcParams: Check all Matplotlib rcParams settings.
    - get_ploted_figure: Get the current figure and axis objects.
    - general_modify_band_plot: Modify the appearance of a band structure plot.
    - generate_gradient_colors: Generate a list of gradient colors.
    - general_modify_line: Modify the appearance of lines in a plot.
    - general_add_vlines_hlines: Add vertical and horizontal lines to a plot.
    - general_modify_legend: Modify the appearance of a legend.
    - general_margin_bin: Adjust the margins and binning for the x and y axes.
    - general_font: Customize global plotting settings for Matplotlib.
    - general_axes: Modifies axis properties like ticks, labels, and limits.
    - general_subplots_adjust: Adjusts the figure's subplots, size, and spacing.
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

def check_axes_size(fig, ax1):
     # 1. 获取 ax1 在 Figure 中的位置和大小（归一化坐标，单位 0~1）
    bbox = ax1.get_position()
    print("ax1 in figure (relative):", bbox)

    # 2. 获取 figure 的物理尺寸（英寸）
    fig_w, fig_h = fig.get_size_inches()

    # 3. 计算 ax1 的实际尺寸（英寸）
    ax1_w_inch = bbox.width * fig_w
    ax1_h_inch = bbox.height * fig_h
    print(f"ax1 size in inches: {ax1_w_inch:.3f} x {ax1_h_inch:.3f}")

    # 4. 如果需要换算为像素（可选）
    dpi = fig.dpi
    print(f"ax1 size in pixels: {ax1_w_inch * dpi:.1f} x {ax1_h_inch * dpi:.1f}")

def check_all_rcParams():
    for key, value in plt.rcParams.items():
        print(f"{key:40} = {value}")

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
    general_add_vlines_hlines(ax, xticks)
    return ax

def generate_gradient_colors(start_color: str = 'red', end_color: str = 'blue', num_colors: int = 2,
                             if_cmap_color: bool = False, cmap_color: str = 'coolwarm'):
    if if_cmap_color:
        cmap = plt.get_cmap(cmap_color)
        values = np.linspace(0, 1, num_colors)
        color_list = cmap(values)
    else:
        cmap = mcolors.LinearSegmentedColormap.from_list("gradient", [start_color, end_color])
        color_list = [cmap(i / (num_colors - 1)) for i in range(num_colors)]
    return color_list

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
        color_list = generate_gradient_colors(gradient_colors[0], gradient_colors[1], len(after_removed_lines),
                                              if_cmap_color, cmap_color)
    color_cycle = cycle(color_list)
    if if_change_color:
        for index, line in enumerate(after_removed_lines):
            line.set_color(next(color_cycle))
    return ax

def general_add_vlines_hlines(  ax: plt.Axes,
                                vlines: list = [],
                                hlines: list = [],
                                vlines_colors: list = ['gray'],
                                hlines_colors: list = ['gray'],
                                zorder_vlines: list = [-1],
                                zorder_hlines: list = [-1],
                                linestyle: str='--'):
    # add vlines, hlines
    color_cycle = cycle(vlines_colors)
    zorder_cycles=cycle(zorder_vlines)
    for vline in vlines:
        ax.axvline(x=vline, color=next(color_cycle), zorder=next(zorder_cycles), linestyle= linestyle)
    color_cycle = cycle(hlines_colors)
    zorder_cycles=cycle(zorder_hlines)
    for hline in hlines:
        ax.axhline(y=hline, color=next(color_cycle), zorder=next(zorder_cycles),  linestyle= linestyle)

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

def general_font( grid: bool = True, grid_linewidth: float = 0.5, grid_linestyle = '--',
                 if_ploted: bool = False,
                  fontsize: int = 28, markersize: int = 20, linewidth: int = 3, 
                  markeredgewidth: int = 3, legend_fontsize: int = 24, 
                  markerfacecolor: str = 'white', fontname: str = 'Arial',
                  loc: str = 'upper right', frameon: bool = True, 
                  borderpad: float=0.0, labelspacing: float=0.5, columnspacing:float=0.5,
                  constrained_layout: bool = False) -> None:
    """
    Customize global plotting settings for Matplotlib, including font styles, line properties, 
    grid appearance, and legend settings. Optionally, apply these settings to an existing plot.

    Args:
        grid (bool, optional): Whether to display grid lines on the plot. Default is True.
        grid_linewidth (float, optional): Line width for grid lines. Default is 0.5.
        grid_linestyle (str, optional): Style of grid lines (e.g., '--', '-'). Default is '--'.
        if_ploted (bool, optional): Whether to apply the settings to the current plot if one exists. Default is False.
        fontsize (int, optional): Font size for text elements such as labels, title, and ticks. Default is 28.
        markersize (int, optional): Size of markers in the plot. Default is 20.
        linewidth (int, optional): Width of plot lines. Default is 3.
        markeredgewidth (int, optional): Width of marker edges. Default is 3.
        legend_fontsize (int, optional): Font size for the legend text. Default is 24.
        markerfacecolor (str, optional): Color of the marker face. Default is 'white'.
        fontname (str, optional): Name of the font to use for text elements. Default is 'Arial'.
        loc (str, optional): Location of the legend on the plot. Default is 'upper right'.
        frameon (bool, optional): Whether to display a frame around the legend. Default is True.
        borderpad (float, optional): Padding around the legend. Default is 0.0.
        labelspacing (float, optional): Spacing between legend labels. Default is 0.5.
        columnspacing (float, optional): Spacing between legend columns. Default is 0.5.

    Returns:
        None: This function modifies the global plotting settings and optionally the current plot.

    Example:
        >>> general_font(fontsize=20, legend_fontsize=18, if_ploted=True)
    """
    plt.rcParams['font.family'] = fontname                 #
    plt.rcParams['font.size'] = fontsize                   #
    plt.rcParams['axes.linewidth'] = linewidth             #
    plt.rcParams['axes.grid'] = grid                       #
    plt.rcParams['grid.linestyle'] = grid_linestyle        #
    plt.rcParams['grid.linewidth'] = grid_linewidth        #
    plt.rcParams["savefig.transparent"] = 'True'           # don't need change
    plt.rcParams['lines.linewidth'] = linewidth            #
    plt.rcParams['lines.markersize'] = markersize          #
    plt.rcParams['lines.markeredgewidth'] = markeredgewidth #
    plt.rcParams['lines.markerfacecolor'] = markerfacecolor #
    
    # layout
    plt.rcParams['figure.constrained_layout.use'] = constrained_layout

    # legend     will done, because the legend will be re-drawed
    plt.rcParams['legend.loc'] = loc
    plt.rcParams['legend.fontsize'] = legend_fontsize
    plt.rcParams['legend.frameon'] = frameon
    plt.rcParams['legend.borderpad'] = borderpad
    plt.rcParams['legend.labelspacing'] = labelspacing
    plt.rcParams['legend.columnspacing'] = columnspacing

    # Math
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = 'Arial'  # 设置数学字体为 Arial

    if if_ploted:
        fig = plt.gcf()
        axis = plt.gca()
        if isinstance(axis, np.ndarray):
            axes = axis.flatten()
        else:  
            axes = [axis]
        for ax in axes:
            # change the font in tick_params
            plt.xticks(fontsize=fontsize, fontname= fontname)
            plt.yticks(fontsize=fontsize, fontname= fontname)
            # change grid
            ax.grid(linewidth=grid_linewidth, linestyle = grid_linestyle)
            ax.grid(grid)
            # re-generate title/xlabel/ylabel about fontsize, fontname
            ax.set_title(ax.get_title(), size=fontsize, fontname=fontname)
            ax.set_xlabel(ax.get_xlabel(), size=fontsize, fontname=fontname)
            ax.set_ylabel(ax.get_ylabel(), size=fontsize, fontname=fontname)
            # change the linewidth in ax.
            for temp in ['top', 'bottom', 'right', 'left']:
                ax.spines[temp].set_linewidth(linewidth) 
            # change line setting
            for line in ax.get_lines():
                line.set_linewidth(linewidth)
                line.set_markersize(markersize)
                line.set_markeredgewidth(markeredgewidth)
                line.set_markerfacecolor(markerfacecolor)
        
def general_axes(ax,
                labelpad: int = 15, 
                tick_pad: int = 10, 
                xlabel: str = 'X-axis Label',
                ylabel: str = 'Y-axis Label',
                xlim: Union[List[float], None] = None,
                ylim: Union[List[float], None] = None,
                if_set_lim: bool = True,
                if_close_right_top_tick: bool = False,
                if_ax2: bool = False,
                color: str = 'black',):
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

    # change the top/bottom/top/right linewidth in general_font() function.
    for axis in axes:
        axis.minorticks_on()
        axis.xaxis.set_minor_locator(AutoMinorLocator(2))
        axis.yaxis.set_minor_locator(AutoMinorLocator(2))
        axis.tick_params(which='major', direction='in', length=8, width=3.0, pad = tick_pad)
        axis.tick_params(which='minor', direction='in', length=4, width=3.0, pad = tick_pad)

        axis.set_xlabel(xlabel, labelpad=labelpad)
        axis.set_ylabel(ylabel, labelpad=labelpad)

        if if_close_right_top_tick:
            axis.tick_params(axis='y', which='both', right=False)
            axis.tick_params(axis='x', which='both', top=False)

        if if_set_lim:
            axis.set_xlim(xlim)
            axis.set_ylim(ylim)
        
        if if_ax2:
            axis.set_ylabel(ylabel, color=color, labelpad = labelpad)
            axis.tick_params(axis='y', labelcolor=color, direction='in', color = color)
            axis.spines['right'].set_color(color) 
            axis.tick_params(axis='y', which='both', color = color)

def general_subplots_adjust(fig: Figure,
                            one_fig_wh: List[float] = [10.72, 8.205], 
                            fig_subp: List[int] = [1, 1],   
                            axes_height: float = 5.89,
                            axes_width: float = 7.31,
                            left: float = 1.918, 
                            top: float = 0.9517):
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
    fig_wh = [fig_subp[1]*one_fig_wh[0], fig_subp[0]*one_fig_wh[1]]
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
    