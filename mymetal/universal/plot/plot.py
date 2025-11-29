"""
plot submodule

This submodule contains functions for creating customized plots.

Functions:
    - my_plot: Creates a customized matplotlib figure with specific layout adjustments.
    - my_plot_brokenaxed: Creates a broken axes plot with customized layout and legend settings.
    - my_plot_modify_ploted_figure: Modify and customize the appearance of a plotted figure.
    - my_plot_colorbar: Add a colorbar to a plot.
    - general_modify_ploted_figure: Deprecated alias for `my_plot_modify_ploted_figure()`.
"""
from ase import Atoms

import os

import numpy as np

import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.colors import to_rgb
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, MaxNLocator
import matplotlib.colors as mcolors

from brokenaxes import brokenaxes

from typing import List, Tuple, Union

import itertools 

from typing import Callable

from mymetal.universal.plot.general import *
from mymetal.post.general import *

__all__ = [
    'my_plot',
    'my_plot_brokenaxed',
    'my_plot_modify_ploted_figure',
    'my_plot_colorbar',
    'general_modify_ploted_figure']


# For old version plot
def general_modify_ploted_figure(*args, **kwargs):
    """
    Deprecated alias for `my_plot_modify_ploted_figure()`.
     
    Passes all arguments directly to the new function.
    """
    return my_plot_modify_ploted_figure(*args, **kwargs)

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
    grid_linewidth: float = 2.0,
    constrained_layout: bool = False,
    if_keep_wspace_hspace: bool = False,
    wspace: float = 3.41,
    hspace: float = 2.315,
) -> Tuple[Figure, List[Axes]]:
    """Creates a customized matplotlib figure with specific layout adjustments.

    This function sets global matplotlib settings, creates a figure with 
    subplots, and adjusts their layout according to the given parameters. 
    It also enables minor ticks and customizes the axes' appearance, margins, 
    and tick settings.
    
    `if_keep_wspace_hspace` allows to keep the specified width and height
    space between subplots, adjusting the figure size accordingly.

    Args:
        one_fig_wh (List[float]): Width and height of the figure in inches 
            (default: [10.72, 8.205] = [7.31 + 3.41, 5.89 + 2.315]).
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
        grid_linewidth (float): Width of the grid lines (default: 0.5).
        constrained_layout (bool): Whether to use constrained layout for the figure
            (default: False).
        if_keep_wspace_hspace (bool): If True, keeps the specified width and height
            space between subplots (default: False).

    Returns:
        Tuple[Figure, List[Axes]]: A tuple containing:
            - fig (Figure): The created figure object.
            - axes (List[Axes]): A list of axes objects for the subplots. 
              If there's only one subplot, it returns a list with one item.
    """

    general_font(grid, grid_linewidth, constrained_layout= constrained_layout)

    if if_keep_wspace_hspace:
        one_fig_wh = [axes_width + wspace, axes_height + hspace]
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
    general_modify_legend = general_modify_legend
    fig.general_margin_bin = general_margin_bin

    return fig, ax

def my_plot_modify_ploted_figure(axes: plt.Axes,
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
                                xlim: list = None, #[-15, 10],
                                ylim: list = None, #[0, 8],
                                if_close_right_top_tick: bool = True,
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
                                save_path: str = './dos.jpg',
                                ########################### add vlines, hlines
                                vlines: list = [],
                                hlines: list = [],
                                vlines_colors: list = ['gray'],
                                hlines_colors: list = ['gray'],
                                zorder_vlines: list = [-1],
                                zorder_hlines: list = [-1],
                                 linestyle: str = '--'
                                )  -> tuple:
    """
    Modify and customize the appearance of a plotted figure. This function allows for configuring various 
    plot properties such as figure size, axis labels, tick settings, line styles, colors, legends, and more.

    Args:
        axes (plt.Axes): The axes of the figure to modify.
        grid (bool, optional): Whether to display a grid on the plot. Default is True.
        grid_linewidth (float, optional): Line width of the grid lines. Default is 0.5.
        one_fig_wh (list, optional): Size of one figure [width, height]. Default is [10.72, 8.205].
        fig_subp (List[int], optional): Subplot configuration [rows, columns]. Default is [1, 1].
        axes_height (float, optional): Height of each axis. Default is 5.89.
        axes_width (float, optional): Width of each axis. Default is 7.31.
        left (float, optional): Left margin of the subplot. Default is 1.918.
        top (float, optional): Top margin of the subplot. Default is 0.9517.
        labelpad (int, optional): Padding for axis labels. Default is 15.
        tick_pad (int, optional): Padding for axis ticks. Default is 10.
        xlabel (str, optional): Label for the x-axis. Default is 'Energies (eV)'.
        ylabel (str, optional): Label for the y-axis. Default is 'Density of States (a.u.)'.
        xlim (list, optional): Limits for the x-axis. Default is None.
        ylim (list, optional): Limits for the y-axis. Default is None.
        if_close_right_top_tick (bool, optional): Whether to display the right and top ticks. Default is True.
        loc (str, optional): Location of the legend. Default is 'upper right'.
        bbox_to_anchor (tuple, optional): Bounding box to anchor the legend. Default is (0.95, 0.95).
        ncol (int, optional): Number of columns in the legend. Default is 1.
        if_show_legend (bool, optional): Whether to display the legend. Default is True.
        y_margin (float, optional): Margin for the y-axis. Default is 0.1.
        y_nbins (int, optional): Number of bins for the y-axis. Default is 3.
        x_margin (float, optional): Margin for the x-axis. Default is 0.1.
        x_nbins (int, optional): Number of bins for the x-axis. Default is 4.
        prune (str, optional): Which axis to prune. Default is 'both'.
        if_autoscale (bool, optional): Whether to autoscale the plot. Default is True.
        change_margin (bool, optional): Whether to modify the margin. Default is False.
        remove_dashed_line (bool, optional): Whether to remove dashed lines from the plot. Default is True.
        remove_line_style (list, optional): List of line styles to remove (e.g., ['--', ':']). Default is ['--', ':'].
        color_list (list, optional): List of colors for the plot. Default is a predefined list of colors.
        if_gradient_color (bool, optional): Whether to apply gradient coloring. Default is False.
        gradient_colors (list, optional): List of colors for the gradient. Default is ['red', 'blue'].
        if_cmap_color (bool, optional): Whether to use colormap coloring. Default is False.
        cmap_color (str, optional): Colormap to use. Default is 'coolwarm'.
        if_change_color (bool, optional): Whether to change the plot color. Default is False.
        if_band_plot (bool, optional): Whether to display a band plot. Default is False.
        if_save (bool, optional): Whether to save the plot to a file. Default is True.
        save_path (str, optional): Path to save the figure. Default is './dos.jpg'.

    Returns:
        tuple: A tuple containing the modified figure and axes objects (fig, ax).

    Example:
        >>> fig, ax = general_modify_ploted_figure(axes, xlabel="Energy", ylabel="DOS", if_save=True, save_path="dos_plot.jpg")
    """
    # some pre-define value; change the fontsize, fontname, linewidth one by one
    general_font(grid, grid_linewidth, if_ploted=True)
    # useless
    ax = axes
    # pyplot
    fig, ax, xlim0, ylim0 = get_ploted_figure()
    # The type of figure
    print('The type(ax) is: ',type(ax))
    # change figsize, always well done!
    fig_wh = fig_subp[1]*one_fig_wh[0], fig_subp[0]*one_fig_wh[1]
    fig.set_size_inches(fig_wh)
    general_subplots_adjust(fig, one_fig_wh, fig_subp, axes_height, axes_width, left, top)

    # change the tick setting and label content
    if xlim == None or ylim == None:
        xlim = xlim0
        ylim = ylim0
    general_axes(ax, labelpad, tick_pad,xlabel, ylabel, xlim, ylim, True, if_close_right_top_tick)

    # change xlim and ylim to look better
    if change_margin:
        general_margin_bin(ax, y_margin, y_nbins, x_margin, x_nbins, prune, if_autoscale)
    # change the line color and remove the dashed line
    general_modify_line(ax, remove_dashed_line, remove_line_style, color_list, if_gradient_color, gradient_colors, if_cmap_color, cmap_color, if_change_color)
    # For band plot, add the vline
    if if_band_plot:
        general_modify_band_plot(ax)
    # add vlines, hlines
    general_add_vlines_hlines(ax, vlines, hlines, vlines_colors, hlines_colors, zorder_vlines, zorder_hlines, linestyle)
    # re-generate the legend
    current_legend = ax.get_legend()
    ax.legend().remove()

    if if_show_legend:
        if current_legend != None:
            labels = [text.get_text() for text in current_legend.get_texts()]
            general_modify_legend(ax.legend(labels, loc=loc, bbox_to_anchor=bbox_to_anchor, ncol=ncol))
        else:
            general_modify_legend(ax.legend(loc=loc, bbox_to_anchor=bbox_to_anchor, ncol=ncol))

    # save figure
    if if_save:
        plt.savefig(save_path, dpi=300)
    return fig, ax

def my_plot_colorbar(original_figsize: tuple=(10.72, 8.205),
                     axsize: tuple=(7.31, 5.89),
                     colorbar_size: tuple=(0.5, 5.89),
                     pad: float = 0.5,
                     left: float = 1.918, 
                     top: float = 0.5,
                     layout: str="constrained",
                     grid: bool=True,
                     grid_linewidth: float=0.5):
    """
    Creates a figure with a main axis and a right-aligned colorbar using absolute positioning.

    Args:
        original_figsize (tuple): Base figure size before adding colorbar (width, height) in inches.
        axsize (tuple): Size of the main axes (width, height) in inches.
        colorbar_size (tuple): Size of the colorbar axes (width, height) in inches.
        pad (float): Padding between main axis and colorbar in inches.
        left (float): Distance from the figure's left edge to the main axis in inches.
        top (float): Distance from the figure's top edge to the top of the main axis in inches.
        layout (str): Layout engine to use, e.g., 'constrained'.
        grid (bool): Whether to show gridlines on the main axis.
        grid_linewidth (float): Line width of the grid.

    Returns:
        Tuple[Figure, Tuple[Axes, Axes]]: The figure and a tuple containing the main and colorbar axes.

    Usage:
        ```python
        import numpy as np
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors

        fig, axes = my_plot_colorbar(original_figsize = [10.72, 10.72], axsize = [7.31, 7.31], colorbar_size=[0.5, 7.31],
                                    grid = False)
                                    
        fig, axes = my_plot_colorbar(layout="none", grid= False)
        ax = axes[0]   # Main axis for plotting
        ax2 = axes[1]  # Colorbar axis

        # main plot
        keys = np.arange(2.70, 2.91, 0.01).tolist()
        keys = [round(key, 2) for key in keys]
        gradient_colors = generate_gradient_colors(if_cmap_color=True, num_colors=len(keys))
        x = [0, 1/3, 0.5, 5/6, 1]
        labels = ['FCC', r'$\dfrac{1}{3}$', 'DHCP' , r'$\dfrac{5}{6}$','HCP']

        # 创建ScalarMappable，用于颜色条
        norm = mcolors.Normalize(vmin=min(keys), vmax=max(keys))  # 范围对应keys
        sm = plt.cm.ScalarMappable(cmap='coolwarm', norm=norm)  # 使用自定义的颜色映射（如 coolwarm）
        sm.set_array([])  # 空的数组，用于colorbar

        for key, color in zip(keys, gradient_colors):
            y = transition_path_dict[key]["energy"]/transition_path_dict[key]["atoms_number"] # eV/atom
            y = (y-y[0])*1e3 # meV/atom
            ax.plot(x, y, label=f'{key:.2f} Å', marker='o', linestyle=':',
                    color=color)#, markerfacecolor='none')
            ax.set_xlabel(r'Hexagonality $i$')
            ax.set_ylabel(r'$E^{i}-E^{0}$ (meV/atom)')

        # 添加 colorbar
        cbar = fig.colorbar(sm, cax = ax2)
        cbar.set_label('In-plane lattice constant (Å)', rotation=270, labelpad=45)
        ```
    """
    left = left
    axes_width = axsize[0]
    axes_height = axsize[1]
    right = original_figsize[0]-left-axes_width
    top = top
    bottom = original_figsize[1]-top-axes_height
    full_figsize =[original_figsize[0]+pad+colorbar_size[0],
                  original_figsize[1]]
    
    general_font(grid, grid_linewidth)
    fig = plt.figure(figsize=full_figsize, layout=layout)
    # figure
    ax1 = fig.add_axes([left/full_figsize[0], bottom/full_figsize[1],
                         axes_width/full_figsize[0], axes_height/full_figsize[1]])
    general_axes(ax1)
    # colorbar
    ax2 = fig.add_axes([(left+axes_width+pad)/full_figsize[0],
                        bottom/full_figsize[1],
                        colorbar_size[0]/full_figsize[0],
                        colorbar_size[1]/full_figsize[1]])
    ax2.tick_params(direction='in', width=3, length=8)
    return fig, (ax1, ax2)
