"""
energy submodule

This submodule provides functions for plotting energy components from computational simulations. It includes a function to
generate plots for various energy terms, fit polynomial curves, and compare energy differences between datasets.

Functions:
    - my_plot_energy_components: Generate plots for energy components with polynomial fits and comparisons.

"""


import numpy as np

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.colors import to_rgb


from typing import List, Tuple

from mymetal.universal.plot.general import general_modify_legend
from mymetal.universal.plot.plot    import my_plot




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
        general_modify_legend(ax.legend(loc='upper center', bbox_to_anchor=bbox_to_anchor))
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

