#from mymetal.calculate.electronic_structure.plotter import DosPlotter, BSPlotter, BSPlotterProjected, BSDOSPlotter

from mymetal.calculate.electronic_structure.universial import *
from numpy import array
import numpy as np
import re
from mymetal.universial.plot.plot import *
import matplotlib.pyplot as plt

from pymatgen.io.lobster import Icohplist
from pymatgen.core.composition import Element
from pymatgen.electronic_structure.plotter import DosPlotter
from pymatgen.io.vasp import BSVasprun, Vasprun
from pymatgen.io.lobster import Doscar

def my_plot_horizontal_vertical(xlim: list=None, ylim:list=None, one_fig_wh: list = [8.205*1.5, 8.205],
                grid: bool = True,  left: float = 1.918, top: float = 0.3, orientation: str = 'horizontal',
                xlabel: str = r'$E - E_{f}$ (eV)', ylabel: str = 'Density of States (e/eV)') -> tuple:
    """
    Create a horizontal or vertical plot with custom axis limits and labels.

    Args:
        xlim (list, optional): The x-axis limits.
        ylim (list, optional): The y-axis limits.
        one_fig_wh (list, optional): The width and height of the figure.
        grid (bool, optional): Whether to show a grid (default is True).
        left (float, optional): The left margin of the plot.
        top (float, optional): The top margin of the plot.
        orientation (str, optional): Orientation of the plot, 'horizontal' or 'vertical' (default is 'horizontal').
        xlabel (str, optional): Label for the x-axis.
        ylabel (str, optional): Label for the y-axis.

    Returns:
        tuple: The figure and axes objects for the plot.
    """
    if orientation == "horizontal":
        axes_height = one_fig_wh[1]-2.315
        axes_width = one_fig_wh[0]-3.41
        #print(left,top)
        fig, axes = my_plot(fig_subp=[1,1], one_fig_wh=one_fig_wh, axes_width = axes_width,
                            axes_height = axes_height,  grid = grid, left = left, top = top)
        axes.set_xlabel(xlabel)
        axes.set_ylabel(ylabel)
        if xlim != None:
            axes.set_xlim(xlim)
        if ylim != None:
            axes.set_ylim(ylim)
        axes.axvline(x=0, color='gray', linestyle='--')
        axes.axhline(y=0, color='gray', linestyle='--')
    else:
        axes_height = one_fig_wh[1]-2.315
        axes_width = one_fig_wh[0]-3.41
        temp = [one_fig_wh[1], one_fig_wh[0]]
        top = top
        bottom = one_fig_wh[1]-top-axes_height
        left = left
        right = one_fig_wh[0]-left-axes_width
        wspace = (left+right)/axes_width  # axes_width 7.31 inch
        hspace = (top+bottom)/axes_height  # axes_height 5.89 inch
        fig, axes = my_plot(fig_subp=[1,1], one_fig_wh=temp, axes_width = axes_height,
                            axes_height = axes_width,  grid = grid, left = bottom, top = right)     
        axes.set_ylabel(xlabel)
        axes.set_xlabel(ylabel)
        if xlim != None:
            axes.set_ylim(xlim)
        if ylim != None:
            axes.set_xlim(ylim)
        axes.axhline(y=0, color='gray', linestyle='--')
        axes.axvline(x=0, color='gray', linestyle='--')
    return fig, axes

def my_plot_orientation(ax=None, orientation:str='horizontal', x: list=None, y: list=None, color: str='darkblue', label: str=None,
                        if_fill: bool=False, if_show_nbands_nelect: bool=False, nbands: int=None, nelect: int= None, marker: str='o', linestyle: str='-')-> plt.Axes:
    """
    Plot data with specified orientation, filling options, and optional lines for band and electron numbers.

    Args:
        ax (matplotlib axis, optional): The axis to plot on.
        orientation (str, optional): Plot orientation, 'horizontal' or 'vertical'.
        x (list, optional): Data for the x-axis.
        y (list, optional): Data for the y-axis.
        color (str, optional): Line color.
        label (str, optional): Label for the plot.
        if_fill (bool, optional): Whether to fill under the curve (default is False).
        if_show_nbands_nelect (bool, optional): Whether to show the nbands and nelect lines (default is False).
        nbands (int, optional): Number of bands.
        nelect (int, optional): Number of electrons.
        marker (str, optional): Marker style.
        linestyle (str, optional): Line style.

    Returns:
        matplotlib axis: The axis object with the plot.
    """
    if orientation == "horizontal":
        ax.plot(x, y, color=color, label=label, marker=marker, linestyle=linestyle)
        colors = np.array([np.linspace(1, 0, len(y)), np.linspace(1, 0, len(y)), np.linspace(1, 1, len(y))]).T
        if if_fill:
            ax.fill_between(x, y, 0, color=color, alpha=0.2)
        if if_show_nbands_nelect:
            ax.axhline(y=nelect, color='gray', linestyle='--')
            ax.axhline(y=nbands, color='gray', linestyle='--')
    else:
        ax.plot(y, x, color=color, label=label, marker=marker, linestyle=linestyle)
        if if_fill:
            ax.fill_betweenx(x,0,y, color=color, alpha=0.2)
        if if_show_nbands_nelect:
            ax.axvline(x=nelect, color='gray', linestyle='--')
            ax.axvline(x=nbands, color='gray', linestyle='--')
    return ax

def my_plot_complete_dos(file:str = "./vasprun.xml", xlim: list=None, ylim:list=None, label:str="Total DOS", zero_at_efermi: bool = True, stack: bool = False,  sigma: float | None = None,
                one_fig_wh: list = [8.205*1.5, 8.205], if_save: bool = True, save_path: str = './dos.jpg',
                grid: bool = True,  left: float = 1.918, top: float = 0.3, orientation: str = 'horizontal', 
                show_spin: bool = False, color_list: list = ['#e41a1c', '#1b9e77','#d95f02','#7570b3','#66a61e','#e6ab02','#a6761d','#666666'],
                if_spd: bool = True, marker='', linestyle='-', if_element: bool = False, bbox_to_anchor: tuple=(0.97, 0.97),
                if_show_nbands_nelect: bool=False, if_fill: bool=True) -> tuple:
    """
    Plot the complete density of states (DOS) from a vasprun.xml file. 测试对多元素的dos, Done.

    Args:
        file (str): Path to the vasprun.xml file.
        xlim (list, optional): The x-axis limits.
        ylim (list, optional): The y-axis limits.
        label (str, optional): Label for the total DOS.
        zero_at_efermi (bool, optional): Whether to align the energy at the Fermi level.
        stack (bool, optional): Whether to stack the DOS (default is False).
        sigma (float, optional): Broadening parameter for DOS (default is None).
        one_fig_wh (list, optional): Width and height of the figure.
        if_save (bool, optional): Whether to save the plot (default is True).
        save_path (str, optional): Path to save the plot image.
        grid (bool, optional): Whether to show the grid (default is True).
        left (float, optional): Left margin.
        top (float, optional): Top margin.
        orientation (str, optional): Plot orientation, 'horizontal' or 'vertical'.
        show_spin (bool, optional): Whether to show spin-up and spin-down channels (default is False).
        color_list (list, optional): List of colors for the DOS components.
        if_spd (bool, optional): Whether to show SPD decomposition (default is True).
        marker (str, optional): Marker style for the plot.
        linestyle (str, optional): Line style.
        if_element (bool, optional): Whether to show element-wise decomposition (default is False).
        bbox_to_anchor (tuple, optional): Bounding box for legend positioning.
        if_show_nbands_nelect (bool, optional): Whether to show band and electron numbers (default is False).
        if_fill (bool, optional): Whether to fill under the DOS curve.

    Returns:
        tuple: The figure and axes objects for the plot.

    Example:
        >>> workdir = "./electronic-structure/bulk/fcc/"
        >>> ds_path = workdir + "dos/vasprun.xml"
        >>> # Plot the complete DOS
        >>> my_plot_complete_dos(ds_path, save_path=workdir + "dos.jpg", orientation="horizontal")
        >>> # Plot the integrated DOS
        >>> my_plot_idos(ds_path, orientation='horizontal', save_path=workdir + "idos.jpg")
    """

    dos = Vasprun(file, parse_potcar_file=False)
    tdos = dos.tdos
    idos = dos.idos
    pdos = dos.pdos
    nbands = dos.parameters['NBANDS']*2
    nelect = dos.parameters['NELECT']
    plotter = DosPlotter()
    plotter.add_dos(label, tdos)
    cdos = dos.complete_dos
    if if_spd:
        spd_dos = cdos.get_spd_dos()
        plotter.add_dos_dict(spd_dos)
    elif if_element:
        element_dos = cdos.get_element_dos()
        plotter.add_dos_dict(element_dos)
    dos_dict = plotter.get_dos_dict()
    #print(dos_dict.keys())

    fig, axes = my_plot_horizontal_vertical(xlim, ylim, one_fig_wh, grid,  left, top, orientation)
    # plot line count for color
    count = 0
    # loop for ["Total DOS", "s", "p", "d"]
    # loop for ["Total DOS", "Si", "O", "Au"]
    for key, value in dos_dict.items():
        #print(key)
        temp_dos = dos_dict[key]
        energies = array(temp_dos['energies'], dtype=float)  # dont minus efermi
        densities = []
        if len(temp_dos['densities'].items()) == 1:
            density = temp_dos['densities']['1']
            #print(density)
            densities.append(density)
            label_list = [key]
            if show_spin:
                label_list = [key + r' $\uparrow$']
        else:
            density1 = temp_dos['densities']['1']
            density0 = temp_dos['densities']['0']
            densities.append(density0)
            densities.append(density1)
            label_list = [key]
            if show_spin:
                label_list =  [key + r' $\uparrow$', key +r' $\downarrow$']
        # plot block
        ax = axes
        # loop for spin up and down
        for i, density in enumerate(densities):
            my_plot_orientation(ax, orientation, energies, density, color_list[count], label_list[i], if_fill, if_show_nbands_nelect, 
                                nbands, nelect, marker, linestyle)
            count += 1
    fig.general_modify_legend(ax.legend(loc='upper right', bbox_to_anchor= bbox_to_anchor))
    if if_save:
        plt.savefig(save_path, dpi=300)
    return fig, axes


def my_plot_idos(file: str = None, xlim: list=None, ylim:list=None, label:str="Integrated DOS", zero_at_efermi: bool = True, 
                 one_fig_wh: list = [8.205*1.5, 8.205], if_save: bool = True, save_path: str = './idos.jpg',
                grid: bool = True,  left: float = 1.918, top: float = 0.3, orientation: str = 'horizontal', color: str = 'darkblue',
                bbox_to_anchor: tuple=(0.97, 0.03), if_show_nbands_nelect: bool=True, if_fill:bool=True) -> tuple:
    """
    Plot the integrated density of states (IDOS) from a vasprun.xml file.

    Args:
        file (str): Path to the vasprun.xml file.
        xlim (list, optional): The x-axis limits.
        ylim (list, optional): The y-axis limits.
        label (str, optional): Label for the integrated DOS.
        zero_at_efermi (bool, optional): Whether to align the energy at the Fermi level.
        one_fig_wh (list, optional): Width and height of the figure.
        if_save (bool, optional): Whether to save the plot (default is True).
        save_path (str, optional): Path to save the plot image.
        grid (bool, optional): Whether to show the grid (default is True).
        left (float, optional): Left margin.
        top (float, optional): Top margin.
        orientation (str, optional): Plot orientation, 'horizontal' or 'vertical'.
        color (str, optional): Line color for the IDOS.
        bbox_to_anchor (tuple, optional): Bounding box for legend positioning.
        if_show_nbands_nelect (bool, optional): Whether to show band and electron numbers (default is True).
        if_fill (bool, optional): Whether to fill under the curve.

    Returns:
        tuple: The figure and axes objects for the plot.

    Example:
        >>> workdir = "./electronic-structure/bulk/fcc/"
        >>> ds_path = workdir + "dos/vasprun.xml"
        >>> # Plot the complete DOS
        >>> my_plot_complete_dos(ds_path, save_path=workdir + "dos.jpg", orientation="horizontal")
        >>> # Plot the integrated DOS
        >>> my_plot_idos(ds_path, orientation='horizontal', save_path=workdir + "idos.jpg")
    """

    from pymatgen.electronic_structure.core import Magmom, Orbital, OrbitalType, Spin
    dos = Vasprun(file, parse_potcar_file=False)
    idos = dos.idos
    pdos = dos.pdos
    nbands = dos.parameters['NBANDS']*2
    nelect = dos.parameters['NELECT']
    print(nelect)
    x = idos.energies - dos.efermi
    y = idos.densities[Spin.up]
    fig, axes = my_plot_horizontal_vertical(xlim, ylim, one_fig_wh, grid,  left, top, orientation, xlabel= r'$E - E_{f}$ (eV)', ylabel = 'Integrated Density of States (e)')
    ax = axes
    my_plot_orientation(ax, orientation, x, y, color, label, if_fill, if_show_nbands_nelect, nbands, nelect, marker='', linestyle='-')
    fig.general_modify_legend(ax.legend(loc='lower right', bbox_to_anchor=bbox_to_anchor))
    if if_save:
        plt.savefig(save_path, dpi=300)
    return fig, axes

def my_plot_element_spd_dos(file:str=None):
    # not finished
    from pymatgen.electronic_structure.core import Magmom, Orbital, OrbitalType, Spin
    dos = Vasprun(file, parse_potcar_file=False)
    cdos = dos.complete_dos
    structure = cdos.structure
    elements = structure.composition.elements
    elements_spd_dos = []
    for element in elements:
        element_spd_dos = cdos.get_element_spd_dos(element)
        
    print(elements)
    pdos = cdos.pdos
    
    #print(structure)
    # for site, p in pdos.items():
    #     print(site.specie, site.coords)
    for orbital, dos in element_spd_dos.items():
        print(orbital)

# workdir="./electronic-structure/bulk/fcc/"
# ds_path = workdir+"dos/vasprun.xml"

# my_plot_complete_dos(ds_path, save_path=workdir+"dos.jpg", orientation="horizontal")
# my_plot_idos(ds_path, orientation='horizontal', save_path=workdir+"idos.jpg")

# workdir="./electronic-structure/bulk/hcp/"
# ds_path = workdir+"dos/vasprun.xml"

# my_plot_complete_dos(ds_path, save_path=workdir+"dos.jpg", orientation="horizontal")
# my_plot_idos(ds_path, orientation='horizontal', save_path=workdir+"idos.jpg")

# ds_path = workdir+"test-dos/vasprun.xml"
# my_plot_element_spd_dos(ds_path)
# my_plot_idos(workdir+"test-dos/vasprun.xml", save_path=workdir+"idos2.jpg")
# my_plot_dos(workdir+"test-dos/vasprun.xml", save_path=workdir+"dos3.jpg", orientation="horizontal", if_spd=False, if_element=True)
# my_plot_dos(workdir+"test-dos/vasprun.xml", save_path=workdir+"dos4.jpg", orientation="vertical")