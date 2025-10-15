"""
workflow plot submodule

This submodule provides functions for plotting various workflow-related data, including convergence tests, NEB trajectories,
and stretch analysis. It includes functions to visualize energy convergence, NEB energy and force profiles, atomic trajectories,
and stretch energy curves.

Functions:
    - my_plot_convergence: Plot convergence test results for energy cutoff or k-point grids.
    - my_plot_neb: Plot NEB (Nudged Elastic Band) energy and force profiles.
    - my_plot_neb_xy: Plot atomic trajectories in 2D with customizable visualization options.
    - my_plot_stretch: Fit quadratic energy-stretch curve, extract equilibrium values, and plot energy and c/a ratio.
"""

import os

import itertools 

import numpy as np

import pandas as pd

import matplotlib.pyplot as plt

from typing import Callable

from mymetal.post.general import my_ployfit
from mymetal.universal.plot.plot import my_plot
from mymetal.universal.plot.general import general_adjust_text, general_modify_legend


# For workflow convergence test
def my_plot_convergence(x: list = None, 
                        y: list = None,

                        encuts: bool=False, 
                        kpoints: bool=False,

                        if_scalex_auto: bool=True,
                        scalex: float=1.0,

                        if_difference: bool=False,
                        encuts_labels_dif: list = ['Energy cutoff (eV)', r'$E-E^{\text{Ref}}$ (meV per atom)'],
                        kpoints_labels_dif: list = ['K-point grid', r'$E-E^{\text{Ref}}$ (meV per atom)'],
                        if_text: bool=False,
                        format: str='.1f',
                        xytext: tuple = (0, 10),
                        text_threshold: float=2.0,
                        text_fontsize: int=20,
                        if_mask: bool=False,
                        mask_condition: Callable[[list], list] = lambda y: np.abs(y) > 1,
                        mask_color: str='red',

                        encuts_labels: list = ['Energy cutoff (eV)', 'Energy (eV per atom)'],
                        kpoints_labels: list = ['K-point grid', 'Energy (eV per atom)'],
                        kpoints_ylabel_rotation: int=45,
                        kpoints_ylabel_fontsize: int=20,

                        marker: str='o',
                        linestyle: str='-',
                        color: str='b',
                        label: str='Total Energy',
                        if_save: bool=False,
                        savefile: str='p_post_convergence.pdf',
                        dpi: int=300,
                        ) -> tuple:
    """
    Plot convergence test results for energy cutoff or k-point grids.

    This function visualizes convergence data, such as total energies or 
    energy differences (relative to reference) versus energy cutoffs or 
    k-point grid sizes. It supports labeling, annotation, automatic x-axis 
    scaling, and highlighting specific points with a mask.

    Args:
        x (list, optional): X-axis data, e.g., cutoff values or k-point grids.
        y (list, optional): Y-axis data, e.g., total energy (eV/atom).
        encuts (bool): If True, interpret x as energy cutoffs.
        kpoints (bool): If True, interpret x as k-point grids.
        if_scalex_auto (bool): Whether to auto-scale the plot width based on x.
        scalex (float): Manual scaling factor for plot width.
        if_difference (bool): Plot energy differences relative to reference.
        encuts_labels_dif (list): X and Y labels for encut difference plot.
        kpoints_labels_dif (list): X and Y labels for k-point difference plot.
        if_text (bool): If True, annotate data points with text.
        format (str): Format for annotation text (e.g., '.1f').
        xytext (tuple): Offset for annotation text in points.
        text_threshold (float): Only annotate points with abs(y) below this.
        text_fontsize (int): Font size for annotations.
        if_mask (bool): Whether to highlight points matching `mask_condition`.
        mask_condition (Callable): A function that returns a boolean mask on y.
        mask_color (str): Color for masked points.
        encuts_labels (list): X and Y labels for encut plot.
        kpoints_labels (list): X and Y labels for k-point plot.
        kpoints_ylabel_rotation (int): Rotation angle for x-tick labels (k-points).
        kpoints_ylabel_fontsize (int): Font size for x-tick labels (k-points).
        marker (str): Marker style for data points.
        linestyle (str): Line style.
        color (str): Line color.
        label (str): Label for legend.

    Returns:
        tuple: (fig, ax) Matplotlib Figure and Axes objects.
    """
    # y units: eV per atom
    if if_scalex_auto:
        scalex = max(1, len(x)/10)
    axes_width = 7.31 * scalex
    one_fig_width = 1.918 + axes_width + 1.492
    fig, axes = my_plot(axes_width = axes_width, one_fig_wh=[one_fig_width, 8.205],)
    ax = axes
    
    if encuts:
        if if_difference:
            xlabel = encuts_labels_dif[0]
            ylabel = encuts_labels_dif[1]
        else:
            xlabel = encuts_labels[0]
            ylabel = encuts_labels[1]
    elif kpoints:
        if if_difference:
            xlabel = kpoints_labels_dif[0]
            ylabel = kpoints_labels_dif[1]
        else:
            xlabel = kpoints_labels[0]
            ylabel = kpoints_labels[1]
        labels = []
        for temp in x:
            tempx = r"${0}\times{1}\times{2}$".format(int(temp[0]), int(temp[1]), int(temp[2]))
            labels.append(tempx)
        x = np.arange(len(labels))
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=kpoints_ylabel_rotation, fontsize=kpoints_ylabel_fontsize)
    
    if if_difference:
        y = y - y[-1] # make the last value zero
        y = y * 1e3   # meV per atom
        if if_text:
            for i in range(len(x)):
                if abs(y[i]) < text_threshold:
                    ax.annotate(f"{y[i]:{format}}", 
                                (x[i], y[i]),   
                                textcoords="offset points",
                                xytext=xytext,
                                va = 'bottom',
                                ha = 'center',
                                fontsize = text_fontsize)

    ax.plot(x, y, marker=marker, linestyle=linestyle, color=color, label=label)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if if_mask:
        mask = mask_condition(y)
        ax.plot(x[mask], y[mask], marker=marker, color=mask_color, linestyle='None')
    
    if if_save:
        plt.savefig(savefile, dpi=dpi)
    return fig, ax


# for workflow NEB trajectory
def my_plot_neb(nebdf: pd.DataFrame = None, nebefdf: pd.DataFrame =  None, 
                spline_df: pd.DataFrame =  None, extsdf: pd.DataFrame =  None, 
                mysplinedf: pd.DataFrame = None,
                natoms: int = 1,
                if_save: bool=False,
                savefile: str='p_post_neb.pdf',
                dpi: int=300,) -> tuple:
    """Plot NEB (Nudged Elastic Band) energy and force profiles.

    This function visualizes energy and force data from NEB calculations,
    including original data, spline fits, and optionally user-defined splines.

    Args:
        nebdf (pd.DataFrame): DataFrame containing original NEB coordinates and forces.
        nebefdf (pd.DataFrame): DataFrame with NEB energy and max force data.
        spline_df (pd.DataFrame): DataFrame with spline-interpolated energy and force.
        extsdf (pd.DataFrame): DataFrame with extrema point data.
        mysplinedf (pd.DataFrame): DataFrame with custom spline energy and force.
        natoms (int): Number of atoms for per-atom normalization.
        if_save (bool): Whether to save the figure as a file.
        savefile (str): Filename to save the plot.
        dpi (int): Resolution of the saved plot.

    Returns:
        tuple: (fig, axes) of the generated matplotlib figure and axes array.
    """

    # read data
    # orginal data
    image = nebdf['image'] 
    coord = nebdf['dist_cum(Å)']
    coord_scale = coord / coord.iloc[-1] # scale to 1
    energy = nebefdf['energy(eV)'] 
    energy_rel = ( energy - energy[0] ) * 1e3 / natoms # meV per atom
    force = nebdf['force_b(eV/Å)']
    maxforce = nebefdf['force_max(eV/Å)']
    # spline data
    spline_coord = spline_df['dist_cum(Å)']
    spline_coord_scale = spline_coord / spline_coord.iloc[-1] # scale to 1
    spline_energy = spline_df['energy_rel(eV)'] * 1e3 / natoms # meV per atom
    spline_force = spline_df['force_b(eV/Å)']
    # exts data
    exts_coord = extsdf['image'] #/ coord.iloc[-1]
    exts_energy = extsdf['energy_rel(eV)'] * 1e3 / natoms # meV per atom
    # my spline data
    my_spline_coord_scale = mysplinedf['dist_cum(Å)'] / mysplinedf['dist_cum(Å)'].iloc[-1] # scale to 1
    my_spline_energy_rel = ( mysplinedf['energy(eV)'] - mysplinedf['energy(eV)'][0]) * 1e3 / natoms # meV per atom
    my_spline_force = mysplinedf['force_b(eV/Å)']

    # plot
    fig, axes = my_plot(fig_subp=[2,3],fig_sharex =False)
    ax = axes[0, 0]
    ax.plot(coord_scale, energy_rel, marker = 'o', linestyle = '')
    ax.plot(spline_coord_scale, spline_energy)
    ax.set_xlabel('Scaled coordinate')
    ax.set_ylabel('Energy (meV per atom)')

    ax = axes[0, 1]
    ax.plot(coord, energy_rel, marker = 'o', linestyle = '')
    ax.plot(spline_coord, spline_energy)
    ax.set_xlabel('Coordinate (Å)')
    ax.set_ylabel('Energy (meV per atom)')

    ax = axes[0, 2]
    ax.plot(coord_scale, force, marker = 'o', linestyle = '')
    ax.plot(spline_coord_scale, spline_force)
    ax.set_xlabel('Scaled coordinate')
    ax.set_ylabel('Force (eV/Å)')

    ax = axes[1, 0]
    ax.plot(coord, maxforce, marker = 'o')
    ax.set_xlabel('Coordinate (Å)')
    ax.set_ylabel('Max force (eV/Å)')

    ax = axes[1, 1]
    ax.plot(coord_scale, energy_rel, marker = 'o', linestyle = '')
    ax.plot(my_spline_coord_scale, my_spline_energy_rel)
    ax.set_xlabel('Scaled coordinate')
    ax.set_ylabel('Energy (meV per atom)')

    ax = axes[1, 2]
    ax.plot(coord_scale, force, marker = 'o', linestyle = '')
    ax.plot(my_spline_coord_scale, my_spline_force)
    ax.set_xlabel('Scaled coordinate')
    ax.set_ylabel('Force (eV/Å)')

    if if_save:
        plt.savefig(savefile, dpi=dpi)

    return fig, axes


def my_plot_neb_xy(xy_list: list = None, z_list: list = None, figxy_wh_lim: list = None, fig_subp: list = None,if_every_atom: bool =True, if_every_frame: bool = False,
               if_every_atom_filename: str = 'neb_trajectory_xy_diff_atoms.png', if_every_frame_filename: str = 'neb_trajectory_xy_diff_frames.png',
               cellxypoints_special: list = None, if_save: bool = True, save_dir: str = './', 
               color_list: list = ['red', 'blue', 'green', 'orange', 'purple', 
                                                    'brown', 'pink', 'gray', 'olive', 'cyan', 'lime', 
                                                    'teal', 'lavender', 'tan', 'salmon', 'gold', 'lightcoral', 
                                                    'skyblue', 'darkblue', 'darkred', 'darkgreen', 'darkorange', 
                                                    ],
                boundary_color: str = 'red',
                special_points_color: str = 'blue',
                num_ax_legend: int = None, loc: str ='center', bbox_to_anchor: tuple =(0.5, 0.5), ncol: int =2,
                xlabel: str = 'X (Å)', ylabel: str = 'Y (Å)',
                alpha: list = [1, 1, 0.7],
                if_show_frame_overlap: bool = False,
                overlap_groups: list = None,
                xytext: list = (0, 10)
                                                    ) -> tuple:
    """Plots atomic trajectories in 2D with customizable visualization options.

    Args:
        xy_list: List containing [positions_xy, delta_xy, dist_xy] and positions_z arrays.
        z_list: Z-coordinate data (unused in 2D plot).
        figxy_wh_lim: Figure dimensions and limits [cell_points, [width, height], axes_height, lims].
        fig_subp: Subplot grid layout [rows, cols].
        if_every_atom: If True, plots each atom's trajectory across frames.
        if_every_frame: If True, plots each frame's atomic configuration.
        if_every_atom_filename: Filename for saving atom-wise plots.
        if_every_frame_filename: Filename for saving frame-wise plots.
        cellxypoints_special: Special points in the xy plane for highlighting.
        if_save: If True, saves the plot to a file.
        save_dir: Directory to save the plot.
        color_list: List of colors for plotting.
        boundary_color: Color for the boundary of the cell.
        special_points_color: Color for special points in the xy plane.
        num_ax_legend: Number of axes in the legend.
        loc: Location of the legend.
        bbox_to_anchor: Bounding box to anchor the legend.
        ncol: Number of columns in the legend.
        xlabel: Label for the x-axis.
        ylabel: Label for the y-axis.
        alpha: List of alpha values for different plot elements.
        if_show_frame_overlap: If True, shows overlapping groups in the plot.
        overlap_groups: List of groups to highlight in the plot.
        xytext: Offset for text annotations.
        **kwargs: Additional plotting parameters.

    Returns:
        tuple: (fig, axes) matplotlib figure and axes objects.

    Raises:
        ValueError: If neither if_every_atom nor if_every_frame is True.

    Example:
        ```python
        # Example usage of my_plot_xy for visualizing atomic trajectories

        # 1. Prepare atomic data from extxyz file
        atomlist = extxyz_to_atomlist(file)
        cell = np.array(atomlist[0].get_cell())
        cellxy = cell[:2, :2]  # Extract 2D cell vectors

        # 2. Convert special points from direct to Cartesian coordinates
        cellxypoints_special_direct = [[0, 0], [1/3, 2/3], [2/3, 1/3]]
        cellxypoints_special = np.dot(np.array(cellxypoints_special_direct), cellxy)

        # 3. Calculate trajectory data
        [positions_xy_list, delta_xy_list, dist_xy_list], [positions_z_list] = get_delta_dist_list(atomlist)

        # 4. Get figure dimensions and limits
        cellxypoints_xy, [figxy_width, figxy_height], axes_height_xy, figxy_lim = get_figxy_wh_lim(atomlist)

        # 5. Identify overlapping atoms
        overlap_pairs, overlap_groups = save_overlap_pairs_groups(
            positions_xy_list, 
            distance_threshold=distance_threshold, 
            save_dir=save_dir
        )

        # 6. Plot trajectories
        if if_plot:
            # Please see mymetal.post.neb.analyze_neb_trajectory for more details on the parameters.
            # Example 1: Plot trajectories by atom (each subplot shows one atom's path through frames)

            # Example 2: Plot trajectories by frame (each subplot shows all atoms in one frame)

            # Example 3: Basic frame-by-frame plot without overlap annotations
        ```
    """
    [positions_xy_list, delta_xy_list, dist_xy_list], [positions_z_list] = xy_list, z_list
    [cellxypoints_xy, [figxy_width, figxy_height], axes_height_xy, figxy_lim] = figxy_wh_lim
    boundary = cellxypoints_xy
    special_points = cellxypoints_special
    num_atoms = len(positions_xy_list)
    num_frames = len(positions_xy_list[0])
    color_cycle = itertools.cycle(color_list)

    if if_every_atom:
        num_ax = num_atoms
        num_plot = num_frames
        filename = if_every_atom_filename
        fig_subp = [1, num_atoms+1] if fig_subp is None else fig_subp
        title = 'Atom'
        num_ax_legend = num_atoms + 1 if num_ax_legend is None else num_ax_legend
    elif if_every_frame:
        num_ax = num_frames
        num_plot = num_atoms # num_points_in_one_ax = num_atoms
        filename = if_every_frame_filename
        fig_subp = [1, num_frames+1] if fig_subp is None else fig_subp
        title = 'Frame'
        num_ax_legend = num_frames + 1 if num_ax_legend is None else num_ax_legend
    else:
        raise ValueError('Either if_every_atom or if_every_frame must be True.')
    
    if if_show_frame_overlap:
        my_plot_group_point = []
        # every frame loop
        for overlap_group in overlap_groups:
            # if len(overlap_group) == 0:
            #     my_plot_point.append([])
            #     continue
            # every group loop in the frame
            # can deal with the case that there is no overlap group in the frame
            my_plot_group_point.append([np.max(group) for group in overlap_group])

    fig_xy, axes_xy = my_plot(one_fig_wh = [figxy_width, figxy_height], axes_height = axes_height_xy, fig_subp=fig_subp, fig_sharex=False, top = 0.5)
    axlist = axes_xy.ravel() # 1-D list, (2,2) => 4
    if len(axlist) < num_ax+1:
        raise ValueError(f'Not enough axes in the figure. Expected at least {num_ax} + 1, got {len(axlist)}.')
    axmainlist = axlist[:num_ax]

    for i, ax in enumerate(axmainlist):
        
        ax.plot(boundary[:,0], boundary[:,1], linestyle='-', color=boundary_color, zorder=10, alpha=alpha[0])
        ax.plot(special_points[:, 0], special_points[:, 1], marker = 'x', linestyle='', color=special_points_color, zorder = 11, alpha=alpha[1])

        if if_show_frame_overlap:
            my_plot_group_point = []
            # every frame loop
            overlap_group = overlap_groups[i]
            # every group loop in the frame
            # can deal with the case that there is no overlap group in the frame
            special_group_points = [np.max(group) for group in overlap_group]

        texts = []
        for j in range(num_plot):
            
            # positions_xy_list[i][j], ith atom, jth frame
            if if_every_atom:
                # i is the atom index, j is the frame index
                positions_xy = positions_xy_list[i]
                label = f'Frame {j+1}'
            elif if_every_frame:
                # i is the frame index, j is the atom index
                positions_xy = np.array([positions_xy_list[j][i] for j in range(num_plot)])
                label = f'Atom {j+1}'
            color = next(color_cycle)
            ax.plot(positions_xy[j, 0], positions_xy[j, 1], marker='o', linestyle='', label = label, color = color, alpha = alpha[2])

            if if_show_frame_overlap:
                for temp_i, point in enumerate(special_group_points):
                    if point == j:
                        group = overlap_group[temp_i]
                        group_str = '-'.join(str(atom+1) for atom in group)
                        text = ax.text(positions_xy[j, 0], positions_xy[j, 1], group_str, 
                                        ha='center', va='center', zorder = 100)
                        texts.append(text)
                        # ax.annotate(group_str, xy=(positions_xy[j, 0], positions_xy[j, 1]),
                        #             xytext=xytext, textcoords='offset points',)
                        #print(group_str)

        ax.set_xlim(figxy_lim[0])
        ax.set_ylim(figxy_lim[1])

        if if_show_frame_overlap:
            general_adjust_text(
                texts = texts,
                ax=ax,
                x = positions_xy[:, 0],
                y = positions_xy[:, 1],
                ensure_inside_axes=True, 
                #iter_lim= 1000, 
                #expand= (1.3, 1.3),
                #force_static=100,          # 控制点与标签的排斥力
                #force_text=0.3,            # 控制标签间的排斥力
                if_avoid_markers=True,
            )
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(f'{title} {i+1}')

    # For adding legend
    axblanklist = axlist[num_ax:]
    for ax in axblanklist:
        ax.axis('off')
    ax_legend = axlist[num_ax_legend]
    handles, labels = axmainlist[-1].get_legend_handles_labels()
    general_modify_legend(ax_legend.legend(handles, labels, loc=loc, bbox_to_anchor=bbox_to_anchor, ncol=ncol))

    if if_save:
        plt.savefig(os.path.join(save_dir, filename), dpi=300, facecolor='white')
        plt.savefig(os.path.join(save_dir, filename.replace('.png', '.pdf')), dpi=300)
    
    return fig_xy, axes_xy



# for workflow Stretch trajectory
def my_plot_stretch(
                jobn: list = None,
                Etot: list = None,
                natoms: int = 1,
                stretch_type: str = '',
                rvectors_ref = None,
                lca: list = None,
                if_save: bool=False,
                savefile: str='p_post_stretch.pdf',
                dpi: int=300,) -> tuple:
    """
    Fit quadratic energy-stretch curve, extract equilibrium values, and plot energy and c/a ratio.

    Args:
        jobn (list): List of stretch factors or job identifiers.
        Etot (list): Total energies corresponding to jobn (eV).
        natoms (int): Number of atoms in the system.
        stretch_type (str): Type of stretch ('x', 'y', 'z', 'xy', etc.).
        rvectors_ref (array): Reference row vectors of the lattice.
        lca (list): c/a ratios for each stretch.
        if_save (bool): Whether to save the plot. Default False.
        savefile (str): Filename for saved plot. Default 'p_post_stretch.pdf'.
        dpi (int): DPI for saved plot. Default 300.

    Returns:
        tuple: 
            - (fig, axes): matplotlib Figure and Axes objects.
            - (coeffs, extr_x, extr_y, extr_rvectors): Fitted quadratic coefficients, equilibrium factor, energy, and scaled row vectors.
    """
    fig, axes = my_plot(fig_subp=[2,1],fig_sharex =False)

    x = np.array([float(job) for job in jobn])
    y = np.array(Etot) / natoms  # eV per atom
    coeffs, _ = my_ployfit(x, y, deg=2)
    p = np.poly1d(coeffs)
    a, b, c = coeffs
    x_fit = np.linspace(min(x), max(x), 100)
    y_fit = p(x_fit)
    extr_x = -b / (2 * a)
    extr_y = p(extr_x)
    y = (y - extr_y) * 1000  # meV per atom
    y_fit = (y_fit - extr_y) * 1000  # meV per atom

    # row vector [a1, a2, a3]
    extr_rvectors = rvectors_ref * extr_x


    ax = axes[0]
    ax.plot(x, y, marker = 'o', linestyle = '')
    ax.plot(x_fit, y_fit)
    ax.set_xlabel('Stretch factor (-)')
    ax.set_ylabel(r'$E - E_0$ (meV per atom)')
    textstr = (
        f'Type: {stretch_type}\n' +
        fr'$a_0$={extr_rvectors},' + "\n" +
        fr'$E_0$={extr_y:.8f} eV/atom'
    )
    ax.text(0.5, 0.95, textstr, transform=ax.transAxes, ha='center', va='top', fontsize = 20)

    ax = axes[1]
    ax.plot(x, lca, marker = 'o', linestyle = '')
    ax.set_xlabel('Stretch factor (-)')
    ax.set_ylabel(r'$\dfrac{c}{a}$', rotation=0)
    # if any(np.array(lca) < 2) and any(np.array(lca) > 0):
    #     ax.set_ylim(0, 2)

    if if_save:
        plt.savefig(savefile, dpi=dpi)

    return (fig, axes) , (coeffs, extr_x, extr_y, extr_rvectors) 
