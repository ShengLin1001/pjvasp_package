"""
atominfo submodule

This submodule provides functions for plotting various atomic structure properties, including interlayer distances,
z-positions, and radial distribution functions (RDF). It also includes functions for reading and writing stretch analysis
results.

Functions:
    - my_plot_interlayer_distance: Calculate and plot interlayer distances.
    - my_plot_zpositions: Plot z-positions of atoms and display material thickness.
    - my_plot_rdf: Plot the radial distribution function (RDF) of an atomic structure.
"""

from ase import Atoms

import os

import numpy as np

import pandas as pd

import matplotlib.pyplot as plt

from mymetal.universal.plot.plot import my_plot
from mymetal.universal.atom.neighbor import get_neighbor_distances

def my_plot_interlayer_distance(atoms: Atoms= None, if_plot: bool = True, 
                                if_save: bool = True, 
                                if_save_txt: bool = True,
                                save_plot_path: str = './p_post_inter_distance.jpg', 
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
    if if_save_txt:
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

def my_plot_rdf( atoms: Atoms = None, bins: float = 50, cutoff: float = 10, title: str = '') -> tuple:
    
    distances, _ = get_neighbor_distances(atoms, cutoff=cutoff)

    fig, axes = my_plot()
    ax = axes
    ax.hist(distances, bins=bins, density=True)
    ax.set_xlim(0, cutoff)
    ax.set_xlabel('r (Å)')
    ax.set_ylabel('g(r)')
    ax.set_title(title)
    return fig, axes

