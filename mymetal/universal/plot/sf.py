"""
Module for plotting sf functions.

This module provides functions to visualize Gaussian functions used in radial symmetry function construction.

Functions:
    - my_plot_gaussian_functions: Plot Gaussian functions with specified parameters.
    - my_plot_infos: Plot parameter relationships and derived quantities of Gaussian functions.

"""


import numpy as np
from mymetal.universal.plot.plot import my_plot
import matplotlib.pyplot as plt
import pandas as pd

def get_gaussian(eta: float = None, r_shift: float = None, r: np.array = None) -> np.array:
    return np.exp(-eta * (r - r_shift)**2)

def my_plot_gaussian_functions( leta: list = None,
                    lrs: list = None,
                    lrc: list = None,
                    fig_subp: tuple = None,
                    if_save: bool = True, save_path: str = 'gaussian_functions.pdf'):
    """
    Plot Gaussian functions used in radial symmetry function construction.

    Args:
        leta (list): List of η parameters.
        lrs (list): List of rₛ (center) parameters.
        lrc (list): List of r_c (cutoff) parameters.
        fig_subp (tuple): Subplot grid dimensions, e.g., (n_rows, n_cols).
        if_save (bool): Whether to save the plot to a file.
        save_path (str): File path to save the figure.

    Returns:
        None: Displays and optionally saves Gaussian function plots.

    Notes:
        - Each subplot shows a Gaussian function f(r) = exp(-η(r - rₛ)²).
        - Dashed lines indicate ±σ from the center position.
    """
    if fig_subp is None:
        raise ValueError('fig_subp must be provided, e.g., (2,3) for 2 rows and 3 columns.')

    lsigma = np.sqrt(1/(2*leta))

    fig, axes = my_plot(fig_subp=fig_subp, fig_sharex=False)
    i = 0
    j = 0
    index = 0

    for eta, r_shift, rcut, sigma in zip(leta, lrs, lrc, lsigma):
        #print(i, j)
        ax = axes[i, j]

        r = np.linspace(0, rcut, 500)
        f_r = get_gaussian(eta=eta, r_shift=r_shift, r=r)
        ax.plot(r, f_r, label=f'eta={eta:.4f}, rs={r_shift:.2f}')
        ax.axvline(r_shift-sigma, linestyle='--')
        ax.axvline(r_shift+sigma, linestyle='--')

        f_r_sigma1 = get_gaussian(eta=eta, r_shift=r_shift, r=r_shift+sigma)
        f_r_sigma2 = get_gaussian(eta=eta, r_shift=r_shift, r=r_shift-sigma)
        # when r = r_s + r_c / 2
        f_r_temp1 = get_gaussian(eta=eta, r_shift=r_shift, r=r_shift + rcut/2)
        f_r_temp2 = get_gaussian(eta=eta, r_shift=r_shift, r=r_shift - rcut/2)
        # when r = r_c
        f_r_rc = get_gaussian(eta=eta, r_shift=r_shift, r=rcut)

        ax.set_xlim(0, rcut)
        ax.set_ylim(0, 1.1)
        ax.set_xlabel('r')
        ax.set_ylabel('f(r)')
        ax.set_title(rf'$e^{{-{eta} (r - {r_shift})^2}}$')
        ax.text(0.95, 0.05, f'$r_{{c}}$={rcut}\n\
                        $r_{{s}}$={r_shift}\n\
                        $\eta$={eta}\n\
                        $\sigma$={sigma:.4f}={sigma/rcut:.3f}$r_c$=1/{rcut/sigma:.3f}$r_c$\n\
                        $log_{{10}}\sigma$={np.log10(sigma):.3f}\n\
                        $ln\sigma$={np.log(sigma):.3f}\n\
                        $f_{{r=r_c}}$={f_r_rc:.3f}', 
                transform=ax.transAxes, verticalalignment='bottom', horizontalalignment='right')

        index = index + 1
        j =  index % fig_subp[1]
        i = index // fig_subp[1]

    if if_save:
        plt.savefig(save_path)



def my_plot_infos(  leta: list = None,
                    lrs: list = None,
                    lrc: list = None,):
    """
    Plot parameter relationships and derived quantities of Gaussian functions.

    Args:
        leta (list): List of η parameters.
        lrs (list): List of rₛ (center) parameters.
        lrc (list): List of r_c (cutoff) parameters.

    Returns:
        tuple: (
            (etas, r_shifts, rcuts, sigmas),
            (sigma_log, sigma_log10),
            (f_r_rc, f_r_rc_log, f_r_rc_log10),
            rc_sigma_ratio
        )
        where each component corresponds to computed and plotted quantities.

    Notes:
        - Displays ln(σ), log₁₀(σ), f(r=r_c), and their logarithmic transforms.
        - Helps analyze how η, rₛ, and r_c influence Gaussian widths and magnitudes.
    """

    df = {
        'eta': leta,
        'rs (Å)': lrs,
        'rc (Å)': lrc
    }
    df = pd.DataFrame(df)
    df_sorted = df.sort_values(by=["rc (Å)", "rs (Å)", "eta"], ascending=[True, True, True])
    etas = df_sorted['eta']
    r_shifts = df_sorted['rs (Å)']
    rcuts = df_sorted['rc (Å)']
    sigmas = np.sqrt(1/(2*leta))
    sigma_log = np.log(sigmas)
    sigma_log10 = np.log10(sigmas)
    f_r_rc = get_gaussian(eta=etas, r_shift=r_shifts, r=rcuts)
    f_r_rc_log = np.log(f_r_rc)
    f_r_rc_log10 = np.log10(f_r_rc)
    rc_sigma_ratio = rcuts / sigmas

    fig, axes = my_plot(fig_subp=(1,5), fig_sharex=False)

    # 对数变换
    ax = axes[0]
    ax.plot(range(len(etas)), sigma_log, marker='o', linestyle='-')
    ax.set_xlabel('Index (-)')
    ax.set_ylabel(r'$ln(\sigma)$')

    ax = axes[1]
    ax.plot(range(len(etas)), sigma_log10, marker='o', linestyle='-')
    ax.set_xlabel('Index (-)')
    ax.set_ylabel(r'$log_{10}(\sigma)$')

    # f_r=rc
    ax = axes[2]
    ax.plot(range(len(etas)), f_r_rc, marker='o', linestyle='-')
    ax.set_xlabel('Index (-)')
    ax.set_ylabel(r'$f_{r=r_c}$')

    ax = axes[3]
    ax.plot(range(len(etas)), f_r_rc_log, marker='o', linestyle='-')
    ax.set_xlabel('Index (-)')
    ax.set_ylabel(r'$ln(f_{r=r_c})$')

    ax = axes[4]
    ax.plot(range(len(etas)), f_r_rc_log10, marker='o', linestyle='-')
    ax.set_xlabel('Index (-)')
    ax.set_ylabel(r'$log_{10}(f_{r=r_c})$')

    return (etas, r_shifts, rcuts, sigmas), (sigma_log, sigma_log10), (f_r_rc, f_r_rc_log, f_r_rc_log10), rc_sigma_ratio


