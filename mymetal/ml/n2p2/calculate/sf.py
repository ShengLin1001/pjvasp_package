"""
symmetry functions (SFs) submodule.

This module provides functions to calculate radial and angular symmetry functions (SFs) for atomic structures using ASE.

Classes:
    - mysfparams: Class to handle symmetry function parameters and calculations.

Methods:
    - cal_sf: Calculate SFs for a list of structures.
    - cal_g2: Compute radial (g2) SF for a single structure.
    - cal_g3_g9: Compute angular (g3 or g9) SF for a single structure.
    - cal_cutoff: Evaluate the cutoff function for a given distance.
    - write_settings_overview: Generate a summary of SF settings.
    - write_parameter_strings: Generate formatted SF definitions.
    - plot: Plot Gaussian functions for the SFs (supported 2, 3, 9).
    
Functions:
    - load_from_npp_data: Load SF parameters from an NNP input file.
    - get_radial_pairs: Generate unique pairs of elements for radial SFs.
    - get_angular_pairs: Generate unique triplets of elements for angular SFs.

Change Log:
    - 2025.10.16: Fixed a bug in cal_g3_g9 where rc parameter was not passed, and added rc offset to Gaussian peaks.
    - 2025.10.15: Updated to support multi-element systems, but doesn't test.
    - 2025.10.14: Initial implementation for single-element systems.
    - 2025.10.21: Move the plot_g2() funtion to class mysfparams, and add the load_from_nnp_data() function.

Note:
    - The cal_sf() method in class mysfparams is computationally expensive for angular SFs,
      as it involves triple nested loops and scales poorly with system size.
"""

import numpy as np
import pandas as pd
from ase import Atoms
from ase.neighborlist import NeighborList
from itertools import combinations_with_replacement, combinations
from mymetal.universal.plot.plot import my_plot
import matplotlib.pyplot as plt

# For one pair of (center, neighbor: list = [neighbor1] or [neighbor1, neighbor2])  atoms
class mysfparams:
    """Handle symmetry function parameters and compute SFs for atomic structures.

    Supports radial (g2) and angular (g3/g9) symmetry functions.

    Args:
        g (str): Type of symmetry function. Options: 'g2', 'g3', 'g9'.
        center_atoms (str): Symbol of central atom.
        neighbor_atoms (list): Symbols of neighbor atoms (1 for g2, 2 for g3/g9).
        leta (list, optional): List of eta parameters.
        lrs (list, optional): List of r_shift parameters.
        lrc (list, optional): List of cutoff radii.
        llambd (list, optional): Lambda parameters for g3/g9, lambda must be -1 or 1 ? .
        lzeta (list, optional): Zeta parameters for g3/g9.
        cutoff_type_id (int, optional): ID for cutoff function type. Default 6 (CT_POLY2).
        cutoff_alpha (float, optional): Inner radius scaling factor for cutoff.

    Attributes:
        g (str): Type of symmetry function.
        center_atoms (str): Symbol of central atom.
        neighbor_atoms (list): Symbols of neighbor atoms.
        leta (list): List of eta parameters.
        lrs (list): List of r_shift parameters.
        lrc (list): List of cutoff radii.
        llambd (list): Lambda parameters for g3/g9.
        lzeta (list): Zeta parameters for g3/g9.
        cutoff_type_id (int): ID for cutoff function type.
        cutoff_alpha (float): Inner radius scaling factor for cutoff.
    """

    def __init__(self, g: str = 'g3', center_atoms: str = 'Au', neighbor_atoms: list = ['Au', 'Au'],   # ['Au', 'Ag'] same ['Ag', 'Au']
                 leta: list = None, lrs: list = None, lrc: list = None, llambd: list = None, lzeta: list = None,
                 cutoff_type_id: int = 6, cutoff_alpha: float = 0.0,):

        # Check input
        if g not in ['g2', 'g3', 'g9']:
            raise ValueError('Not support this sf_type')
        if g == 'g2':
            if len(leta) != len(lrs) or len(leta) != len(lrc):
                raise ValueError('For g2, length of leta, lrs, lrc must be equal')
            if len(neighbor_atoms) != 1:
                raise ValueError('For g2, length of neighbor_atoms must be 1.')
        elif g == 'g3' or g == 'g9':
            if lrs is None or len(lrs) == 0:
                lrs = [0.0] * len(leta)
            if llambd is None or lzeta is None:
                raise ValueError('For g3/g9, llambd and lzeta must be provided')
                # lambda = 1, zeta = 1,4
                #llambd = [1.0] * len(leta)
                #lzeta = [1.0, 4.0] * (len(leta) // 2) + [1.0] * (len(leta) % 2)
            if len(leta) != len(lzeta) or len(leta) != len(llambd) or len(leta) != len(lrc):
                raise ValueError('For g3/g9, length of leta, lzeta, llambd, lrc must be equal')
            if len(neighbor_atoms) != 2:
                raise ValueError('For g3/g9, length of neighbor_atoms must be 2.')
            
        self.g = g
        self.center_atoms = center_atoms
        self.neighbor_atoms = neighbor_atoms
        self.leta = leta
        self.lrs = lrs
        self.lrc = lrc
        self.llambd = llambd
        self.lzeta = lzeta
        self.cutoff_type_id = cutoff_type_id
        self.cutoff_alpha = cutoff_alpha
    
    ## calculate Radial / Angular SF
    def cal_sf(self, latoms: list = None, if_ave_in_one_strurture: bool = False) -> list:
        """Calculate symmetry functions for a list of structures.

        Args:
            latoms (list of ase.Atoms): List of atomic structures.
            if_ave_in_one_strurture (bool): If True, average SFs over atoms in each structure.

        Returns:
            np.ndarray: SF array of shape (M, N) if averaged, else (total_atoms, N).
        """
        # M * N
        # M - number of structures
        # N - number of SFs
        leta = self.leta
        lrs = self.lrs
        lrc = self.lrc
        lzeta = self.lzeta
        llambda = self.llambd

        sfs = []
        m = len(latoms)
        n = len(leta)
        m_if_not_ave = 0

        # LOOP ON M
        for atoms in latoms:
            # LOOP ON N
            sf_in_one_structure = []
            m_if_not_ave += len(atoms)

            if self.g == 'g2':
                for eta, rs, rc in zip(leta, lrs, lrc):
                    # append sf list to sf_in_one_structure
                    # x atom in structure, 1 sf: return 1*x list
                    sf_in_one_structure.append(self.cal_g2(atomsc = atoms, eta = eta, rs = rs, rc = rc))
            elif self.g == 'g3' or self.g == 'g9':
                for eta, zeta, my_lambda, rc, rs in zip(leta, lzeta, llambda, lrc, lrs):
                    sf_in_one_structure.append(self.cal_g3_g9(atomsc = atoms, eta = eta, zeta = zeta, my_lambda = my_lambda, 
                                                              rc = rc,rs = rs, if_wide = False if self.g == 'g3' else True))
            else:
                raise ValueError('Not support this sf_type')
            
            # sf_in_one_structure: N * x
            sf_in_one_structure = np.array(sf_in_one_structure).T   # x * N
            if if_ave_in_one_strurture:
                sf_in_one_structure = np.mean(sf_in_one_structure, axis=0)  # 1 * N
            # append list sf_in_one_structure to sfs
            sfs.append(sf_in_one_structure)

        # Check shape of sfs
        sfs = np.vstack(sfs)  # M * N if if_ave_in_one_strurture else ( M * x ) rol * N, x is not same in different structure

        if if_ave_in_one_strurture:
            desired_shape = (m, n)
        else:
            desired_shape = (m_if_not_ave, n)

        if sfs.shape != desired_shape:
            raise ValueError('Shape of sfs not equal to (m, n)')
        
        return sfs

    # For single structure, return a g2 for specified (eta, rs, rc)
    def cal_g2(self, atomsc: Atoms = None, eta: float = None, rs: float = None, rc: float =   None,
                ) -> list:
        """Compute radial (g2) SF for a single structure and given parameters.

        Args:
            atomsc (ase.Atoms): Atomic structure.
            eta (float): Gaussian width parameter.
            rs (float): Radial shift.
            rc (float): Cutoff radius.

        Returns:
            list: SF values for all atoms in the structure.
        """

        atoms = atomsc.copy()
        cutoffs = [rc / 2] * len(atoms)  # NeighborList needs each atom's radius
        nl = NeighborList(cutoffs, self_interaction=False, bothways=True, skin=0.0)
        nl.update(atoms)

        # LOOP over every center atom
        lsf = [0.0] * len(atoms)
        for i in range(len(atoms)):

            if atoms[i].symbol != self.center_atoms:
                continue

            indices, offsets = nl.get_neighbors(i)   # offset - shift the periodic image
            num_neighbors = len(indices)

            # LOOP over every neighbor j of center atom i
            for j, offset in zip(indices, offsets):
                
                if atoms[j].symbol not in self.neighbor_atoms:
                    continue

                r_ij_vec = atoms.positions[j] + np.dot(offset, atoms.get_cell()) - atoms.positions[i]
                r_ij = np.linalg.norm(r_ij_vec)

                # Check r_ij
                if r_ij > rc:
                    raise ValueError('r_ij > rc')
                if r_ij == 0:
                    raise ValueError('r_ij == 0')

                lsf[i] += np.exp(-eta * (r_ij - rs) ** 2)\
                                        * self.cal_cutoff( rc, r_ij)

        # list
        return lsf

    # For single structure, return a g3/g9 for specified (eta, zeta, lambda, rc)
    def cal_g3_g9(self, atomsc: Atoms = None, eta: float = None, zeta: float = None, my_lambda: float = None, rc: float =   None, rs: float = None,
                 if_wide: bool = False) -> list:
        """Compute angular (g3 or g9) SF for a single structure and given parameters.

        Args:
            atomsc (ase.Atoms): Atomic structure.
            eta (float): Gaussian width parameter.
            zeta (float): Angular exponent.
            my_lambda (float): Angular parameter.
            rc (float): Cutoff radius.
            rs (float): Radial shift.
            if_wide (bool): If True, compute g9; else compute g3.

        Returns:
            list: SF values for all atoms in the structure.
        """
        atoms = atomsc.copy()
        
        cutoffs = [rc / 2] * len(atoms)  # NeighborList needs each atom's radius
        nl = NeighborList(cutoffs, self_interaction=False, bothways=True, skin=0.0)
        nl.update(atoms)

        # LOOP over every center atom
        lsf = [0.0] * len(atoms)

        lnneigh = [0.0] * len(atoms)

        # i - center atom
        for i in range(len(atoms)):

            if atoms[i].symbol != self.center_atoms:
                continue

            indices, offsets = nl.get_neighbors(i)   # offset - shift the periodic image
            num_neighbors = len(indices)
            lnneigh[i] = num_neighbors
            
            # LOOP over every pair of neighbors (i, j, k), while j < k to avoid double counting
            # Ensure for the same center atom, the same (i,j,k) is only calculated once
            neighbor_j_k_list = list(combinations(list(range(num_neighbors)), 2))  # 0 to num_neighbors-1, all (j,k) pairs

            if len(neighbor_j_k_list) != (num_neighbors * (num_neighbors - 1)) // 2:
                raise ValueError('len(j_k_list) not equal to (num_neighbors * (num_neighbors - 1)) // 2')

            for neighbor_j, neighbor_k in neighbor_j_k_list:
                # j from 0 to num_neighbors-1
                # k from j+1 to num_neighbors-1
                # 首先一个原子对只被考虑了一次，另外 neighbor_atoms 中的顺序不影响结果
                # 例如 ['Au', 'Ag'] 和 ['Ag', 'Au'] 是一样的
                j = indices[neighbor_j]
                k = indices[neighbor_k]
                offsets_j = offsets[neighbor_j]
                offsets_k = offsets[neighbor_k]

                if atoms[j].symbol not in self.neighbor_atoms or atoms[k].symbol not in self.neighbor_atoms:
                    continue

                r_ij_vec = atoms.positions[j] + np.dot(offsets_j, atoms.get_cell()) - atoms.positions[i]
                r_ij = np.linalg.norm(r_ij_vec)
                    
                r_ik_vec = atoms[k].position + np.dot(offsets_k, atoms.get_cell()) - atoms.positions[i]
                r_ik = np.linalg.norm(r_ik_vec)
                
                r_jk_vec = r_ik_vec - r_ij_vec
                r_jk = np.linalg.norm(r_jk_vec)

                # Check r_ij, r_ik
                if r_ij > rc or r_ik > rc:
                    raise ValueError('r_ij or r_ik > rc')
                if r_ij == 0 or r_ik == 0 or r_jk == 0:
                    raise ValueError('r_ij or r_ik or r_jk == 0')
                
                cos_theta = np.dot(r_ij_vec, r_ik_vec) / (r_ij * r_ik)

                # important: minus rs, not plus
                # Shift
                r_ijs = r_ij - rs
                r_iks = r_ik - rs
                r_jks = r_jk - rs

                if if_wide:
                    # g9, wide angular
                    cutoff_product =   self.cal_cutoff(rc, r_ij)\
                                        * self.cal_cutoff(rc, r_ik)\
                    
                    gaussian = np.exp(-eta * (r_ijs**2 + r_iks**2))

                else:
                    # g3, narrow angular
                    cutoff_product =   self.cal_cutoff(rc, r_ij)\
                                        * self.cal_cutoff(rc, r_ik)\
                                        * self.cal_cutoff(rc, r_jk)
                                
                    gaussian = np.exp(-eta * (r_ijs**2 + r_iks**2 + r_jks**2))

                scale_factor = 2 ** (1 - zeta)
                lambda_item = (1 + my_lambda * cos_theta) ** zeta
                lsf[i] +=  scale_factor * lambda_item * gaussian * cutoff_product

        return lsf

    def cal_cutoff(self, rc: float = None, r: float = None) -> float:
        """Evaluate the cutoff function for a given distance.

        Args:
            rc (float): Cutoff radius.
            r (float): Distance between atoms.

        Returns:
            float: Cutoff function value in [0,1].
        """
        type_id = self.cutoff_type_id
        alpha = self.cutoff_alpha

        cutoff_type_map = {
            0: "CT_HARD",
            1: "CT_COS",
            2: "CT_TANHU",
            3: "CT_TANH",
            4: "CT_EXP",
            5: "CT_POLY1",
            6: "CT_POLY2",
            7: "CT_POLY3",
            8: "CT_POLY4"
        }

        if type_id not in cutoff_type_map:
            raise ValueError(f'Not support cutoff function type: {type_id}')
        
        if rc is None or r is None:
            raise ValueError('rc and r must be provided')
        
        type_str = cutoff_type_map[type_id]

        # inner radius
        rci = alpha * rc
        
        if r >= rc:
            return 0.0
        elif r >= 0 and r < rci:
            return 1.0
        elif r >= rci and r < rc:
            x = (r - rci) / (rc - rci)

            if type_str == 'CT_HARD':
                return 1.0
            
            elif type_str == 'CT_COS':
                # f(x) = 0.5 [cos(pi*x) + 1]
                return 0.5 * (np.cos(np.pi * (r - rci) / (rc - rci)) + 1)
            
            elif type_str == 'CT_TANHU':
                # f_c(r) = tanh^3(1 - r/rc)
                return np.tanh(1 - r / rc) ** 3
            
            elif type_str == 'CT_TANH':
                # f_c(r) = c * tanh^3(1 - r/rc), 且 f(0)=1
                # 由 f(0)=1 → c = 1 / tanh^3(1)
                c = 1 / (np.tanh(1) ** 3)
                return c * np.tanh(1 - r / rc) ** 3

            elif type_str == 'CT_EXP':
                # f(x) = exp(1 - 1/(1 - x^2))
                return np.exp(1 - 1 / (1 - x ** 2))

            elif type_str == 'CT_POLY1':
                # f(x) = (2x - 3)x^2 + 1
                return (2 * x - 3) * x ** 2 + 1

            elif type_str == 'CT_POLY2':
                # f(x) = ((15 - 6x)x - 10)x^3 + 1
                return ((15 - 6 * x) * x - 10) * x ** 3 + 1

            elif type_str == 'CT_POLY3':
                # f(x) = (x((20x - 70)x + 84) - 35)x^4 + 1
                return (x * ((20 * x - 70) * x + 84) - 35) * x ** 4 + 1

            elif type_str == 'CT_POLY4':
                # f(x) = (x(x((315 - 70x)x - 540) + 420) - 126)x^5 + 1
                return (x * (x * ((315 - 70 * x) * x - 540) + 420) - 126) * x ** 5 + 1

            else:
                raise ValueError(f'Not support cutoff function type: {type_str}')

    def write_settings_overview(self):
        content = ''
        if self.g == 'g2':
            content += '#' * 100 + '\n'
            content += f'# Radial symmetry functions set, for center atom: {self.center_atoms}, neighbor atom: {self.neighbor_atoms[0]}\n'
            content += '#' * 100 + '\n'
            header = f'{"#  id":>5} {"eta":>12} {"rs":>12} {"rc":>12}\n'
            content += header

            for i, (eta, rs, rc) in enumerate(zip(self.leta, self.lrs, self.lrc)):
                content += f'#{i:4d} {eta:12.4e} {rs:12.4e} {rc:12.4e}\n'
        
        if self.g == 'g3' or self.g == 'g9':
            if self.g == 'g3':
                g_type = 'Narrow angular symmetry functions'
            else:
                g_type = 'Wide angular symmetry functions'
            content += '#' * 100 + '\n'
            content += f'# {g_type} set, for center atom: {self.center_atoms}, neighbor atoms: {self.neighbor_atoms[0]}, {self.neighbor_atoms[1]}\n'
            content += '#' * 100 + '\n'
            header = f'{"#  id":>5} {"eta":>12} {"lambda":>12} {"zeta":>12} {"rc":>12} {"rs":>12}\n'
            content += header

            for i, (eta, my_lambda, zeta, rc, rs) in enumerate(zip(self.leta, self.llambd, self.lzeta, self.lrc, self.lrs)):
                content += f'#{i:4d} {eta:12.4e} {my_lambda:12.4e} {zeta:12.4e} {rc:12.4e} {rs:12.4e}\n'
        return content

    def write_parameter_strings(self):
        """Generate a formatted string describing the symmetry functions.

        Returns:
            str: Formatted SF definitions for the given parameters.
        """
        content = ''
        if self.g=='g2':
            for eta, rs, rc in zip(self.leta, self.lrs, self.lrc):
                # <element-central> 2 <element-neighbor> <eta> <rshift> <rcutoff>
                content += 'symfunction_short %3s 2 %3s %12.4e %12.4e %12.4e\n'%(self.center_atoms,self.neighbor_atoms[0],\
                                                                eta,rs,rc)

        if self.g=='g3' or self.g=='g9':
            if self.g=='g3':
                type_id = '3'
            else:
                type_id = '9'
            for eta, zeta, my_lambda, rc, rs in zip(self.leta, self.lzeta, self.llambd, self.lrc, self.lrs):
                    # <element-central> 3 <element-neighbor1> <element-neighbor2> <eta> <lambda> <zeta> <rcutoff> <<rshift>>
                    # <element-central> 9 <element-neighbor1> <element-neighbor2> <eta> <lambda> <zeta> <rcutoff> <<rshift>>
                    content += 'symfunction_short %3s %s %3s %3s %12.4e %3d %12.4e %12.4e %12.4e\n'%(self.center_atoms,type_id,self.neighbor_atoms[0],self.neighbor_atoms[1],\
                                                            eta, my_lambda, zeta, rc, rs)
                    
        return content
                
    # Gaussian functions plotting
    def plot(   self,
                fig_subp: tuple = None,
                if_save: bool = False, save_path: str = 'symmetry_functions.pdf') -> tuple:

        if self.g == 'g2':
            parameters = [self.leta, self.lrs, self.lrc]
        elif self.g == 'g3' or self.g == 'g9':
            parameters = [self.leta, self.lrc, self.lrs, self.llambd, self.lzeta]
        else:
            raise ValueError(f'Not support this sf_type: {self.g}')

        if fig_subp == None:
            fig_subp = [2, len(self.leta) // 2 + len(self.leta) % 2]

        fig, axes = my_plot(fig_subp=fig_subp, fig_sharex=False)
        i = 0
        j = 0
        index = 0

        for params in zip(*parameters):

            ax = axes[i, j]

            if self.g == 'g2':
                eta, rs, rc = params
                sigma = np.sqrt(1/(2*eta))
                # Input: eta, rs, rc, r_ij
                r = np.linspace(0, rc, 500)
                y = np.exp(-eta * (r - rs) ** 2)
                    #                    * self.cal_cutoff( rc, r)
                xlabel = r'$r$ (Å)'
                ylabel = r'$g_2$'
                xlim = (0, rc)
                ylim = (-0.05, 1.1)
                formula = r'$e^{-\eta (r - r_s)^2}$'
                x = r

                text = (
                    rf"{formula}" +"\n"
                    rf"$r_{{c}}$ = {rc}" +"\n"
                    rf"$r_{{s}}$ = {rs}" +"\n"
                    rf"$\eta$ = {eta}" +"\n"
                    rf"$r_c / \sigma$ = {rc/sigma:.3f}"
                )


            elif self.g == 'g3' or self.g == 'g9':
                eta, rc, rs, my_lambda, zeta = params
                sigma = np.sqrt(1/(2*eta))
                # Don't multiply cutoff here
                # Input: eta, zeta, my_lambda, rc, rs, r_ij, r_ik, r_jk, theta
                # Here only consider various theta.
                theta = np.linspace(0, 2 * np.pi, 500)  # 0 to pi
                y = 2 ** (1 - zeta) *\
                      (1 + my_lambda * np.cos(theta)) ** zeta
                    # * np.exp(-eta * ( (r_ij - rs) ** 2 + (r_ik - rs) ** 2 + (r_jk - rs) ** 2 ))\
                    # * self.cal_cutoff( rc, r_ij)\
                    # * self.cal_cutoff( rc, r_ik)\
                    # * self.cal_cutoff( rc, r_jk)
                xlabel = r'$\theta (degree)$'
                ylabel = r'$g_3$' if self.g == 'g3' else r'$g_9$'
                xlim = (0, 360)
                ylim = (-0.05, 2.1)
                formula = r'$2^{1-\eta}(1+\lambda\cos\theta_{ijk})^\zeta$'
                x = theta * 360 / (2 * np.pi)
                
                text = rf"{formula}" +"\n"     +\
                    rf"$r_{{c}}$ = {rc}" +"\n" +\
                    rf"$r_{{s}}$ = {rs}" +"\n" +\
                    rf"$\eta$ = {eta}" +"\n"   +\
                    rf"$r_c / \sigma$ = {rc/sigma:.3f}" +"\n" +\
                    rf"$\zeta$ = {zeta}" +"\n" +\
                    rf"$\lambda$ = {my_lambda}"

            ax.plot(x, y)
            ax.axvline(rs-sigma, linestyle='--')
            ax.axvline(rs+sigma, linestyle='--')
            ax.set_xlim(xlim[0], xlim[1])
            ax.set_ylim(ylim[0], ylim[1])
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            #print(ax.get_ylim())
            ax.text(0.95, 0.05, text, 
                    transform=ax.transAxes, verticalalignment='bottom', horizontalalignment='right')

            index = index + 1
            j =  index % fig_subp[1]
            i = index // fig_subp[1]

        if if_save:
            plt.savefig(save_path)

        return fig, axes


def load_from_npp_data(nnp_file: str = 'input.nn') -> tuple:
    """
    Args:
        nnp_file (str): Path to the NNP input file containing symmetry function definitions.
    
    Returns:
        tuple: Three dictionaries for g2, g3, and g9 symmetry functions.
            Each dictionary maps (center_atom, neighbor_atoms) to a DataFrame of parameters.
    
    Example:
        ```python
        g2sfs, g3sfs, g9sfs = load_from_npp_data('input.nn')

        for (centor, neighbors), df in g2sfs.items():
            temp_dict = {
                'center_atoms': df['lcenter_atoms'].iloc[0],
                'neighbor_atoms': df['lneighbor_atoms'].iloc[0],
                'leta': df['leta'].tolist(),
                'lrs': df['lrs'].tolist(),
                'lrc': df['lrc'].tolist()
            }
            g2sfs = mysfparams(g = 'g2', **temp_dict)

            g2sfs.plot()


        for (centor, neighbors), df in g3sfs.items():
            temp_dict = {
                'center_atoms':  df['lcenter_atoms'].iloc[0],
                'neighbor_atoms': df['lneighbor_atoms'].iloc[0],
                'leta': df['leta'].tolist(),
                'llambd': df['llambd'].tolist(),
                'lzeta': df['lzeta'].tolist(),
                'lrc': df['lrc'].tolist()
            }
            g3sfs = mysfparams(g = 'g3', **temp_dict)

            g3sfs.plot()
        ```
    """

    # For saving parameters
    g2 = {
        'lcenter_atoms': [],
        'lneighbor_atoms': [],
        'leta': [],
        'lrs': [],
        'lrc': []
    }

    g3 = {
        'lcenter_atoms': [],
        'lneighbor_atoms': [],
        'leta': [],
        'llambd': [],
        'lzeta': [],
        'lrc': [],
        'lrs': [],
    }

    g9 = {
        'lcenter_atoms': [],
        'lneighbor_atoms': [],
        'leta': [],
        'llambd': [],
        'lzeta': [],
        'lrc': [],
        'lrs': []

    }

    with open(nnp_file, 'r') as fd:
        lines = fd.readlines()

    for line in lines:
        if 'symfunction_short' in line:
            l = line.split()
            if l[2] == '2':
                # symfunction_short <element-central> 2 <element-neighbor> <eta> <rshift> <rcutoff>
                g2['lcenter_atoms'].append(l[1])
                g2['lneighbor_atoms'].append([l[3]])
                g2['leta'].append(float(l[4]))
                g2['lrs'].append(float(l[5]))
                g2['lrc'].append(float(l[6]))
            elif l[2] == '3':
                # symfunction_short <element-central> 3 <element-neighbor1> <element-neighbor2> <eta> <lambda> <zeta> <rcutoff> <<rshift>>
                g3['lcenter_atoms'].append(l[1])
                g3['lneighbor_atoms'].append([l[3], l[4]])
                g3['leta'].append(float(l[5]))
                g3['llambd'].append(float(l[6]))
                g3['lzeta'].append(float(l[7]))
                g3['lrc'].append(float(l[8]))
                # lrs is optional
                if len(l) == 10:
                    g3['lrs'].append(float(l[9]))
                else:
                    g3['lrs'].append(float(0.0))

            elif l[2] == '9':
                g9['lcenter_atoms'].append(l[1])
                g9['lneighbor_atoms'].append([l[3], l[4]])
                g9['leta'].append(float(l[5]))
                g9['llambd'].append(float(l[6]))
                g9['lzeta'].append(float(l[7]))
                g9['lrc'].append(float(l[8]))
                # lrs is optional
                if len(l) == 10:
                    g9['lrs'].append(float(l[9]))
                else:
                    g9['lrs'].append(float(0.0))

    # Group parameters by center atom and neighbor atoms
    lg_groups = []
    for g in [g2, g3, g9]:
        df = pd.DataFrame(g)
        df['lneighbor_atoms'] = df['lneighbor_atoms'].apply(tuple)  # Convert lists to tuples for grouping
        # iterator of (group_name, group_df)
        dfgroups = df.groupby(['lcenter_atoms', 'lneighbor_atoms'])
        # tuple() to convert iterator to tuple, (((lcentor, lneighbor), sub_df1),...), then dict() to convert to dict
        g_groups = dict(tuple(dfgroups))

        # Check correctness
        for (centor, neighbors), sub_df in g_groups.items():

            center_col = sub_df['lcenter_atoms']
            neigh_col = sub_df['lneighbor_atoms'].apply(tuple)

            if center_col.nunique() != 1 or center_col.iloc[0] != str(centor):
                raise ValueError(f"Center atom mismatch! Expected {centor}, got {center_col.unique()}")

            if neigh_col.nunique() != 1 or neigh_col.iloc[0] != tuple(neighbors):
                raise ValueError(f"Neighbor atoms mismatch! Expected {neighbors}, got {neigh_col.unique()}")

        lg_groups.append(g_groups)

    #      g2,           g3,         , g9
    return lg_groups[0], lg_groups[1], lg_groups[2]

def get_radial_pairs(elements: list = None) -> list:
    """Generate unique radial atom pairs [center, neighbor] for symmetry functions.

    Args:
        elements (list): List of element symbols.

    Returns:
        list of tuple: Sorted list of unique pairs (center, neighbor) including self-pairs.

    Note:
        ['Cu', 'S'] and ['S', 'Cu'] both needs to apply for SFs.
    """
    uni_elements = list(set(elements))
    return sorted(list(combinations_with_replacement(uni_elements, 2)))

def get_angular_pairs(elements: list = None) -> list:
    """Generate unique angular atom triplets for symmetry functions.

    Each triplet is [center, neighbor1, neighbor2] with neighbor1 <= neighbor2.

    Args:
        elements (list): List of element symbols.

    Returns:
        list of list: Sorted list of triplets covering all unique center-neighbor combinations.
    """
    uni_elements = list(set(elements))
    num_ele = len(uni_elements)
    triplets = []
    # [comb[0], comb[1]]  -> [ele, comb[0], comb[1]]
    count = 0
    for comb in combinations_with_replacement(uni_elements, 2):
        for ele in uni_elements:
            triplets.append([ele, comb[0], comb[1]])
            count += 1
    if count != num_ele * (num_ele * (num_ele + 1) / 2):
        raise ValueError('count not equal to num_ele * (num_ele * (num_ele + 1) / 2)')
    return sorted(triplets)
