"""
symmetry functions (SFs) submodule.

This module provides functions to calculate radial and angular symmetry functions (SFs) for atomic structures using ASE.

Classes:
    - mysfparams: Class to handle symmetry function parameters and calculations.

Functions:
    - get_radial_pairs: Generate unique pairs of elements for radial SFs.
    - get_angular_pairs: Generate unique triplets of elements for angular SFs.

Change Log:
    - 2025.10.16: Fixed a bug in cal_g3_g9 where rc parameter was not passed, and added rc offset to Gaussian peaks.
    - 2025.10.15: Updated to support multi-element systems, but doesn't test.
    - 2025.10.14: Initial implementation for single-element systems.

Note:
    - The cal_sf() method in class mysfparams is computationally expensive for angular SFs,
      as it involves triple nested loops and scales poorly with system size.
"""

import numpy as np
from ase import Atoms
from ase.neighborlist import NeighborList
from itertools import combinations_with_replacement, combinations


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
            if lrs is None:
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
