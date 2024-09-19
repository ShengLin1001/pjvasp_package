from numpy import array, ndarray, argsort
from ase.build import add_adsorbate
from ase import Atoms
import numpy as np
from ase.build import add_adsorbate
from numpy import argsort
from typing import List, Tuple, Optional

# def add_hydro_atoms(atoms:  Atoms = None,
#                     add_symbol: chr = 'H', 
#                     added_symbol: chr = 'O', 
#                     surf_range: ndarray = array([0.1, 11.1]), 
#                     shift_distance: ndarray = array([0.0, 0.0, 2.5]),
#                     surf_direction: int = 2) -> Atoms:
#     """
#     Add hydrogen atoms to the surface of the atoms object.
#     """
#     atoms = atoms.copy()
#     #collect = atoms.copy()
#     for atom in atoms:
#         if atom.position[surf_direction] > surf_range[0] and atom.position[surf_direction] < surf_range[1] and atom.symbol == added_symbol:
#             #print(atom.index, atom.symbol, atom.position)
#             atoms.append(Atom(symbol=add_symbol, position=atom.position + shift_distance))
#     atoms = atoms[atoms.numbers.argsort()]
#     return atoms

def add_hydro_atoms(my_atoms:  Atoms = None,
                    add_symbol: chr = 'H', 
                    added_symbol: chr = 'O', 
                    surf_range: ndarray = array([0.1, 11.1]), 
                    shift_distance: ndarray = array([0.0, 0.0, 2.5]),
                    surf_direction: int = 2,
                    reverse_every: bool = False,
                    reverse_length: int = 1,
                    fix_directions: list = [2],
                    sort_direction: int = 0,
                    offset = None,
                    mol_index = 0,
                    manual: bool = False) -> Atoms:
    """
    Adds hydrogen atoms to the surface of the specified atomic structure. For crystal systems.
    
    Args:
        my_atoms (Atoms, optional): The atomic structure to which hydrogen atoms will be added.
        add_symbol (str, optional): The symbol of the atom to add, defaults to 'H' for hydrogen.
        added_symbol (str, optional): The symbol of the surface atom to which hydrogen atoms are added, defaults to 'O' for oxygen.
        surf_range (ndarray, optional): A 1D array specifying the range along the surface direction within which atoms are considered for hydrogen addition.
        shift_distance (ndarray, optional): A 3D vector specifying the shift in position to place the hydrogen atom relative to the surface atom.
        surf_direction (int, optional): The axis index (0, 1, or 2 corresponding to x, y, or z) along which the surface is oriented.
        reverse_every (bool, optional): A flag to reverse the shift direction after each addition, defaults to False.
        fix_directions (list, optional): A list of indices of the axes where the shift distance should not be reversed.
        sort_direction (int, optional): The axis index to sort the surface atoms by their coordinates.
        offset (optional): Offset for placing the added atom, can be used for complex surface structures.
        mol_index (int, optional): Molecular index to specify which molecule to modify, useful in multi-molecular systems.

    Returns:
        Atoms: The modified atomic structure with added hydrogen atoms.

    Raises:
        TypeError: If `my_atoms` is not an instance of `Atoms`.
        TypeError: If `add_symbol` or `added_symbol` are not instances of `str`.
        ValueError: If `surf_range` or `shift_distance` do not contain exactly two or three elements.
        ValueError: If `surf_direction` or `sort_direction` are not one of the valid indices (0, 1, or 2).

    Example usage:
        added_atoms = add_hydro_atoms(temp, add_symbol='H', added_symbol='O', surf_range=[16.7, 17], shift_distance=[0.9441063, 0.2422, 0], 
                                      reverse_every=True)

        added_atoms = add_hydro_atoms(temp, add_symbol='H', added_symbol='O', surf_range=[16.7, 17], shift_distance=[0, 0, 0.9584])

        added_atoms = add_hydro_atoms(added_atoms, add_symbol='H', added_symbol='O', surf_range=[7.4, 7.6], shift_distance=[0, 0, -0.9584])

    Example process:
        # INPUT the SiO2 model
        atoms = read_vasp('SiO2.poscar')
        print(atoms)
        
        # SLAB model
        temp = atoms.repeat([1, 1, 2])
        print(temp)

        # Remove atoms based on z-position
        del temp[[atom.index for atom in temp if atom.position[2] < 0.9 or atom.position[2] > 10.18]]
        temp = temp[temp.numbers.argsort()]
        temp.center(vacuum=5, axis=2)
        print(temp)

        # Add hydrogen atoms
        added_atoms = add_hydro_atoms(temp, add_symbol='H', added_symbol='O', surf_range=[14.2, 15], shift_distance=[0, 0, 0.9584])
        added_atoms = add_hydro_atoms(added_atoms, add_symbol='H', added_symbol='O', surf_range=[4.9, 5.1], shift_distance=[0, 0, -0.9584])
        added_atoms.center(vacuum=5, axis=2)
        print(added_atoms)
        view(added_atoms)
    
        added_atoms = my_add_hydro_atoms(temp, add_symbol = 'H', added_symbol = 'O', surf_range=[16.7, 17], shift_distance = [0.9441063, 0.2422, 0], 
                                 reverse_every=True)
        added_atoms = add_hydro_atoms(temp, add_symbol = 'H', added_symbol = 'O', surf_range=[16.7, 17], shift_distance = [0, 0, 0.9584])
        added_atoms = add_hydro_atoms(added_atoms, add_symbol = 'H', added_symbol = 'O', surf_range=[7.4, 7.6], shift_distance = [0, 0, -0.9584])
    """
    if not isinstance(my_atoms, Atoms):
        raise TypeError("Expected my_atoms to be an instance of Atoms")
    if not isinstance(add_symbol, str) or not isinstance(added_symbol, str):
        raise TypeError("Expected add_symbol, added_symbol to be the instance of chr")
    if len(surf_range) != 2:
        raise ValueError("surf_range must contain exactly two elements")
    if len(shift_distance) != 3:
        raise ValueError("shift_distance must contain exactly three elements")
    if surf_direction not in [0, 1, 2] or sort_direction not in [0, 1, 2]:
        raise ValueError("surf_direction and sort_direction must be 0, 1, or 2")
    #print(surf_range)
    surf_range = array(surf_range)
    shift_distance = array(shift_distance)
    atoms = my_atoms.copy()
    atoms.set_pbc(True)
    shift = shift_distance
    collect = []
    # test top position
    positions = atoms.get_positions()
    position = positions[:,fix_directions]
    position.sort()
    max_position = float(max(position))

    for atom in atoms:
        if atom.position[surf_direction] > surf_range[0] and atom.position[surf_direction] < surf_range[1] and atom.symbol == added_symbol:
            collect.append((atom.index, atom.position))

    sorted_indices = argsort([item[sort_direction] for index,item in collect])
    sorted_collect = [collect[i] for i in sorted_indices]
    #print(sorted_collect)
    count = 0
    for index,item in sorted_collect:
        item = array(item)
        count = count + 1
        if manual:
            # Copy atoms for visualization and highlight the specific atom
            #view(atoms)
            # Manual input for direction change
            # Request user input for direction modifiers
            direction_modifiers = input(f"For this adsorbate of {index}\nEnter direction modifiers as a list (e.g., [1, 1, -1]): ")
            if direction_modifiers.strip() == '':
                direction_modifiers = [1, 1, 1]
            else:
                direction_modifiers = eval(direction_modifiers)  # Converts input string to a list

            shift = shift_distance * array(direction_modifiers)
            add_adsorbate(slab = atoms, adsorbate = add_symbol, height = shift[2], position = tuple(shift[:2] +item[:2]),
                        offset = offset, mol_index = mol_index)
        else:
            #print(shift[2])
            #print(tuple(shift[:2] +item[:2]))
            add_adsorbate(slab = atoms, adsorbate = add_symbol, height = item[2] - max_position + shift[2], position = tuple(shift[:2] +item[:2]),
                        offset = offset, mol_index = mol_index)
            if reverse_every and count % reverse_length == 0:
                shift = -1. * shift
                for fix_direction in fix_directions:
                    shift[fix_direction] = shift_distance[fix_direction]
        #atoms.append(Atom(symbol=add_symbol, position=atom.position + shift))
    atoms = atoms[atoms.numbers.argsort()]
    atoms.wrap()
    return atoms

def get_neighbors(atom_index: int = None, 
                  positions: ndarray = None, 
                  cell: ndarray = None,
                  pbc: List[bool] = None,
                  cutoff: float = None) -> Tuple[List[int], List[ndarray], List[float]]:
    """
    Get the list of neighbors for a given atom, along with their relative displacements 
    and distances, considering periodic boundary conditions (PBC).

    Args:
        atom_index (int): The index of the atom to find neighbors for.
        positions (ndarray): Array of atom positions, shape should be (N, 3).
        cell (ndarray): 3x3 matrix representing the unit cell.
        pbc (list of bool): A list of 3 booleans indicating whether PBC is applied in each 
            of the x, y, z directions.
        cutoff (float): The cutoff distance within which atoms are considered neighbors.

    Returns:
        tuple:
            neighbors (list of int): List of indices of neighboring atoms.
            offsets (list of ndarray): List of displacement vectors from the central atom 
                to each neighbor.
            distances (list of float): List of distances from the central atom to each neighbor.

    Raises:
        ValueError: If `positions`, `cell`, or `pbc` have incorrect dimensions, or if the cutoff 
            is not positive.
        IndexError: If `atom_index` is out of bounds.
    """
    # Input validation
    if not isinstance(atom_index, int) or atom_index < 0:
        raise TypeError(f"atom_index must be a non-negative integer, but got {atom_index}.")
    
    if  positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError(f"positions must be a numpy array of shape (N, 3), but got {positions.shape}.")
    
    if  cell.shape != (3, 3):
        raise ValueError(f"cell must be a numpy array of shape (3, 3), but got {cell.shape}.")
    
    if  len(pbc) != 3 or not all(isinstance(b, bool) for b in pbc):
        raise ValueError(f"pbc must be a list or array of 3 booleans, but got {pbc}.")
    
    if  cutoff <= 0:
        raise ValueError(f"cutoff must be a positive number, but got {cutoff}.")

    if atom_index >= len(positions):
        raise IndexError(f"atom_index {atom_index} is out of bounds for positions array of length {len(positions)}.")


    neighbors = []
    offsets = []
    distances = []
    
    # 获取当前原子的坐标
    ref_pos = positions[atom_index]

    # 获取晶胞的逆矩阵，用于处理周期性边界条件
    inv_cell = np.linalg.inv(cell)
    
    # 遍历其他所有原子
    for i, pos in enumerate(positions):
        if i == atom_index:
            continue  # 跳过自身
        
        # 计算原子间的位移向量
        displacement = pos - ref_pos
        
        # 考虑周期性边界条件
        for j in range(3):  # x, y, z 方向
            if pbc[j]:
                displacement[j] -= np.round(displacement[j] / cell[j, j]) * cell[j, j]
        
        # 将位移向量映射回原胞内
        displacement_frac = np.dot(displacement, inv_cell)
        displacement_cart = np.dot(displacement_frac, cell)
        
        # 计算原子间距离
        distance = np.linalg.norm(displacement_cart)
        
        # 如果距离小于截断半径，则该原子是邻居
        if distance < cutoff:
            neighbors.append(i)
            offsets.append(displacement_cart)  # 保存邻居的位移向量
            distances.append(distance)  # 保存邻居距离
    
    return neighbors, offsets, distances



def find_matching_atom_in_bulk(slab_atom_position: np.ndarray = None, 
                               bulk_positions: np.ndarray = None, 
                               tolerance: float = 1e-5) -> Optional[int]:
    """
    Finds the corresponding atom in the bulk based on the position of an atom in the slab.

    Args:
        slab_atom_position (np.ndarray): Coordinates of the atom in the slab.
        bulk_positions (np.ndarray): Coordinates of all atoms in the bulk.
        tolerance (float, optional): The tolerance for position comparison, default is 1e-5.

    Returns:
        Optional[int]: Index of the matching atom in the bulk. If no match is found, returns -1.

    Raises:
        ValueError: If `slab_atom_position` does not have a shape of (3,) or if `bulk_positions` 
                    is not a 2D array with shape (N, 3).
    """
    # Input validation
    if slab_atom_position.shape != (3,):
        raise ValueError(f"slab_atom_position must have shape (3,), but got shape {slab_atom_position.shape}.")
    
    if bulk_positions.ndim != 2 or bulk_positions.shape[1] != 3:
        raise ValueError(f"bulk_positions must have shape (N, 3), but got shape {bulk_positions.shape}.")

    # Iterate through all atoms in bulk to find the matching one
    for i, bulk_atom_position in enumerate(bulk_positions):
        # Compare the positions using np.allclose with the specified tolerance
        if np.allclose(slab_atom_position, bulk_atom_position, atol=tolerance):
            return i  # Return the index of the matching atom
    
    return -1  # Return -1 if no match is found

def find_unsaturated_neighbors(offset_bulk: list = None,
                               offset_slab: list = None,
                               tolerance: float=1e-3) -> list:
    """
    Identifies missing neighbors by comparing displacement vectors between bulk 
    and slab neighbors, indicating dangling bonds in the slab.

    Args:
        offset_bulk (list): List of displacement vectors for neighbors in the bulk.
        offset_slab (list): List of displacement vectors for neighbors in the slab.
        tolerance (float): A tolerance value for comparing displacement vectors, default is 1e-3.

    Returns:
        list: A list of displacement vectors that represent missing neighbors 
              in the slab, which exist in the bulk.
    """
    unsaturated_neighbors = []
    
    for bulk_offset in offset_bulk:
        found = False
        for slab_offset in offset_slab:
            # 计算两个位移向量之间的差异
            difference = np.linalg.norm(np.array(bulk_offset) - np.array(slab_offset))
            
            # 如果差异在容差范围内，认为是相同的邻居
            if difference < tolerance:
                found = True
                break
        
        # 如果在 slab 中找不到对应的邻居，则认为是缺失的
        if not found:
            unsaturated_neighbors.append(bulk_offset)
    
    return unsaturated_neighbors

def passivate_surface_custom(bulk: Atoms = None, 
                             slab: Atoms = None, 
                             adsorbates: dict = {'Si':'O','O':'H'},
                             cutoff: float = 1.9,
                             weights: dict = {'Si': 1.0,'O': 0.5},
                             slab_pbc: list = [True, True, True],
                             bulk_pbc: list = [True, True, True],
                             vacuum: float = 7.5,
                             axis: int = 2,
                             if_print: bool = True) -> Atoms:
    """
    Adds adsorbates (e.g., hydrogen) to passivate dangling bonds on the surface of a slab.
    
    Args:
        bulk (Atoms): The bulk crystal structure as an ASE Atoms object.
        slab (Atoms): The slab structure, cut from the bulk, as an ASE Atoms object.
        adsorbates (dict): A dictionary specifying the adsorbates to add. The key is the surface 
            atom symbol (e.g., 'Si'), and the value is the adsorbate symbol (e.g., 'H' for hydrogen).
        cutoff (float): The cutoff distance used for neighbor detection (default is 1.9).
        weights (dict): A dictionary specifying the weight factors for positioning the adsorbates.
            For example, {'Si': 1.0, 'O': 0.5}.
        slab_pbc (list of bool): A list of booleans specifying whether periodic boundary conditions
            should be applied along each axis for the slab (e.g., [True, True, True]).
        bulk_pbc (list of bool): A list of booleans specifying whether periodic boundary conditions
            should be applied along each axis for the bulk (e.g., [True, True, True]).
        vacuum (float): The amount of vacuum to add to the slab in the specified axis direction.
        axis (int): The axis (0, 1, or 2) along which vacuum should be added (default is 2, for the z-axis).
        if_print (bool): Control if print the amount of the dangling bonds.
        enlarge_size: 

    Returns:
        Atoms: The passivated slab structure with adsorbates added.

    Raises:
        ValueError: If the `bulk` or `slab` is not provided as an ASE Atoms object.
        ValueError: If `adsorbates`, `weights`, `slab_pbc`, or `bulk_pbc` are incorrectly formatted.
        ValueError: If `cutoff`, `vacuum`, or `axis` are not valid.
    """
    
    # Input validation
    if not isinstance(bulk, Atoms) or not isinstance(slab, Atoms):
        raise ValueError("Both `bulk` and `slab` must be instances of ASE Atoms.")
    
    if not isinstance(adsorbates, dict) or not all(isinstance(k, str) and isinstance(v, (str, Atoms)) for k, v in adsorbates.items()):
        raise ValueError("`adsorbates` must be a dictionary with atom symbols or Atoms.")
    
    if not isinstance(weights, dict) or not all(isinstance(k, str) and isinstance(v, (int, float)) for k, v in weights.items()):
        raise ValueError("`weights` must be a dictionary with atom symbols as keys and numeric values as weights.")
    
    if len(slab_pbc) != 3 or not all(isinstance(b, bool) for b in slab_pbc):
        raise ValueError("`slab_pbc` must be a list of three boolean values (e.g., [True, True, True]).")
    
    if len(bulk_pbc) != 3 or not all(isinstance(b, bool) for b in bulk_pbc):
        raise ValueError("`bulk_pbc` must be a list of three boolean values (e.g., [True, True, True]).")
    
    if not isinstance(cutoff, (int, float)) or cutoff <= 0:
        raise ValueError("`cutoff` must be a positive number.")
    
    if not isinstance(vacuum, (int, float)) or vacuum <= 0:
        raise ValueError("`vacuum` must be a positive number.")
    
    if not isinstance(axis, int) or axis not in [0, 1, 2]:
        raise ValueError("`axis` must be an integer and one of [0, 1, 2].")

    # 1. 获取slab的表面原子
    positions = slab.get_positions()
    z_max = max(positions[:, 2])
    z_min = min(positions[:, 2])
    
    # 定义一个简单的z轴阈值，判断顶部和底部的表面原子
    surface_atoms_top = [atom.index for atom in slab if atom.position[2] > z_max - 2.0]
    surface_atoms_bottom = [atom.index for atom in slab if atom.position[2] < z_min + 2.0]

    #print(slab.get_pbc())
    #print(bulk.get_pbc())
    # 2. 查找存在悬挂键的原子
    #print(surface_atoms_top + surface_atoms_bottom)
    adsorbate_info = []  # 用来存储要添加的吸附物信息
    for idx in surface_atoms_top + surface_atoms_bottom:
        atom_symbol = slab[idx].symbol
        if atom_symbol not in adsorbates:
            continue
        
        # slab 中原子的配位数
        indices_slab, offset_slab, distance_slab = get_neighbors(idx, slab.positions, array(slab.get_cell()), slab_pbc,cutoff)
        num_neighbors_slab = len(indices_slab)
        #print(idx, num_neighbors_slab, indices_slab, offset_slab, distance_slab)
        
        # bulk 中相同位置原子的配位数（使用相同 index）
        idx_bulk = find_matching_atom_in_bulk(slab.positions[idx], bulk.positions)
        if idx_bulk:
            indices_bulk, offset_bulk, distance_bulk = get_neighbors(idx_bulk, bulk.positions, bulk.get_cell(), bulk_pbc, cutoff)
            # 检查配位数是否一致
            # 对比目标配位数, also can compare with the bulk's
            # target_num_neighbors = target_coordination[atom_symbol]
            if num_neighbors_slab < len(indices_bulk):
                # 找出缺失的邻居
                unsaturated_neighbors = find_unsaturated_neighbors(offset_bulk, offset_slab, tolerance=1e-3)
                #print(f"原子 {idx} 缺失的邻居: {unsaturated_neighbors}")
                # 计算悬挂键方向
                # 根据 unsaturated_neighbors 的方向，添加钝化原子
                if unsaturated_neighbors:
                    for unsaturated_neighbor in unsaturated_neighbors:
                        if atom_symbol in adsorbates:
                            adsorbate_type = adsorbates[atom_symbol]
                            adsorbate_info.append((idx, atom_symbol, adsorbate_type, slab.positions[idx], unsaturated_neighbor))
        else:
            print(f"Error: Cannot find matching bulk atom for slab atom {idx}")
    #print(adsorbate_info)
    # 将氢和氢氧根加到 slab 中
    slab_with_passivation = slab.copy()

    counts = weights.copy()
    for symbol in counts:
        counts[symbol] = 0
    for adsorbate in adsorbate_info:
        idx = adsorbate[0]
        atom_symbol = adsorbate[1]
        adsorbate_type = adsorbate[2]
        atom_position = adsorbate[3]
        offset = adsorbate[4]
        adsorbate_position = atom_position + weights[atom_symbol] * offset
        add_adsorbate(slab_with_passivation, adsorbate_type, adsorbate_position[2]-z_max, adsorbate_position[:2])
        counts[atom_symbol] = counts[atom_symbol] + 1
    if if_print:
        for symbol in counts:
            print(f'Totally, dangling bonds of {symbol}: {counts[symbol]}')
    slab_with_passivation.wrap()
    slab_with_passivation.center(vacuum = vacuum, axis = axis)
    slab_with_passivation = slab_with_passivation[slab_with_passivation.numbers.argsort()]
    return slab_with_passivation




