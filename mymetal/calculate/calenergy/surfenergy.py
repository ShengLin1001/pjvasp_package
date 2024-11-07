"""
surfenergy module

This module provides functions for calculating surface energy of materials. 
Surface energy is an important property for materials science, especially 
in the study of material interfaces, thin films, and catalysis.

Functions:
    - cal_surface_energy: Calculate the surface energy of a given material structure.
"""

from ase import Atoms

def cal_surface_energy(bulk_energy: float = None,
                       bulk_atoms_number: int = None, 
                       relaxed_surface_energy: float = None,
                       surface_atoms_number: int = None,
                       area: float = None,
                       energy_unit: str = 'eV',
                       factor: int = 2) -> float:
    """
    Calculate surface energy based on bulk and relaxed surface energy data.
    
    Args:
        bulk_energy (float): Total bulk energy in the specified unit.
        bulk_atoms_number (int): Number of atoms in the bulk structure (must be a positive integer).
        relaxed_surface_energy (float): Total relaxed surface energy in the specified unit.
        surface_atoms_number (int): Number of atoms in the surface structure (must be a positive integer).
        area (float): Surface area in Å² (must be a non-negative float).
        energy_unit (str): Unit of energy, either 'eV' (electron volts) or 'J' (joules). Default is 'eV'.
        
    Returns:
        float: Surface energy in eV/Å² or J/m², depending on the selected energy unit.

    Raises:
        ValueError: If the specified energy unit is not 'eV' or 'J'.
        ValueError: If bulk_atoms_number or surface_atoms_number is not a positive integer.
        ValueError: If area is not a non-negative float.
    """
    
    # 检查 bulk_atoms_number 和 surface_atoms_number 是否是正整数
    if not isinstance(bulk_atoms_number, int) or bulk_atoms_number <= 0:
        raise ValueError(f"bulk_atoms_number must be a positive integer. Got: {bulk_atoms_number}")
    if not isinstance(surface_atoms_number, int) or surface_atoms_number <= 0:
        raise ValueError(f"surface_atoms_number must be a positive integer. Got: {surface_atoms_number}")
    
    # 检查 area 是否为非负数
    if not isinstance(area, (int, float)) or area < 0:
        raise ValueError(f"area must be a non-negative number. Got: {area}")
    
    # 能量转换系数: 1 eV/Å² = 16.021766 J/m²
    conversion_factor = 1.0
    if energy_unit == 'J':
        conversion_factor = 16.021766  # 将 eV/Å² 转换为 J/m²
    elif energy_unit != 'eV':
        raise ValueError(f"Invalid energy unit: {energy_unit}. Supported units are 'eV' and 'J'.")

    # 计算每个原子的体相能量
    bulk_energy_per_atom = bulk_energy / bulk_atoms_number  # eV/atom 或 J/atom
    
    # 计算表面能
    surface_energy = (relaxed_surface_energy - bulk_energy_per_atom * surface_atoms_number) / (factor * area)
    
    # 根据单位进行转换
    surface_energy_converted = surface_energy * conversion_factor
    
    return surface_energy_converted