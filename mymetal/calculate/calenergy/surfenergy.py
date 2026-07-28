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
    """Calculate surface energy from bulk and relaxed-slab energies.

    The implemented expression is

    ``(E_slab - E_bulk / N_bulk * N_slab) / (factor * area)``.
    Energy inputs are interpreted as eV and ``area`` as Å². ``energy_unit``
    selects the return unit; it does not change the input-energy unit.
    
    Args:
        bulk_energy: Total bulk-reference energy in eV.
        bulk_atoms_number: Number of atoms in the bulk reference.
        relaxed_surface_energy: Total relaxed-slab energy in eV.
        surface_atoms_number: Number of atoms in the slab.
        area: Area of one surface in Å².
        energy_unit: Return ``"eV"`` for eV/Å² or ``"J"`` for J/m².
        factor: Number of equivalent surfaces represented by the excess
            energy. A symmetric slab normally uses 2; verify this assumption
            for the model being analysed.
        
    Returns:
        Surface energy in eV/Å² when ``energy_unit="eV"`` or J/m² when
        ``energy_unit="J"``.

    Raises:
        ValueError: If an atom count, ``area``, or ``factor`` is not positive,
            or if ``energy_unit`` is neither ``"eV"`` nor ``"J"``.
    """
    
    # 检查 bulk_atoms_number 和 surface_atoms_number 是否是正整数
    if not isinstance(bulk_atoms_number, int) or bulk_atoms_number <= 0:
        raise ValueError(f"bulk_atoms_number must be a positive integer. Got: {bulk_atoms_number}")
    if not isinstance(surface_atoms_number, int) or surface_atoms_number <= 0:
        raise ValueError(f"surface_atoms_number must be a positive integer. Got: {surface_atoms_number}")
    
    # 面积和表面数进入分母；在这里拒绝 0，避免把输入错误拖成除零错误。
    if not isinstance(area, (int, float)) or area <= 0:
        raise ValueError(f"area must be a positive number. Got: {area}")
    if not isinstance(factor, int) or factor <= 0:
        raise ValueError(f"factor must be a positive integer. Got: {factor}")
    
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
