from ase import Atoms
import numpy as np

def calculate_strain_matrix(initial_atoms: Atoms, deformed_atoms: Atoms):
    """_summary_

    :param initial_atoms: _description_
    :type initial_atoms: Atoms
    :param deformed_atoms: _description_
    :type deformed_atoms: Atoms
    :return: _description_
    :rtype: _type_
    """
    # 获取初始和变形后的晶格矩阵
    initial_cell = initial_atoms.cell
    deformed_cell = deformed_atoms.cell

    # 计算变形矩阵 F = B * A^(-1)
    deformation_matrix = np.dot(deformed_cell, np.linalg.inv(initial_cell))

    # 计算应变矩阵 ε = 1/2 * (F^T * F - I)
    strain_matrix = 0.5 * (np.dot(deformation_matrix.T, deformation_matrix) - np.eye(3))

    return strain_matrix

# 创建示例Atoms对象
# 假设这里有两个Atoms对象，一个是未变形的，另一个是变形后的
initial_atoms = Atoms('Si2', positions=[[0, 0, 0], [1.35, 1.35, 1.35]], cell=[2.7, 2.7, 2.7], pbc=[1, 1, 1])
deformed_atoms = Atoms('Si2', positions=[[0, 0, 0], [1.4, 1.4, 1.4]], cell=[2.8, 2.8, 2.8], pbc=[1, 1, 1])

# 计算应变矩阵
strain_matrix = calculate_strain_matrix(initial_atoms, deformed_atoms)
print("Strain Matrix:\n", strain_matrix)
