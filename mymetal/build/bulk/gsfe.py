"""
gsfe module

This module provides functions to create generalized stacking fault energy (GSFE) models for different crystal structures,
including FCC and HCP. The main function `create_gsfe_model` allows users to specify the type of GSFE model, lattice constants,
and supercell size to generate an ASE Atoms object suitable for GSFE calculations.

Functions:
    - create_gsfe_model(gsfe_type: str, a: float, c: float, size: np.array) -> Atoms: Creates a GSFE model based on the specified parameters.

"""

from mymetal.build.bulk.create import create_fcc_111, create_hcp_basal, create_hcp_prism1, create_hcp_prism2
from mymetal.build.bulk.create import vasp_create_hcp_pyr1, vasp_create_hcp_pyr2, vasp_create_fcc_100
from ase.io import write, read
from myvasp import vasp_func as vf
import numpy as np
import os
from ase import Atoms
from ase.constraints import FixCartesian, FixScaled


def create_gsfe_model(gsfe_type: str = None, a: float = None, c: float = None, size: np.array = None) -> Atoms:
    """创建指定滑移体系的 GSFE 原子模型。

    根据 `gsfe_type` 选择 FCC 或 HCP 晶体取向，并调用相应的结构构建函数生成
    用于广义层错能（GSFE）计算的 ASE Atoms 对象。当未显式指定 `size` 时，
    使用内置的默认超胞尺寸。

    Args:
        gsfe_type: GSFE 模型类型。可选值包括：
            "FCC_111"、"FCC_100"、"HCP_basal"、"HCP_prism1w"、
            "HCP_pyr1w" 和 "HCP_pyr2"。
        a: 晶格常数 a。
        c: HCP 晶格常数 c。用于计算轴比 c/a。
        size: 超胞尺寸。若为 None，则根据 `gsfe_type` 使用默认尺寸。

    Returns:
        生成的 GSFE 原子模型。
    """

    # size 默认值字典，针对不同的 gsfe_type 提供默认超胞尺寸
    dict_gsfe_type = {
        "FCC_111": [1, 1, 7],
        "FCC_100": [1, 1, 5],
        "HCP_basal": [1, 1, 10],
        "HCP_prism1w": [1, 1, 5],
        "HCP_pyr1w": [1, 1, 10],
        "HCP_pyr2": [1, 1, 5]
    }

    if gsfe_type in dict_gsfe_type.keys() and size is None:
        size = dict_gsfe_type[gsfe_type]

    if a is not None and c is not None:
        ca = c / a
    size = np.array(size)

    if gsfe_type == "HCP_basal":  # hcp
        # 2n
        atoms = create_hcp_basal(a=a, c=c, size = size)
    elif gsfe_type == "HCP_prism1w":
        # 4n
        atoms = create_hcp_prism1(a=a, c=c, size = size)
    elif gsfe_type == "HCP_pyr1w":
        # 2n
        atoms = vasp_create_hcp_pyr1(a, ca, size, bp=33)
    elif gsfe_type == "HCP_pyr2":
        # 4n
        atoms = vasp_create_hcp_pyr2(a, ca, size)
    elif gsfe_type == "FCC_111":  # fcc
        # 3n
        atoms = create_fcc_111(a=a, size=size)
    elif gsfe_type == "FCC_100":
        # 4n
        atoms = vasp_create_fcc_100(a, size)
    return atoms
