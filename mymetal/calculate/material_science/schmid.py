"""calculate schmid factors for only FCC slip systems (now)."""

import numpy as np
import pandas as pd

def cal_fcc_schmid_factors(normal_orientation: list=[1,1,6]) -> pd.DataFrame:
    """Calculates Schmid factors for FCC slip systems.

    This function computes the Schmid factors for FCC (Face-Centered Cubic) 
    slip systems given a specified loading direction. The calculation considers 
    four primary slip planes {111} and six linearly independent Burgers vectors 
    <110>.

    Args:
        normal_orientation (list, optional): The loading direction represented 
            as a three-element list of crystallographic indices. Defaults to [1,1,6].

    Returns:
        pd.DataFrame: A DataFrame containing Schmid factors and related parameters 
        for each valid slip system. The columns include:
            - `slip_plane`: The normal vector of the slip plane.
            - `dislocation`: The Burgers vector of the slip system.
            - `force_direction`: The applied loading direction.
            - `schmid_factor`: The computed Schmid factor.
            - `cos_lambda`: Cosine of the angle between force direction and slip direction.
            - `cos_phi`: Cosine of the angle between force direction and slip plane normal.

    Example:
        >>> df = cal_fcc_schmid_factors([1,1,6])
        >>> print(df.head())
    """
    # FCC晶体中只有4个滑移面
    slip_plane = np.array([
        [1, 1, 1],
        [-1, 1, 1],
        [1, -1, 1],
        [1, 1, -1]
    ])

    # FCC中共有12种完美位错的柏氏矢量, 但实际上线性无关的只有6个
    perfect_dislocation = 0.5 * np.array([
        [1, 0, 1],
        [-1, 1, 0],
        [0, 1, 1],
        [-1, 0, 1],
        [0, -1, 1],
        [1, 1, 0]
    ])


    schmid_factors_dict = {}
    schmid_factors = []
    planes = []
    dislocations = []
    cos_lambda_list = []
    cos_phi_list = []
    force_direction = []
    # 遍历所有滑移面上的所有可滑位错的柏氏矢量
    for i, vector in enumerate(slip_plane):
        for j, burgers in enumerate(perfect_dislocation):
            # 确保位错矢量在对应的滑移面上
            if abs(np.dot(vector, burgers)) < 1e-3:
                cos_lambda = np.dot(normal_orientation, burgers) / (np.linalg.norm(normal_orientation) * np.linalg.norm(burgers))
                cos_phi = np.dot(normal_orientation, vector) / (np.linalg.norm(normal_orientation) * np.linalg.norm(vector))
                schmid_factor = abs(cos_lambda * cos_phi)
                schmid_factors.append(schmid_factor)
                planes.append(vector)
                dislocations.append(burgers)
                cos_lambda_list.append(cos_lambda)
                cos_phi_list.append(cos_phi)
                force_direction.append(normal_orientation)
    schmid_factors_dict['slip_plane'] = planes
    schmid_factors_dict['dislocation'] = dislocations
    schmid_factors_dict['force_direction'] = force_direction
    schmid_factors_dict['schmid_factor'] = schmid_factors
    schmid_factors_dict['cos_lambda'] = cos_lambda_list
    schmid_factors_dict['cos_phi'] = cos_phi_list

    df = pd.DataFrame(schmid_factors_dict)
    return df