import numpy as np

def my_sort(refx: list=[], Etot: list=None, Eent: list = None, pres: list = None, reverse: bool = False) -> tuple:
    """
    Sort convergence data based on input parameter (cutoff or k-point grid).

    Args:
        refx (list): List or array of input values used for sorting (1D or 2D).
        Etot (list): Total energy values.
        Eent (list): Electronic entropy values.
        pres (list): Pressure values.
        reverse (bool): If True, sort in descending order.

    Returns:
        tuple:
            x_sorted (np.ndarray): Sorted input values.
            Etot_sorted (np.ndarray): Sorted total energies.
            Eent_sorted (np.ndarray): Sorted entropy values.
            pres_sorted (np.ndarray): Sorted pressure values.
    """
    refx = np.array(refx)
    Etot = np.array(Etot)
    Eent = np.array(Eent)
    pres = np.array(pres)
    x = refx
    if x.ndim == 1:
        sort_idx = np.argsort(x)        # 1D 
    else:
        sort_idx = np.argsort(x[:, 0])  # 2D 

    if reverse:
        sort_idx = sort_idx[::-1]

    x_sorted = x[sort_idx]
    Etot_sorted = np.array(Etot)[sort_idx] 
    Eent_sorted = np.array(Eent)[sort_idx] 
    pres_sorted = np.array(pres)[sort_idx]
    return x_sorted, Etot_sorted, Eent_sorted, pres_sorted
