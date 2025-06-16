"""Some universal functions for electronic structure calculations."""

import numpy as np
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.plotter import BSPlotter

def get_n_band(
    bs: BandStructureSymmLine = None, ibands: int = 0, spin: Spin = Spin.up
) -> tuple[list, list]:
    """Get the k-point distances and corresponding energy values for a specified band.

    This function extracts the k-point distances and energy values for a given band index
    from the band structure data. It handles spin-polarized band structures and returns
    two lists: one for the k-point distances and one for the corresponding energy values.

    Args:
        bs (BandStructureSymmLine): The band structure object from Pymatgen.
        ibands (int): The index of the band to extract. Defaults to 0.
        spin (Spin): The spin direction (Spin.up or Spin.down). Defaults to Spin.up.

    Returns:
        tuple[list, list]: A tuple containing two lists:
            - k_distances (list): The k-point distances.
            - energies (list): The energy values for the specified band and spin.

    Examples:
        >>> k_distances, energies = get_n_band(bs, 2, Spin.up)
        >>> plt.plot(k_distances, energies)
    """
    from pymatgen.electronic_structure.plotter import BSPlotter

    if spin == Spin.up:
        spin = "1"
    else:
        spin = "2"  # Assuming "2" is used for Spin.down

    bsplot = BSPlotter(bs)
    data = bsplot.bs_plot_data()

    k_distances = list()
    energies = list()

    for xpath, epath in zip(data["distances"], data["energy"][spin]):
        k_distances.extend(xpath)
        energies.extend(epath[ibands])

    return k_distances, energies



def summarize_band_structure_info(bs: BandStructureSymmLine, spin: Spin = Spin.up, receive: list = None) -> tuple:
    """Summarize and print key information from the band structure.

    This function extracts and prints essential details from the band structure data,
    including distances, energy shapes, valence band maximum (VBM), conduction band minimum (CBM),
    tick labels, zero energy, metallicity, and band gap.

    Args:
        bs (BandStructureSymmLine): The band structure object from Pymatgen.
        spin (Spin): The spin direction (Spin.up or full). Defaults to Spin.up. if you use full, the full array will be returned.

    Returns:
        tuple: A tuple containing the extracted data:
            - distances (np.ndarray): The k-point distances.
            - energy (np.ndarray): The energy values.
            - vbm (float): The valence band maximum.
            - cbm (float): The conduction band minimum.
            - ticks (list): The tick labels.
            - zero_energy (float): The zero energy.
            - is_metal (bool): Whether the material is metallic.
            - band_gap (float): The band
    
    Examples:
        >>> summarize_band_structure_info(bs)
        Band Structure Data Summary:
        ---------------------------
        Keys in data: dict_keys(['distances', 'energy', 'vbm', 'cbm', 'ticks', 'zero_energy', 'is_metal', 'band_gap'])
        Shape of distances array: (1, 301)
        Shape of energy array for spin '1': (301,)
        Valence Band Maximum (VBM): 3.5516
        Conduction Band Minimum (CBM): 4.017
        Ticks: [0.00, 0.50, 1.00, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00]
               ['GAMMA', 'K', 'K', 'GAMMA', 'M', 'M', 'K', 'GAMMA', 'A']
        Zero Energy: 4.017
        Is Metal: False
        Band Gap: 0.4654
    """
    # Initialize the BSPlotter and get plot data
    bsplot = BSPlotter(bs)
    data = bsplot.bs_plot_data()

    # Gather information
    distance = data["distances"]
    if spin == Spin.up:
        energy = data["energy"]['1']
    else:
        energy = list(data["energy"].values())
    distances_shape = np.array(distance).shape
    energy_shape = np.array(energy).shape
    vbm = data['vbm']
    cbm = data['cbm']
    ticks = data['ticks']
    zero_energy = data['zero_energy']
    is_metal = data['is_metal']
    band_gap = data['band_gap']

    ticks['label'] = [label.replace('GAMMA', '$\\Gamma$') for label in ticks['label']]
    # Print the information
    print("Band Structure Data Summary:")
    print("---------------------------")
    print(f"Keys in data: {data.keys()}")
    print(f"Shape of distances array: {distances_shape}")
    print(f"Shape of energy array for spin '1': {energy_shape}")
    print(f"Valence Band Maximum (VBM): {vbm}")
    print(f"Conduction Band Minimum (CBM): {cbm}")
    print(f"Ticks: {[f'{dist:.2f}' for dist in ticks['distance']]}")
    print(f"       {ticks['label']}")
    print(f"Zero Energy: {zero_energy}")
    print(f"Is Metal: {is_metal}")
    print(f"Band Gap: {band_gap}")
    full = (distance, energy, vbm, cbm, ticks, zero_energy, is_metal, band_gap)
    receive = full
    return bsplot


def get_n_kpoints_band(bs: BandStructureSymmLine, ibands: int = 0, spin: Spin = Spin.up) -> tuple:
    """Get the k-points and corresponding energy values for a specified band.

    This function extracts the k-points and energy values for a given band index
    from the band structure data. It handles spin-polarized band structures and returns
    a list of k-points and corresponding energy values.

    Args:
        - bs (BandStructureSymmLine): The band structure object from Pymatgen.
        - ibands (int): The index of the band to extract. Defaults to 0.
        - spin (Spin): The spin direction (Spin.up or Spin.down). Defaults to Spin.up.
    
    Returns:
        - list: A list containing the k-points and corresponding energy values for the specified band and spin.
    
    Examples:
        >>> get_n_kpoints_band(bs, 2, Spin.up)
    """
    n = 0
    k = []
    for kpoints, e in zip(bs.kpoints, bs.bands[spin][ibands, :]):
        n += 1
        k.append([kpoints.frac_coords, e])
        print(
            "kx = %5.3f  ky = %5.3f  kz = %5.3f   eps(k) = %8.4f"
            % (tuple(kpoints.frac_coords) + (e,))
        )
    return k

