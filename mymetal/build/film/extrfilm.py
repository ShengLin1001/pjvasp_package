"""
extrfilm submodule

This submodule provides functions for extracting information from thin film
structures, including the lattice parameters, thickness, total energy, and
surface area.

Functions:
    - my_extr_lattice: Extracts the in-plane lattice parameter from a thin film
        structure.
    - my_extr_thick: Extracts the thickness of a thin film structure.
    - my_extr_etot: Extracts the total energy from a VASP OUTCAR file.
    - cal_area: Calculates the surface area of a thin film structure.
"""

from ase.utils import reader
from ase import Atoms
from mymetal.universial.data.datachange import list_to_char
from mymetal.post.newmain import myfindall
from numpy import array, cross

##############################################################
############# calculate the lattice
##############################################################

@reader
def my_extr_lattice(atoms: Atoms = None, stretch_type: str = None, tolerance = 1e-6) -> float:
    """
    For primitive cell of film, FCC / HCP\n
    extract the inplane cell length
    """
    my_cell = atoms.cell.cellpar()
    lattice = 0.0
    if stretch_type == 'a1':
        lattice = my_cell[0]
    elif stretch_type == 'a2':
        lattice = my_cell[1]
    elif stretch_type == 'a1a2':
        if abs(my_cell[0] - my_cell[1]) < tolerance:
            lattice = my_cell[1]
        else:
            raise TypeError(f'{my_cell[0]} and {my_cell[1]} is an unsame')
    return lattice

@reader
def my_extr_thick(atoms: Atoms = None, direction: int = 2) -> float:
    """For film, extract the thickness perpendiculat to the plane, almost along z-direction"""
    positions = atoms.get_positions()
    position = positions[:,direction]
    position.sort()
    thick = float(max(position) - min(position))
    return thick

@reader
def my_extr_etot(filename: str = None,
                 etot_config = [{'char': "energy(sigma->0)", 'basic_pattern': r'energy\(sigma->0\) =\s*', 'match_type': 'float',
                                'switch': 1, 'basic_content1':'','myindex1':0  , 'search_once' : False}]
                ) -> float:
    """extract etot from OUTCAR"""

    fd = filename
    etot = []
    if fd is not None:
        for config in etot_config:
            myfindall(fd, **config, mycontent=etot)
    etot = float(list_to_char(etot))
    return etot

# calculate the surface area in xy plane
def cal_area(atoms: Atoms = None) -> float:
    lattice = array(atoms.get_cell())
    a = lattice[0,:]
    b = lattice[1,:]
    area = abs(cross(a, b)[2])
    return area



