from ase import Atoms
from math import sqrt, pow
from mymetal.build.extrfilm import cal_area
def find_cubic( atom: Atoms = None, type: str = '2D') -> Atoms:
    """
    Find the cubic cell of the atoms(2D/3D).
    """
    atoms = atom.copy()
    cell = atoms.get_cell()
    area = cal_area(atoms)
    volume = atoms.get_volume()
    if type == '2D' or '2d' or '2':
        cell[0,0] = sqrt(area)
        cell[1,1] = sqrt(area)
        cell[0,1] = 0
        cell[1,0] = 0
    elif type == '3D' or '3d' or '3':
        cell[0,0] = pow(volume, 1/3)
        cell[1,1] = pow(volume, 1/3)
        cell[2,2] = pow(volume, 1/3)
        cell[0,1] = 0
        cell[0,2] = 0
        cell[1,0] = 0
        cell[1,2] = 0
        cell[2,0] = 0
        cell[2,1] = 0
    else:
        print('Wrong type!')
    atoms.set_cell(cell)
    return atoms