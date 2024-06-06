from ase import Atoms
from math import sqrt
from mymetal.build.film.extrfilm import cal_area
def find_cubic( prim: Atoms = None, type: str = 'hcp') -> Atoms:
    """
    Find the xy plane cubic cell of the primitive atoms(hcp/fcc film).\n
    input [a,b,c,any, any, 90]
    """
    atoms = prim.copy()
    cell = atoms.get_cell()
    area = cal_area(atoms)
    #volume = atoms.get_volume()
    if type == 'hcp' or 'fcc':
        cell[0,0] = sqrt(area/sqrt(3.0))
        cell[1,1] = sqrt(3.0)*cell[0,0]
        cell[0,1] = 0
        cell[1,0] = 0
    else:
        print('Wrong type!')
    atoms.set_cell(cell)
    return atoms