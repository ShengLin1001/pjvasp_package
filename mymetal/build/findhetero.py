from mymetal.build.extrfilm import cal_area, my_extr_thick
from ase import Atoms
import math

# set_length Done
# find_hetero Not Done!!!!!!!!!!!!!!!!!!!!

def find_hetero( up_layer: Atoms = None,
                bot_layer: Atoms = None,
                inter_dis: float = None,
                fix_bot_lat: bool = True) -> Atoms:
    """
    Find the hetero interface between two layers.\n
    We think the bot_layer is a bulk-like material.\n
    so the lattice constant of it is fixed now.\n
    """
    up_layer.center(vacuum=0, axis=2)
    bot_layer.center(vacuum=0, axis=2)
    area_up = cal_area(up_layer)
    area_bot = cal_area(bot_layer)
    thick_up = my_extr_thick(up_layer)
    thick_bot = my_extr_thick(bot_layer)
    #view(up_layer)
    print(up_layer)
    print(bot_layer)
    print(area_bot/area_up)

def set_length( atoms: Atoms = None,
                length: float = None,
                axis: int = 2) -> Atoms:
    """
    Set the length of the atoms along the axis.\n
    """
    cell = atoms.get_cell()
    cell[axis] = cell[axis]/magnitude(cell[axis])*length
    atoms.set_cell(cell)
    return atoms

def magnitude(vector: list = None) -> float:
    """
    function definition to compute magnitude o f the vector
    """
    return math.sqrt(sum(pow(element, 2) for element in vector))

# find_hetero(gold_fcc_film, sio2_film_prim)
# set_length(sio2_film_prim, 50)
# print(sio2_film_prim)