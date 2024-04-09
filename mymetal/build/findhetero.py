from mymetal.build.extrfilm import cal_area, my_extr_thick
from ase import Atoms
import math
from ase.build import make_supercell
from ase.build.tools import sort
from numpy import array, zeros, concatenate
from ase.visualize import view

# set_length Done
# find_hetero Not Done!!!!!!!!!!!!!!!!!!!!
# construct_hetero is Done!!!!!!!!!!!!

# test1
# find_hetero(gold_fcc_film, sio2_film_prim)
# set_length(sio2_film_prim, 50)
# print(sio2_film_prim)

# test2
# from hetbuilder.algorithm import CoincidenceAlgorithm
# from hetbuilder.plotting import InteractivePlot
# from ase.io.vasp import read_vasp, write_vasp
# bottom = read_vasp('CONTCAR_SiO_film')
# top = read_vasp('CONTCAR_FCC_film')
# # we set up the algorithm class
# alg = CoincidenceAlgorithm(bottom, top)
# # we run the algorithm for a choice of parameters
# results = alg.run(Nmax = 5, Nmin = 0, 
#                   tolerance = 2, weight = 0)
# j = results[1]
# print(type(j))
# print(j.bottom)
# print(j.top)
# print(j.stack)
# print(j.strain)
# view(bottom)
# view(j.stack)
# view(build_supercells(bottom, top, j.M, j.N, j.angle, j._weight))
# if j.stack is same as the supercells 


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

# Not Done
################################
# DONE 

def build_supercells(primitive_bottom: Atoms = None, 
                         primitive_top: Atoms = None,
                         M: array = None,
                         N: array = None,
                         angle_z: float = None, 
                         weight: float = 0.5,
                         distance: float = 3.5,
                         vacuum: float = 15,
                         pbc: list = [True, True, False], 
                         reorder=True ) -> Atoms:
    """
    For construct supercell for two primitive cell known transition matrix and rotation angle of top layer around z-dir\n
    Input: must be primitive cell\n
           M, N: bottom, top matrix\n
           angle_z: units: degree, +z direction\n
           distance, vacuum: units: Angstron\n
           weight: must be 0~1, fixing bottom~top\n
           pbc, reorder: the default value is good\n
    Output: 
           heterostructure
    """
    bottom = primitive_bottom.copy()
    top = primitive_top.copy()
    #print(top)
    # make_supercell
    bottom_sup = make_supercell(bottom, M)
    top_sup = make_supercell(top, N)
    #print(bottom_sup)

    # rotate top layer around z, units: degree
    top_sup.rotate(angle_z, 'z', rotate_cell=True)
    #print(top_sup)
    # translate from hetbuilder (atom_functions.cpp/.h) stack_atoms function
    stack = stack_atoms(bottom_sup, top_sup, weight, distance, vacuum, pbc, reorder)
    #print(1)
    return stack

def stack_atoms(bottom_sup: Atoms = None, top_sup: Atoms = None, weight: float = 0.5, distance: float = 3.5, vacuum: float = 15, 
                pbc: list = [True, True, False], reorder: bool=True, shift_tolerance: float = 1e-5) -> Atoms:
    """
    After searching the coincidenced supecell, we get the transition matrix(M for bottom, N for top), and make two supercells, and then try to match them.\n
    Additionally, bottom, top, weight, distance, vacuum.\n
    usage: updating...
    Noted: the cell C = A + weight * [B - A], A - bottom, B - top
    """
    # get the max, min z-position 
    bottom = bottom_sup.copy()
    top = top_sup.copy()
    min_z1, max_z1 = bottom.positions[:,2].min(), bottom.positions[:,2].max()
    min_z2, max_z2 = top.positions[:,2].min(), top.positions[:,2].max()
    bottom_thickness = max_z1 - min_z1
    top_thickness = max_z2 - min_z2
    #print(bottom_thickness)
    #print(distance)
    shift = bottom_thickness + distance

    # add distance
    bottom.positions[:,2] -= min_z1
    bottom.positions[:,2] += shift_tolerance
    top.positions[:,2] -= min_z2
    top.positions[:,2] += shift
    top.positions[:,2] += shift_tolerance

    # generate lattice
    new_lattice = zeros((3,3))
    
    # C = A + weight * [B - A]
    for i in range(2):  
        for j in range(2):
            #print(i,j)
            new_lattice[i, j] = bottom.cell[i, j] + weight * (top.cell[i, j] - bottom.cell[i, j])
    new_lattice[2,2] = bottom_thickness + distance + top_thickness + vacuum

    # combine the position and symbols information
    all_positions = concatenate([bottom.positions, top.positions], axis=0)
    all_symbols = bottom.get_chemical_symbols() + top.get_chemical_symbols()

    # create new Atoms
    stack = Atoms(symbols=all_symbols, positions=all_positions, cell=new_lattice, pbc=pbc)
    stack.center()
    if reorder:
        stack = sort(stack)
    return stack


def scale_cell_xy(atoms_origin: Atoms = None, new_cell: array = None):
    """
    After calculating C - supercell, try to scall xy of cell, but fixed z.\n
    It's important for fixing someone's lattice constant, otherwise, the z-dir constant'll changed!
    """
    atoms = atoms_origin.copy()
    # 获取原始的笛卡尔坐标
    original_cart_pos = atoms.get_positions()
    
    # 获取按新晶胞缩放的分数坐标
    scaled_positions = atoms.get_scaled_positions()

    # 更新晶胞，注意这里假设new_cell是一个numpy数组或类似的可直接赋值给atoms.cell的结构
    atoms.set_cell(new_cell)

    # 将缩放的分数坐标转换回笛卡尔坐标，这里直接使用Atoms对象的.set_scaled_positions方法
    atoms.set_scaled_positions(scaled_positions)

    # 保持z方向上的笛卡尔坐标不变
    updated_cart_pos = atoms.get_positions()
    updated_cart_pos[:, 2] = original_cart_pos[:, 2]

    # 更新原子的位置
    atoms.set_positions(updated_cart_pos)


