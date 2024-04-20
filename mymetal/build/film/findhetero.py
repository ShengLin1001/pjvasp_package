from mymetal.build.film.extrfilm import cal_area, my_extr_thick
from ase import Atoms
import math
from ase.build import make_supercell, stack
from ase.build.tools import sort
from numpy import array, zeros, concatenate
from ase.visualize import view
from hetbuilder.algorithm import CoincidenceAlgorithm
from hetbuilder.plotting import InteractivePlot
import matplotlib
from mymetal.universial.atom.moveatom import move_atoms

# set_length Done
# find_hetero Done!!!!!!!!!!!!!!!!!!!!
# construct_hetero is Done!!!!!!!!!!!!
# test see below

# take from hetbuilder/algorithm.py Line 333 def run()
def find_hetero(bottom: Atoms = None,
                top: Atoms = None,
                Nmax: int = 10,
                Nmin: int = 0,
                angles: list = [],
                angle_limits: tuple = (0, 90),
                angle_stepsize: float = 1.0,
                tolerance: float = 0.1,
                weight: float = 0.5,
                distance: float = 3.5,
                vacuum: float = 15,
                standardize: bool = False,
                no_idealize: bool = False,
                symprec: float = 1e-5,
                angle_tolerance: float = 5,
                verbosity: int = 0,
                plot: bool = True,
                plot_type: str = 'TkAgg' , 
                move_list: list = [0.1, 0.1, 0.1],
                center: bool = True,
                pbc: list = [True, True, False], 
                reorder=True
                ) -> list:
    """
    Find a heterostucture using hetbuilder API.\n
    ----------
    Executes the coincidence lattice algorithm.
    ----------

    Args:
        Nmax (int): Maximum number of translations. Defaults to 10.
        Nmin (int): Minimum number of translations. Defaults to 0.
        angles (list): List of angles in degree to search. Takes precedence over angle_limits and angle_stepsize.
        angle_limits (tuple): Lower and upper bound of angles to look through with given step size by angle_stepsize. Defaults to (0, 90) degree.
        angle_stepsize (float): Increment of angles to look through. Defaults to 1.0 degree.
        tolerance (float): Tolerance criterion to accept lattice match. Corresponds to a distance in Angström. Defaults to 0.1.
        weight (float): The coincidence unit cell is C = A + weight * (B-A). Defaults to 0.5.
        distance (float): Interlayer distance of the stacks. Defaults to 4.0 Angström.
        vacuum (float): Thickness of the vacuum layer of the stacks. Defaults to 15.0 Angström.
        standardize (bool): Perform spglib standardization. Defaults to true.
        no_idealize (bool): Does not idealize unit cell parameters in the spglib standardization routine. Defaults to False.
        symprec (float): Symmetry precision for spglib. Defaults to 1e-5 Angström.
        angle_tolerance (float): Angle tolerance fo the spglib `spgat` routines. Defaults to 5.
        verbosity (int): Debug level for printout of Coincidence Algorithm. Defaults to 0.
        plot_type: Default is TkAgg or Qt5Agg, to interactive plot using jupyter notebook.
        move_list: Default is [0.1, 0.1, 0.1], avoid the value error of 0.0...
        center: if center, default is True.
        pbc: [True, True, False]
        reorder: True

    Returns:
        list : A list of :class:`~hetbuilder.algorithm.Interface`.

    """
    # we set up the algorithm class
    alg = CoincidenceAlgorithm(bottom, top)
    # we run the algorithm for a choice of parameters
    results = alg.run(Nmax, Nmin, angles, angle_limits, angle_stepsize,
                tolerance, weight, distance, vacuum, standardize,
                no_idealize, symprec, angle_tolerance, verbosity)
    if results is not None:
        for result in results:
            #result.stack = build_supercells(result.bottom, result.top, result.M, result.N, result.angle,  result._weight, distance, vacuum, pbc, reorder)
            result.stack = move_atoms(result.stack, move_list)
            if center:
                result.stack.center()
            result.stack = result.stack[result.stack.numbers.argsort()]
        if plot:
            matplotlib.use(plot_type) # if failed,  pip install Qt5Agg{PyQt5  PySide2} TkAgg{tkinter}
            iplot = InteractivePlot(bottom=bottom, top=top, results=results, weight=weight)
            iplot.plot_results()
    else:
        print("Sorry, we didn't find any heterostructure")
    return results

# take from hetbuilder/ploting.py Line 211 class InteractivePlot
def my_plot_results(bottom: Atoms = None,
        top: Atoms = None,
        results: list = None,
        weight: float = 0.5,
        plot_type: str = 'Qt5Agg') -> None:
    """ 
    Interactive visualization of the results via matplotlib. 
    ----------
    
    Args:
        bottom (ase.atoms.Atoms): Lower layer as primitive.
        top (ase.atoms.Atoms): Upper layer as primitive.
        results (list): List of :class:`~hetbuilder.algorithm.Interface` returned from the coincidence lattice search.
        weight (float, optional): Weight of the supercell.
        plot_type: Default is Qt5Agg, to interactive plot using jupyter notebook.
    """
    matplotlib.use(plot_type) # if failed,  pip install PyQt5  PySide2
    iplot = InteractivePlot(bottom=bottom, top=top, results=results, weight=weight)
    iplot.plot_results()

def set_length( atoms: Atoms = None,
                length: float = None,
                axis: int = 2) -> Atoms:
    """
    Set the length of the atoms along the axis.\n
    """
    atoms_copy = atoms.copy()
    cell = atoms_copy.get_cell()
    cell[axis] = cell[axis]/magnitude(cell[axis])*length
    atoms_copy.set_cell(cell)
    return atoms_copy

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
                         M: array = [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                         N: array = [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                         angle_z: float = 0, 
                         weight: float = 0.5,
                         distance: float = 3.5,
                         vacuum: float = 15,
                         pbc: list = [True, True, False], 
                         reorder=True,
                          if_stack = True ) -> Atoms:
    """
    For construct supercell for two primitive cell known transition matrix and rotation angle of top layer around z-dir
    ----------
    
    Input: must be primitive cell\n
           M, N: bottom, top matrix\n
           angle_z: units: degree, +z direction\n
           distance, vacuum: units: Angstron\n
           weight: must be 0~1, fixing bottom~top\n
           pbc, reorder: the default value is good\n
           **if_stack**: if false, return bottom, top after make supercell and rotate it.\n
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
    if if_stack:
        stack = stack_atoms(bottom_sup, top_sup, weight, distance, vacuum, pbc, reorder)
        #print(1)
        return stack
    else:
        return bottom_sup, top_sup

def stack_atoms(bottom_sup: Atoms = None, top_sup: Atoms = None, weight: float = 0.5, distance: float = 3.5, vacuum: float = 15, 
                pbc: list = [True, True, False], reorder: bool=True, shift_tolerance: float = 1e-5) -> Atoms:
    """
    After searching the coincidenced supecell, we get the transition matrix(M for bottom, N for top), and make two supercells, and then try to match them.\n
    Additionally, bottom, top, weight, distance, vacuum.\n
    usage: updating...
    Noted: the cell C = A + weight * [B - A], A - bottom, B - top.\n
    has error, not cartesian, should be scaled positions.\n
    has been fixed!\n
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

    # scale x-y position of new cell and fix z-direction
    bottom = scale_cell_xy(atoms_origin = bottom, new_cell = new_lattice)
    top = scale_cell_xy(atoms_origin = top, new_cell = new_lattice)

    # combine the position and symbols information
    all_positions = concatenate([bottom.positions, top.positions], axis=0)
    all_symbols = bottom.get_chemical_symbols() + top.get_chemical_symbols()

    # create new Atoms
    stack = Atoms(symbols=all_symbols, positions=all_positions, cell=new_lattice, pbc=pbc)
    stack.center()
    if reorder:
        stack = sort(stack)
    return stack


# take from https://github.com/hongyi-zhao/hetbuilder
def scale_cell_xy(atoms_origin: Atoms = None, new_cell: array = None) -> Atoms:
    """
    After calculating C - supercell, try to scall xy of cell, but fixed z.\n
    It's important for fixing someone's lattice constant, otherwise, the z-dir constant'll changed!
    """
    atoms = atoms_origin.copy()
    original_cart_pos = atoms.get_positions()
    scaled_positions = atoms.get_scaled_positions()

    atoms.set_cell(new_cell)
    atoms.set_scaled_positions(scaled_positions)

    updated_cart_pos = atoms.get_positions()
    updated_cart_pos[:, 2] = original_cart_pos[:, 2]

    atoms.set_positions(updated_cart_pos)
    return atoms

def split_model(hetero: Atoms = None, bottom_type: list = None, top_type: list = None, adjust: bool=False, vacuum: float = 15.0, axis: int = 2) -> list:
    """
    Split a heterostructure to substrate and film according to the atom type.
    ----------
    Usage: 
    ----------
        >> `bottom, top = split_model(hetero, ['Si','O','H'], ['Au'])`
    """
    bottom = Atoms()
    bottom.set_cell(hetero.get_cell())
    top = Atoms()
    top.set_cell(hetero.get_cell())
    if Atoms:
        if bottom_type and top_type:
            for atom in hetero:
                if atom.symbol in bottom_type:
                    bottom.append(atom)
                if atom.symbol in top_type:
                    top.append(atom)
        else:
            print("Warning: Atom type not specified in bottom_type or top_type.")
    else:
        print("Warning: hetero is None type")
    if adjust:
        top.center(vacuum, axis)
        bottom.center(vacuum, axis)
    return [bottom, top]


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

