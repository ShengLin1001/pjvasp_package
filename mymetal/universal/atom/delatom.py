"""
delatom submodule

This submodule provides functions for deleting atoms from a structure. It includes functions for deleting atoms based on
their positions, types, or both. These functions are designed to streamline common tasks in materials science simulations
and data handling.

Functions:
    - mydel_pos_type: Delete atoms based on their positions and types.
    - check_position: Check if an atom is within a specified position range.
"""

from ase import Atoms
from numpy import ndarray, array
from mymetal.universal.print.printafter import print_after_blank
from inspect import getargvalues, currentframe, stack
from mymetal.universal.check.checkinput import check_input

def mydel_pos_type(atoms: Atoms = None,
          position_strict: list = None,
          type_strict: list = None) -> Atoms:
    """
    using to delete some atoms in specified region or/and type\n
    ps: [xl, xh, yl, yh, zl, zh]\n
    ts: ['H', 'O', ... ]\n
    temp = mydel_pos_type(temp, [float('-inf'), float('inf'), float('-inf'), float('inf'), float('-inf'), 0.9])\n
    temp = mydel_pos_type(temp, [float('-inf'), float('inf'), float('-inf'), float('inf'), 10.18, float('inf')])\n
    """
    calling_function = stack()[1].function
    #check_input(getargvalues(currentframe()).locals)
    temp = atoms.copy()
    p = temp.positions
    ps = position_strict
    ts = type_strict
    if atoms:
        if ps is not None and ts is None:
            del temp[[atom.index for atom in temp if check_position(ps, atom.position) ]]
        elif ps is None and ts is not None:
            del temp[[atom.index for atom in temp if atom.symbol in ts ]]
        elif ps is not None and ts is not None:
            del temp[[atom.index for atom in temp if ((atom.symbol in ts) and check_position(ps, atom.position) )]]
        else:
            print_after_blank('the position strict and type strict', calling_function)
    else:
        print_after_blank('the input atoms', calling_function)
    temp = temp[temp.numbers.argsort()]
    return temp

def check_position(position_strict: ndarray = None,
                   position: ndarray = None) -> bool:
    """checking the atom is in ps [xo, xh, yo, yh, zo, zh]"""
    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    tag = False
    p = position
    ps = position_strict
    if position is None:
        print_after_blank('the atom position', calling_function)
    else:
        if p[0] > ps[0] and p[0] < ps[1] and p[1] > ps[2] and p[1] < ps[3] and p[2] > ps[4] and p[2] < ps[5]:
            tag = True
    return tag

