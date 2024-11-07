"""
printafter submodule

This submodule contains functions for printing messages after a certain condition is met.

Functions:
    - print_after_read: Prints a message if all the specified items have been read.
    - print_after_blank: Prints a message if the specified item is blank.
    - print_after_not_supported: Prints a message if the specified item is None.
    - print_after_cant_read: Prints a message if the specified item is None.
"""

from inspect import getargvalues, currentframe, stack
from mymetal.universial.check.checkinput import check_input

##############################################################
############# note after read
##############################################################

def print_after_read(count_read, stretch_factor = [], special_name = '', name = ''):
    """
    note if read all {special_name} - OUTCAR\n
    {name} is class name {self.name}\n
    nead a count_read variable\n
    """

    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    if count_read:
        if count_read == len(stretch_factor):
            print(f'Read all {special_name} - {name}')
        else:
            print(f'Not all {special_name} are readed - {name}')
    else:
        print_after_blank('count_read', calling_function, [])

def print_after_blank(char_name = '', calling_function = '', specified_blank = ''):
    """
    combined two lines\n
            print(f"the value is None! - {calling_function}")\n
            char_out = ''\n
    """
    print(f"the || {char_name} || is blanked! - {calling_function}")
    return specified_blank

def print_after_not_supported(char_name = '', calling_function = '', specified_blank = ''):
    """
    combined two lines\n
            print(f"the value is None! - {calling_function}")\n
            char_out = ''\n
    """
    print(f"the || {char_name} || is None! - {calling_function}")
    return specified_blank

def print_after_cant_read(char_name = '', calling_function = '', specified_blank = ''):
    """
    combined two lines\n
            print(f"the value is None! - {calling_function}")\n
            char_out = ''\n
    """

    print(f"the || {char_name} || is None! - {calling_function}")
    return specified_blank