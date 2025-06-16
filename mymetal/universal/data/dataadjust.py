"""
dataadjust submodule

This submodule provides functions for adjusting data structures. It includes functions for removing blank spaces, combining lists,
normalizing values, and justifying text. These functions are designed to streamline common tasks in materials science simulations
and data handling.

Functions:
    - rm_blank: Remove blank spaces from a list of strings.
    - my_add_list: Combine two lists and sort them.
    - normalize_float_int: Normalize values by converting them to floats and removing trailing zeros.
    - myjust: Justify text to the left or right.
    - my_down_up: Enlarge a list by moving its elements up or down.
"""

from inspect import getargvalues, currentframe, stack
from mymetal.universal.check.checkinput import check_input

##############################################################
############# adjust variable 
##############################################################

def rm_blank(extr_content, remove_string=' ', removed_string=''):
    """
    for 1-D list, remove blank\n
    like ["reached required accuracy","reached required accuracy"]  to \n
         [ "reachedrequiredaccuracy" , "reachedrequiredaccuracy" ]\n
    """

    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    for num in range(len(extr_content)):
        extr_content[num] = extr_content[num].replace(remove_string, removed_string)
        # print(extr_content)
        # print('\n')

def my_add_list(list1 = [], list2 = []):
    """"
    combine two list, and sort it
    """
    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    my_list = []
    for item in list1 + list2:
        if item not in my_list:
            my_list.append(item)
    my_list.sort()
    return my_list

def normalize_float_int(value):
    """
    Helper function to normalize values\n
    Convert to float and remove trailing zeros, if it's a float\n
    """

    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    try:
        float_value = float(value)
        return str(float_value).rstrip('0').rstrip('.')  # Remove trailing zeros and dot if it's an int
    except ValueError:
        return value  # Leave other types as they are

def myjust(left = True, num = 10, output = ''):
    """
    left just or right just of output\n
    num is digit of variable\n
    like 'dig' =>> 'dig       '  =>> '       dig'\n
    return char_out\n
    """
    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    char_out = ''
    if left:
        char_out = str(output).ljust(num)
    else:
        char_out = str(output).rjust(num)
    return char_out

def my_down_up(list, move_list_down_up = [0, 0]):
    """
    enlarge list\n
    like      [1,4], move_list = [-2,3]\n
         =>>  [-1,0,1,2,3,4,5,6,7]  (-1, 7)\n
    """

    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    my_list = []
    lower_move, upper_move = move_list_down_up
    for item in list:
        #my_item = [item + i for i in range(lower_move, upper_move + 1)]
        my_item = [item + i for i in range(lower_move, upper_move + 1)]
        my_list = my_add_list(my_list, my_item)
    return my_list
