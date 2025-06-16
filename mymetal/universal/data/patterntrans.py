"""
patterntrans submodule

This submodule provides functions for transforming data patterns. It includes functions for generating regular expressions
from input data and parameters. These functions are designed to streamline common tasks in materials science simulations
and data handling.

Functions:
    - my_pattern_trans: Generate a regular expression from input data and parameters.
"""

from inspect import getargvalues, currentframe, stack
from mymetal.universal.check.checkinput import check_input
from mymetal.universal.print.printafter import print_after_blank

def my_pattern_trans(match_type, char = ''):
    """
    generate the regular expression\n
    char to regular expression\n
    """

    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    mypattern = ''
    if match_type:
        if match_type == 'int':
            mypattern = r'([+-]?\d+)'
        elif match_type == 'float':
            mypattern = r'([+-]?\d*\.?\d+|[+-]?\d+\.)'
        elif match_type == 'scientific notation':
            mypattern = r'([+-]?\d*\.?\d+[eE][-+]?\d+)'
        elif match_type == 'bool':
            mypattern = r'\b(true|false|True|False|1|0|yes|no|Yes|No|T|F|t|f)\b'
        elif match_type == 'char':
            mypattern = r'(\w+)\s+'
        elif match_type == 'full line':
            mypattern = f'{char}'
        else:
            mypattern = match_type
    else:
        print_after_blank('match_type', calling_function, [])
    return mypattern