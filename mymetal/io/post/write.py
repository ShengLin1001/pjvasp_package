""" 
write submodule

This submodule provides functions for writing content to files for post-processing tasks related to materials simulations designed for mymetal.post subpackage.

Functions:
    - write_content_to_file: Writes content to a file.
"""

from inspect import getargvalues, currentframe, stack
from mymetal.universal.check.checkinput import check_input

##############################################################
#############  output 
##############################################################

def write_content_to_file(content_char, path, write_type='a'):
    """
    write content_char to file controled by write_type.\n
    for example, a-append, w-write overlap\n
    """
    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    try:
        with open(path, write_type) as post_file:
            post_file.write(content_char)
            #print(f"Content written to {path}")
    except FileNotFoundError:
        print(f"i can't open the path: '{path}'! - {calling_function}")