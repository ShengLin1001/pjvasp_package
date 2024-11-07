"""
construct submodule

This submodule provides functions for generating content for post-processing tasks related to materials simulations designed for mymetal.post subpackage.

Functions:
    - create_content: Creates the content required for post-processing tasks.
"""

from inspect import getargvalues, currentframe, stack
from mymetal.universial.check.checkinput import check_input
from mymetal.universial.data.dataadjust import myjust
from mymetal.universial.data.datachange import list_to_char

##############################################################
############# construct the content
##############################################################

def create_content(list = [], tag_split = '   ', tag_end = '', tag_begin = '', left = True, num = 10):
    """
    have a list of char, using myjust() to formatted them.\n
    output a char to write in a specified file\n
    """

    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    mycontent_char = ''
    mycontent = []
    mycontent.append(tag_begin)
    for i in list:
        formatted_i = myjust(left, num, i)
        content_i = f'{formatted_i}' + tag_split
        mycontent.append(content_i)
    mycontent.append(tag_end)
    mycontent_char = list_to_char(mycontent, char_tag = '')
    return mycontent_char
