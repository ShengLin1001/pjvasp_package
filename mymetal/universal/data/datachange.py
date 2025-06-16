"""
datachange submodule

This submodule provides functions for changing data formats. It includes functions for converting lists to strings, strings to dictionaries, 
and dictionaries to strings. These functions are designed to streamline common tasks in materials science simulations and data handling.

Functions:
    - list_to_char: Convert a list of characters or dictionaries to a string.
    - char_to_dic: Convert a string to a dictionary.
    - dic_to_char: Convert a dictionary to a string.
"""

from inspect import getargvalues, currentframe, stack
from mymetal.universal.check.checkinput import check_input
from mymetal.universal.print.printafter import print_after_blank

##############################################################
############# change type
##############################################################

def list_to_char(mycontent = [], list_type = "char", char_tag = '   ',
                 dic_tag = '\n', key_value_interval = ': ', item_interval = ', '):
    """
    consisted of two function\n
    list_of_char_to_char to output the list of char to file\n
    list_of_dic_to_char to output the list of dic to file\n
    """

    # mycontent1 = ['item1', 'item2', 'item3']
    # joined_string = 'item1 item2 item3'
    # joined_string = ' '.join(mycontent1)

    def list_of_char_to_char(mycontent = [], tag = '   '):
        """for 1-D list of char, change it to char to optput to file"""
        
        mycontent_char = ''
        for content in mycontent:
            mycontent_char = mycontent_char  + f'{content}' + tag
        return mycontent_char
    
    def list_of_dic_to_char(list_of_dicts, tag = '\n', key_value_interval = ': ', item_interval = ', '):
        """
        for 1-D list of dic, change it to char\n
                  [{"name":"1","gender":"female"}, {"name":"2","gender":"male"}] \n
         =>>      "name: 1,gender: female\n
                   name: 2, gender: male"\n
        """

        formatted_dicts = [item_interval.join([f"{key}{key_value_interval}{value}" for key, value in d.items()]) + tag for d in list_of_dicts]
        return ''.join(formatted_dicts)
    
    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)

    temp = ''
    if list_type == "char":
        temp = list_of_char_to_char(mycontent, char_tag )
    elif list_type == "dic":
        temp = list_of_dic_to_char(mycontent, dic_tag, key_value_interval, item_interval)
    else:
        print(f"the list_type wanted to transformed is not supported, please fix it! - {calling_function}")
    return temp

def char_to_dic(mycontent_char):
    """"name:n1 gender:male" =>> {"name":"n1", "gender":"male"}"""
    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    mycontent_dic = {}
    # Split the input line by spaces to separate parameters
    mycontent_char_copy = mycontent_char.split()
    # Create a dictionary to store the parsed parameters
    # Iterate through the parameters and parse them into key-value pairs
    for char in mycontent_char_copy:
        dic_name, dic_value = char.split(':')
        mycontent_dic[dic_name] = dic_value
    return mycontent_dic

def dic_to_char(mycontent_dic={}, tag = '   '):
    """f"name:n1 {tag} gender:male" <<= {"name":"n1", "gender":"male"}"""
    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    # Initialize an empty string to store the formatted content
    mycontent_char = ''
    if mycontent_dic:
        # Iterate over key-value pairs in the dictionary
        for key, value in mycontent_dic.items():
            # Append each key-value pair to the formatted string
            mycontent_char = mycontent_char + f'{key}:{value}' + tag
        # Remove the trailing space at the end of the string
        mycontent_char += '\n'
    else:
        print_after_blank('mycontent_dic', calling_function, [])
    return mycontent_char