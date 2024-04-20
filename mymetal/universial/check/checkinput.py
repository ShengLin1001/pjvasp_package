from inspect import stack

##############################################################
############# check input 
##############################################################

def check_input(args = [], check_type = 'none'):
    """
    check_input_blank function is not correct, I think it's useless\n
    first line  - check_input_none\n
    second line - check_input_blank\n
    """

    def check_input_none(args = [], calling_function = ''):
        # print(calling_function)
        for arg_name, arg_value in args.items():
            if arg_value is None:
                raise ValueError(f"In '{calling_function}': Argument '{arg_name}' is missing or None")
        
    # Get the name of the calling function
    calling_function = stack()[1].function
    if check_type == 'none':
        check_input_none(args, calling_function )
    else:
        print(f"this check type {check_type} is not supported")