"""
general module

This module provides a function to read delimited files into a pandas DataFrame with flexible header handling.

Functions:
    - general_read: Reads a delimited file into a DataFrame, allowing for custom headers and separators.
    - general_write: Writes a DataFrame to a formatted plain text file.

Change Log:
    - 2025.10.16: Added general_write function to write DataFrame to a formatted text file.
                  Updated general_read to use index_col and header_row parameters.

"""
import pandas as pd
from typing import Optional, List
import numpy as np

# read a delimited file into a DataFrame with flexible header handling
def general_read(filepath: str = None, has_header: bool = True, header_names: Optional[List[str]] = None, sep: str = r"\s+", comment_char='#',
                 index_col: int = None, header_row: int = 0) -> pd.DataFrame:
    """Reads a delimited file into a DataFrame with flexible header handling.

    Args:
        filepath: Path to the input file. If None, may trigger default behavior.
        has_header: If True, uses first line as column names. Defaults to True. header_row specifies which row to use as header.
        header_names: Custom column names when has_header=False. Ignored if has_header=True.
        sep: Delimiter regex pattern. Defaults to any whitespace (r"\s+").
        comment_char: Lines starting with this character are skipped. Defaults to '#'.
        index_col: Column to set as index. Defaults to None (no index).
        header_row: Row number to use as header when has_header=True. Defaults to 0 (first line).

    Returns:
        pd.DataFrame: Parsed data with appropriate column headers.

    Note:
        - Uses Python engine for regex separator support. For large files, consider specifying engine='c'.
        - has_header = True, index_col = None, header_row = 0 (default) is the most common case.
    """

    if has_header:
        df = pd.read_csv(filepath,  
                         sep=sep,
                         header=header_row,
                         comment=comment_char,
                         engine='python',
                         index_col=index_col)
    else:
        df = pd.read_csv(filepath, 
                         sep=sep,
                         header=None,
                         names=header_names,
                         comment=comment_char,
                         engine='python',
                         index_col=index_col)
    return df

# write dataframe
def general_write(filename: str = None, dfc: pd.DataFrame = None, int_format: str = '>16d', str_format: str = '>16s',
                        bool_format: str = ':>16', float_format: str = '16.10f', 
                        if_write_col_num: bool = False,
                        if_write_row_num: bool = False):
    """Write a DataFrame or array to a formatted plain text file.

    Args:
        filename (str): Output file path.
        dfc (pd.DataFrame or array-like): Data to write.
        int_format (str): Format string for integer columns (default '>16d').
        str_format (str): Format string for string columns (default '>16s').
        bool_format (str): Format string for boolean columns (default ':>16', outputs 1/0).
        float_format (str): Format string for float columns (default '16.10f').
        if_write_col_num (bool): Whether to write column headers (default False).
        if_write_row_num (bool): Whether to write row indices (default False).

    Notes:
        - > right align.
        - < left align.
        - ^ center align.
        - Non-standard types (list, dict) are printed as-is without formatting.
    """
    if type(dfc) != pd.DataFrame:
        dfc = pd.DataFrame(dfc)

    df = dfc.copy()

    for col in df.columns:
        ctype = get_col_type(df, col)
        if ctype == 'int':
            df[col] = df[col].apply(lambda x: f'{x:{int_format}}')
        elif ctype == 'float':
            df[col] = df[col].apply(lambda x: f'{x:{float_format}}')
        elif ctype == 'bool':
            # True => 1, False => 0
            df[col] = df[col].apply(lambda x: f'{int(x):{bool_format}}')
        elif ctype == 'str':
            df[col] = df[col].apply(lambda x: f'{x:{str_format}}')
        else:  # list, dict or other types
            print(f'Column {col} is not formatted.')

    # formatting row index
    if if_write_row_num == False:
        df.index = ['']*len(df)

    with open(filename, 'w') as f:
        f.write(df.to_string(index=True, header= True if if_write_col_num else False))


def get_col_type(df, col):
    dtype = df[col].dtype
    if pd.api.types.is_integer_dtype(dtype):
        return 'int'
    elif pd.api.types.is_float_dtype(dtype):
        return 'float'
    elif pd.api.types.is_bool_dtype(dtype):
        return 'bool'
    elif pd.api.types.is_object_dtype(dtype):
        return 'str'
    else:
        return 'other'