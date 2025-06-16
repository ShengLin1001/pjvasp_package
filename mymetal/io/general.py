"""
general module

This module provides a function to read delimited files into a pandas DataFrame with flexible header handling.

Functions:
    - general_read: Reads a delimited file into a DataFrame, allowing for custom headers and separators.

"""
import pandas as pd
from typing import Optional, List

def general_read(filepath: str = None, has_header: bool = True, header_names: Optional[List[str]] = None, sep: str = r"\s+", comment_char='#') -> pd.DataFrame:
    """Reads a delimited file into a DataFrame with flexible header handling.

    Args:
        filepath: Path to the input file. If None, may trigger default behavior.
        has_header: If True, uses first line as column names. Defaults to True.
        header_names: Custom column names when has_header=False. Ignored if has_header=True.
        sep: Delimiter regex pattern. Defaults to any whitespace (r"\s+").
        comment_char: Lines starting with this character are skipped. Defaults to '#'.

    Returns:
        pd.DataFrame: Parsed data with appropriate column headers.

    Note:
        Uses Python engine for regex separator support. For large files, consider specifying engine='c'.
    """
    if has_header:
        df = pd.read_csv(filepath,  
                         sep=sep,
                         header=0,
                         comment=comment_char,
                         engine='python')
    else:
        df = pd.read_csv(filepath, 
                         sep=sep,
                         header=None,
                         names=header_names,
                         comment=comment_char,
                         engine='python')
    return df
