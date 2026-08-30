from mymetal.universal.print.print import fail
import string
from pathlib import Path

# Check part
# True means pass, False means fail
def check_positive_int(value, name):
    """Validate that *value* is a positive integer (string of digits, no leading zero).

    Args:
        value: The value to check (will be cast to ``str``).
        name (str): Parameter name for the error message.

    Returns:
        int: The validated positive integer.

    Raises:
        SystemExit: If *value* is not a positive integer.
    """
    value = str(value)
    _DIGITS = set(string.digits) # {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'}
    # 确认没负号，小数点，或其他非数字字符；同时确保第一个字符不是 '0'（除非值本身是 '0'，但这里不允许 '0' 作为正整数）。
    if value and value[0] in "123456789" and all(c in _DIGITS for c in value):
        return int(value)
    fail(f"%s must be a positive integer: %s" % (name, value if value else "<unset>"))

def check_none(value, name):
    """Validate that *value* is not ``None``; return *value* on success.

    Args:
        value: The value to check.
        name (str): Parameter name for the error message.

    Returns:
        The validated *value*.

    Raises:
        SystemExit: If *value* is ``None``.
    """
    if value is not None:
        return value
    fail(f"%s 参数不能为空" % name)

def check_basename(value):
    """Validate that *value* is a basename (no path separators).

    Args:
        value (str): The string to check.

    Returns:
        str: The validated basename.

    Raises:
        SystemExit: If *value* contains ``/`` or ``\\``.
    """
    if "/" in value or "\\" in value:
        fail(f"目录名称 {value} 不能包含路径分隔符")
    return value

def check_absolute_path(value):
    """Validate that *value* is an absolute path.

    Uses ``pathlib.Path.is_absolute()``. On Windows, a path like
    ``/work`` is NOT absolute (it needs a drive letter, e.g. ``C:/work``).
    On Linux, ``/work`` is absolute.

    Args:
        value (str | Path): The path to check.

    Returns:
        The validated *value*.

    Raises:
        SystemExit: If *value* is not absolute.
    """
    if not Path(value).is_absolute():
        fail(f"必须是绝对路径: {value}")
    return value