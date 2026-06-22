from mymetal.universal.print.print import fail
import string
from pathlib import Path

# Check part
# True means pass, False means fail
def check_positive_int(value, name):
    value = str(value)
    _DIGITS = set(string.digits) # {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'}
    # 确认没负号，小数点，或其他非数字字符；同时确保第一个字符不是 '0'（除非值本身是 '0'，但这里不允许 '0' 作为正整数）。
    if value and value[0] in "123456789" and all(c in _DIGITS for c in value):
        return int(value)
    fail(f"%s must be a positive integer: %s" % (name, value if value else "<unset>"))

def check_none(value, name):
    if value is not None:
        return value
    fail(f"%s 参数不能为空" % name)

def check_basename(value):
    if "/" in value or "\\" in value:
        fail(f"目录名称 {value} 不能包含路径分隔符")
    return value

def check_absolute_path(value):
    if not Path(value).is_absolute():
        fail(f"必须是绝对路径: {value}")
    return value