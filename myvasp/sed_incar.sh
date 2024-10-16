#!/bin/bash

# bash sed_incar.sh ISTART 1 INCAR
# 第一个参数：匹配的字符（如NSW）
match_string="$1"

# 第二个参数：替换成的字符串（如50）
replacement_value="$2"

# 第三个参数：文件名
filename="${3:-INCAR}"

# 使用sed命令查找包含匹配字符的行并替换为匹配字符=替换值
sed -i "/$match_string/c\\   $match_string=$replacement_value" "$filename"

echo "Replaced line containing '$match_string' with '$match_string=$replacement_value' in $filename"

