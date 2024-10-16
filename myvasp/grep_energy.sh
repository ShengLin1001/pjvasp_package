#!/bin/bash

# ./script.sh /path/to/output.txt 2 3
# ./script.sh -m 1 -M 1       ##当前目录检索
# 设置默认值
output_file="./energy.txt"
mindepth=1
maxdepth=2

# 处理命令行选项
while getopts "o:m:M:" opt; do
  case $opt in
    o) output_file="$OPTARG" ;;  # 指定输出文件
    m) mindepth="$OPTARG" ;;     # 指定最小深度
    M) maxdepth="$OPTARG" ;;     # 指定最大深度
    \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
  esac
done

# 删除旧的输出文件
rm -f $output_file

# 查找目录并提取能量
find ./ -mindepth $mindepth -maxdepth $maxdepth -type d | while read dir; do
    # 指定OUTCAR路径，获取能量值
    energy=$(grep '  without' "$dir/OUTCAR" | tail -n 1 | awk '{print $7}')
    # 将结果写入输出文件
    echo "$dir  $energy" >> "$output_file"
done

echo "Energy extraction complete. Results saved in $output_file"

