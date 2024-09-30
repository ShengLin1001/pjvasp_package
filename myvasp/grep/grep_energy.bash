#!/bin/bash

output_file="./energy_fcc.txt"  # 指定输出文件的位置
rm $output_file

find ./ -mindepth 1 -type d | while read dir; do
    #echo "$dir"
    # 指定OUTCAR路径，获取能量值
    energy=$(grep '  without' "$dir/OUTCAR" | tail -n 1 | awk '{print $7}')
    # 提取 energy  
    echo "$dir  $energy" >> "$output_file"
done

