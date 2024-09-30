#!/bin/bash

# 定义源文件
sub_script="sub.544.sh"

# 递归查找所有三级目录
find ./ -mindepth 1  -type d | while read dir; do
    # 复制INCAR和KPOINTS文件到每个三级目录
    cp "$sub_script" "$dir"
    echo "$dir"
    # 切换到该目录并提交任务
    cd "$dir"
    echo 103 |vaspkit
    sbatch "$sub_script"
    cd - > /dev/null
done

