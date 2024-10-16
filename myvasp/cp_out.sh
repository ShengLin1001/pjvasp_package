#!/bin/bash

# bash cp_out.sh -t 1 -o slurm*.out
# 设置默认的目标文件夹、out文件名称和默认文件夹名称
target_folder="./backup"
out_file="vasp.out"
default_folder="."
default_files=("CONTCAR" "POSCAR" "OUTCAR" "XDATCAR" "INCAR" "KPOINTS")

# 处理命令行选项
while getopts "t:o:d:f:" opt; do
  case $opt in
    t)  target_folder="$OPTARG" ;;  # 目标文件夹
    o)  out_file="$OPTARG" ;;       # out 文件
    f)  default_folder="$OPTARG" ;; # 默认文件夹
    d)  default_files="$OPTARG" ;;
    \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
  esac
done

# 如果目标文件夹不存在，则创建它
if [ ! -d "$target_folder" ]; then
    mkdir -p "$target_folder"
fi

# 检查 out 文件是否存在于默认文件夹中
if [ -e "$default_folder/$out_file" ]; then
    # 复制 out 文件到目标文件夹
    cp "$default_folder/$out_file" "$target_folder/"
    echo "Copied $out_file from $default_folder to $target_folder"
else
    echo "Error: $out_file does not exist in $default_folder"
fi

# 复制 default_files 中的文件
for file in "${default_files[@]}"; do
  if [ -e "$default_folder/$file" ]; then
    cp "$default_folder/$file" "$target_folder/"
    echo "Copied $file from $default_folder to $target_folder"
  else
    echo "Warning: $file does not exist in $default_folder, skipping..."
  fi
done


