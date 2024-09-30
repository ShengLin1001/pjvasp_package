SOURCE_FOLDER="1.000"

# 需要创建的拉伸文件夹列表
STRETCH_FACTORS=(0.960 0.980 1.020 1.040)

# 平衡态时的晶格常数，假设POSCAR第二行初始值为这个常数
EQUILIBRIUM_LATTICE_CONSTANT=2.8485  # 请根据你的实际晶格常数调整

# 要复制的文件列表
FILES_TO_COPY=( "INCAR" "POSCAR" "KPOINTS" "POTCAR" "sub.544.sh")

# 遍历每个拉伸因子
for FACTOR in "${STRETCH_FACTORS[@]}"
do
    # 创建新的文件夹
    NEW_FOLDER=$(printf "%.3f" $FACTOR)
    echo "创建文件夹: $NEW_FOLDER"
    mkdir -p "$NEW_FOLDER"

    # 将文件从源文件夹复制到新文件夹
    for FILE in "${FILES_TO_COPY[@]}"
    do
        echo "复制 $FILE 到 $NEW_FOLDER 文件夹"
        cp "$SOURCE_FOLDER/$FILE" "$NEW_FOLDER/"
    done

    # 修改POSCAR文件中的第二行 (晶格常数)
    POSCAR_FILE="$NEW_FOLDER/POSCAR"
    
    if [ -f "$POSCAR_FILE" ]; then
        # 计算新的晶格常数
        NEW_LATTICE_CONSTANT=$(echo "$EQUILIBRIUM_LATTICE_CONSTANT * $FACTOR" | bc -l)
        
        # 修改POSCAR文件第二行
        echo "修改 $POSCAR_FILE 第二行为: $NEW_LATTICE_CONSTANT"
        sed -i "2s/.*/$NEW_LATTICE_CONSTANT/" "$POSCAR_FILE"
    else
        echo "文件 $POSCAR_FILE 不存在，跳过修改"
    fi

    cd "$NEW_FOLDER"
    sbatch sub.544.sh
    cd ../
done

echo "脚本执行完毕。"

