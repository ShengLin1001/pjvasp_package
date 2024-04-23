#!/bin/bash

# 定义源代码目录和文档目录
SRC_DIR="../mymetal"
DOC_DIR="docs/source/mymetal"

# 创建基础文档目录
mkdir -p $DOC_DIR

# 递归遍历源代码目录中的每个子目录
function generate_docs {
    local current_dir=$1
    # 生成.rst文件
    sphinx-apidoc -o "$DOC_DIR/${current_dir#$SRC_DIR/}" "$current_dir" -f

    # 遍历子目录
    for subdir in "$current_dir"/*; do
        if [ -d "$subdir" ]; then
            generate_docs "$subdir"
        fi
    done
}

# 调用函数，开始生成文档
generate_docs $SRC_DIR

echo "Documentation generation complete."
