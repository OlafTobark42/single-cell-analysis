#!/bin/bash

# 创建一个新的文件夹用于存放按样本分组后的文件夹
mkdir -p output

# 遍历当前目录中的每个文件
for file in *; do
    if [[ -f $file ]]; then
        # 提取样本名称，假设文件名以样本名称结尾
        sample=${file%%_counts.csv.gz}
        sample=${sample%%_metadata.csv.gz}

        # 创建用于存放样本文件的文件夹
        mkdir -p "output/$sample"

        # 将文件移动到对应的样本文件夹中
        mv "$file" "output/$sample/"

    fi
done

