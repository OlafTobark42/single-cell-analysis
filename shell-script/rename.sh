#!/bin/bash

# 创建一个新的文件夹用于存放按样本分组后的文件夹
mkdir -p output

# 遍历当前目录中的每个文件
for file in *; do
    if [[ -f $file ]]; then
        # 提取样本名称，假设文件名以样本名称开头，直到最后一个下划线
        sample=${file%_*}

        # 创建用于存放样本文件的文件夹
        mkdir -p "output/$sample"

        # 根据文件类型重命名并移动文件
        if [[ $file == *"matrix.mtx.gz" ]]; then
            mv "$file" "output/$sample/matrix.mtx.gz"
        elif [[ $file == *"barcodes.tsv.gz" ]]; then
            mv "$file" "output/$sample/barcodes.tsv.gz"
        elif [[ $file == *"features.tsv.gz" ]]; then
            mv "$file" "output/$sample/features.tsv.gz"
        fi
    fi
done

# 将output文件夹中的所有内容移回上一级目录
mv output/* ./

# 删除output文件夹
rmdir output

