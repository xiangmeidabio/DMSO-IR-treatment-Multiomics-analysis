#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
数据预处理脚本
包含数据加载、特征筛选、标准化和训练测试集划分
"""

import pandas as pd
import numpy as np
import sys
import os
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
import pickle

def main():
    print("--- 开始数据预处理 ---")
    
    # --- 1. 配置文件路径（使用相对路径）---
    # 相对于当前脚本所在目录的路径
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)  # 上一级目录
    
    # 数据文件路径
    all_data_file = os.path.join(project_root, "data", "All_Groups_Combined_transposed.csv")
    gene_list_file = os.path.join(project_root, "data", "key_gene.csv")
    
    # 输出目录
    output_dir = os.path.join(project_root, "results")
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"数据文件路径: {all_data_file}")
    print(f"基因列表文件路径: {gene_list_file}")
    print(f"输出目录: {output_dir}")
    
    # --- 2. 加载100个基因列表 ---
    try:
        gene_df = pd.read_csv(gene_list_file, header=None, index_col=0)
        selected_gene_list = gene_df.index.tolist()
        print(f"成功加载 {len(selected_gene_list)} 个统计学筛选的基因。")
    except FileNotFoundError:
        print(f"错误: 找不到基因列表文件: {gene_list_file}")
        sys.exit(1)
    
    # --- 3. 加载所有四组数据 ---
    try:
        data_df = pd.read_csv(all_data_file, index_col=0)
        print(f"成功加载四组数据，总形状: {data_df.shape}")
    except FileNotFoundError:
        print(f"错误: 找不到四组数据文件: {all_data_file}")
        sys.exit(1)
    
    # --- 4. 分离 X 和 y ---
    if 'label' not in data_df.columns:
        print(f"错误: 在 {all_data_file} 中未找到 'label' 列。")
        sys.exit(1)
        
    X_4class = data_df.drop('label', axis=1)
    y_4class = data_df['label']
    
    # 确保 y_4class 是字符串
    y_4class = y_4class.astype(str)
    
    # --- 5. 检查并筛选特征 ---
    available_genes = [gene for gene in selected_gene_list if gene in X_4class.columns]
    missing_genes = [gene for gene in selected_gene_list if gene not in X_4class.columns]
    
    print(f"在数据文件中找到了 {len(available_genes)} / {len(selected_gene_list)} 个基因。")
    
    if missing_genes:
        print(f"警告: 未在数据中找到 {len(missing_genes)} 个基因。")
        print(f"  (例如: {missing_genes[:5]}...)")
        if not available_genes:
            print("错误: 0个基因被匹配上，请检查基因名是否一致 (大小写、空格等)。")
            sys.exit(1)
    
    # 根据找到的基因列表筛选 X
    X_filtered = X_4class[available_genes]
    
    # --- 6. 先划分数据（未标准化） ---
    print("正在划分训练集、验证集和测试集（在标准化前）...")

    # 先分出测试集 (10%) — 在未缩放的数据上进行划分以避免数据泄露
    X_temp, X_test, y_temp, y_test, names_temp, names_test = train_test_split(
        X_filtered.values, y_4class, data_df.index,
        test_size=0.1, random_state=42, stratify=y_4class
    )

    # 再从剩余数据中分出验证集 (20% of remaining = 18% of total)
    X_train, X_val, y_train, y_val, names_train, names_val = train_test_split(
        X_temp, y_temp, names_temp,
        test_size=0.2222, random_state=42, stratify=y_temp  # 0.2222 * 0.9 ≈ 0.2
    )

    # --- 7. 标准化：仅在训练集上 fit，然后 transform 三个集合 ---
    print("正在执行 StandardScaler 标准化 (在训练集上 fit)...")
    scaler = StandardScaler()

    # 在训练集上 fit
    scaler.fit(X_train)

    # transform 三个集合
    X_train = scaler.transform(X_train)
    X_val = scaler.transform(X_val)
    X_test = scaler.transform(X_test)
    
    print(f"\n数据集划分结果:")
    total = X_train.shape[0] + X_val.shape[0] + X_test.shape[0]
    print(f"训练集: {X_train.shape[0]} 样本 ({X_train.shape[0] / total * 100:.1f}%)")
    print(f"验证集: {X_val.shape[0]} 样本 ({X_val.shape[0] / total * 100:.1f}%)")
    print(f"测试集: {X_test.shape[0]} 样本 ({X_test.shape[0] / total * 100:.1f}%)")
    
    # --- 9. 保存处理后的数据 ---
    print("正在保存处理后的数据...")
    
    # 保存训练集
    train_data = {
        'X_train': X_train,
        'y_train': y_train,
        'names_train': names_train
    }
    
    # 保存验证集
    val_data = {
        'X_val': X_val,
        'y_val': y_val,
        'names_val': names_val
    }
    
    # 保存测试集
    test_data = {
        'X_test': X_test,
        'y_test': y_test,
        'names_test': names_test
    }
    
    # 保存其他重要信息
    metadata = {
        'available_genes': available_genes,
        'scaler': scaler,
        'feature_names': available_genes,
        'class_distribution': y_4class.value_counts().to_dict()
    }
    
    # 使用pickle保存
    with open(os.path.join(output_dir, 'train_data.pkl'), 'wb') as f:
        pickle.dump(train_data, f)
    
    with open(os.path.join(output_dir, 'val_data.pkl'), 'wb') as f:
        pickle.dump(val_data, f)
        
    with open(os.path.join(output_dir, 'test_data.pkl'), 'wb') as f:
        pickle.dump(test_data, f)
        
    with open(os.path.join(output_dir, 'metadata.pkl'), 'wb') as f:
        pickle.dump(metadata, f)
    
    # 另外保存为CSV格式（可选）
    train_df = pd.DataFrame(X_train, columns=available_genes)
    train_df.insert(0, 'SampleID', names_train)
    train_df['label'] = y_train.values
    train_df.to_csv(os.path.join(output_dir, 'train_data.csv'), index=False)
    
    val_df = pd.DataFrame(X_val, columns=available_genes)
    val_df.insert(0, 'SampleID', names_val)
    val_df['label'] = y_val.values
    val_df.to_csv(os.path.join(output_dir, 'val_data.csv'), index=False)
    
    test_df = pd.DataFrame(X_test, columns=available_genes)
    test_df.insert(0, 'SampleID', names_test)
    test_df['label'] = y_test.values
    test_df.to_csv(os.path.join(output_dir, 'test_data.csv'), index=False)
    
    print(f"\n--- 预处理完成 ---")
    # X 已经被划分并分别缩放，打印原始过滤后特征矩阵的形状以表示总体特征维度
    print(f"最终 X (特征) 原始形状 (未缩放): {X_filtered.shape}")
    print(f"最终 y (标签) 形状: {y_4class.shape}")
    print(f"特征数量: {len(available_genes)}")
    print(f"类别分布: {y_4class.value_counts().to_dict()}")
    
    print(f"\n数据已保存到:")
    print(f"  - {os.path.join(output_dir, 'train_data.pkl')} (训练集)")
    print(f"  - {os.path.join(output_dir, 'val_data.pkl')} (验证集)")
    print(f"  - {os.path.join(output_dir, 'test_data.pkl')} (测试集)")
    print(f"  - {os.path.join(output_dir, 'metadata.pkl')} (元数据)")
    print(f"  - CSV格式文件也已保存")
    
    return True

if __name__ == "__main__":
    main()
