#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
数据预处理脚本 - 仅使用9个关键保护基因

此脚本执行以下操作：
1. 加载原始的合并数据集。
2. 筛选出指定的9个关键基因作为特征。
3. 将数据划分为训练集、验证集和测试集。
4. 对特征进行标准化处理。
5. 将处理后的数据（训练、验证、测试集和元数据）保存到一个新的子目录 'results/key_genes_subset/' 中，
   以便后续模型训练使用，同时避免覆盖原始的100基因数据集。
"""

import pandas as pd
import numpy as np
import os
import pickle
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, LabelEncoder
import warnings

warnings.filterwarnings('ignore', category=FutureWarning)

def main():
    """主执行函数"""
    print("="*60)
    print("--- 数据预处理 (仅9个关键基因) ---")
    print("="*60)

    # --- 关键基因列表 ---
    key_genes = [
        "Rps2", "Uqcrb", "Gapdh", "Impa2", "Ptma", "Prtn3", "Ramp1", "Snrpf", "Cenpw"
    ]
    print(f"🧬 将使用以下 {len(key_genes)} 个关键基因作为特征:\n{key_genes}")

    # --- 路径设置 ---
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    data_dir = os.path.join(project_root, "data")
    output_dir = script_dir  # 输出到当前脚本所在目录
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"\n📁 输入数据目录: {data_dir}")
    print(f"📁 输出数据将保存在: {output_dir}")

    # --- 1. 加载数据 ---
    print("\n🔄 [Step 1] 正在加载原始合并数据...")
    try:
        # 假设原始数据文件名为 'merged_data.csv'
        data_path = os.path.join(data_dir, 'All_Groups_Combined_transposed.csv')
        df = pd.read_csv(data_path)
        print(f"✅ 原始数据加载成功! 形状: {df.shape}")
    except FileNotFoundError:
        print(f"❌ 错误: 找不到原始数据文件 '{data_path}'。")
        print("请确保 'data/All_Groups_Combined_transposed.csv' 文件存在。")
        return

    # --- 2. 准备特征和标签 ---
    print("\n🔬 [Step 2] 正在准备特征和标签...")
    
    # 检查关键基因是否存在于数据中
    missing_genes = [gene for gene in key_genes if gene not in df.columns]
    if missing_genes:
        print(f"❌ 错误: 以下关键基因在数据中不存在: {missing_genes}")
        return
        
    X = df[key_genes]
    y = df['label']
    sample_names = df['SampleID']

    # 标签编码
    le = LabelEncoder()
    y_encoded = le.fit_transform(y)
    
    print(f"✅ 特征矩阵 (X) 形状: {X.shape}")
    print(f"✅ 标签 (y) 形状: {y_encoded.shape}")
    print(f"✅ 类别分布:\n{y.value_counts()}")

    # --- 3. 划分数据集 ---
    print("\n🔪 [Step 3] 正在划分训练、验证和测试集...")
    # 仿照旧脚本，采用 70% 训练, 20% 验证, 10% 测试 的划分
    # 先分出 10% 作为测试集
    X_temp, X_test, y_temp, y_test, names_temp, names_test = train_test_split(
        X, y_encoded, sample_names, test_size=0.1, random_state=42, stratify=y_encoded
    )
    # 再从剩余的 90% 中分出验证集 (22.22% * 90% ≈ 20% of total)
    X_train, X_val, y_train, y_val, names_train, names_val = train_test_split(
        X_temp, y_temp, names_temp, test_size=0.2222, random_state=42, stratify=y_temp
    )
    
    total_samples = len(X)
    print(f"  - 训练集: {X_train.shape[0]} 个样本 ({X_train.shape[0]/total_samples:.1%})")
    print(f"  - 验证集: {X_val.shape[0]} 个样本 ({X_val.shape[0]/total_samples:.1%})")
    print(f"  - 测试集: {X_test.shape[0]} 个样本 ({X_test.shape[0]/total_samples:.1%})")

    # --- 4. 特征标准化 ---
    print("\n📏 [Step 4] 正在对特征进行标准化...")
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_val_scaled = scaler.transform(X_val)
    X_test_scaled = scaler.transform(X_test)
    print("✅ 标准化完成!")

    # --- 5. 保存处理后的数据 ---
    print("\n💾 [Step 5] 正在保存处理后的数据...")
    
    # 封装数据
    train_data = {'X_train': X_train_scaled, 'y_train': pd.Series(y_train), 'names_train': names_train.values}
    val_data = {'X_val': X_val_scaled, 'y_val': pd.Series(y_val), 'names_val': names_val.values}
    test_data = {'X_test': X_test_scaled, 'y_test': pd.Series(y_test), 'names_test': names_test.values}
    
    # 元数据
    metadata = {
        'feature_names': key_genes,
        'class_distribution': y.value_counts().to_dict(),
        'label_encoder': le,
        'scaler': scaler
    }

    # 保存到文件
    with open(os.path.join(output_dir, 'train_data.pkl'), 'wb') as f:
        pickle.dump(train_data, f)
    with open(os.path.join(output_dir, 'val_data.pkl'), 'wb') as f:
        pickle.dump(val_data, f)
    with open(os.path.join(output_dir, 'test_data.pkl'), 'wb') as f:
        pickle.dump(test_data, f)
    with open(os.path.join(output_dir, 'metadata.pkl'), 'wb') as f:
        pickle.dump(metadata, f)
        
    print(f"✅ 所有数据已成功保存至 '{output_dir}'")
    print("\n🎉 数据预处理流程执行完毕!")

if __name__ == "__main__":
    main()
