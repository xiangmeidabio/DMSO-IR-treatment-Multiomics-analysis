#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Gradient Boosting模型训练脚本 - 仅使用9个关键基因

此脚本基于 '06_gradient_boosting_model.py' 修改而来，专门用于在
仅包含9个关键基因的数据子集上训练和评估Gradient Boosting模型。

主要改动：
- 从脚本所在目录 ('7-best_model/') 加载由 '06_data_preprocessing_key_genes.py' 生成的数据。
- 将所有输出（模型、结果、报告、图表）保存到子目录 'key_genes_model_outputs/' 中，
  以避免与100个基因训练的原始模型结果混淆。
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pickle
import sys
from sklearn.ensemble import GradientBoostingClassifier
from skopt.space import Real, Integer
from skopt import BayesSearchCV
from sklearn.metrics import accuracy_score, roc_auc_score, f1_score, confusion_matrix, classification_report
import traceback

def load_preprocessed_data(data_dir):
    """从指定目录加载预处理好的数据"""
    print(f"正在从 '{data_dir}' 加载预处理好的数据...")
    
    try:
        with open(os.path.join(data_dir, 'train_data.pkl'), 'rb') as f:
            train_data = pickle.load(f)
        with open(os.path.join(data_dir, 'val_data.pkl'), 'rb') as f:
            val_data = pickle.load(f)
        with open(os.path.join(data_dir, 'test_data.pkl'), 'rb') as f:
            test_data = pickle.load(f)
        with open(os.path.join(data_dir, 'metadata.pkl'), 'rb') as f:
            metadata = pickle.load(f)
        
        print("✅ 数据加载成功!")
        print(f"训练集形状: {train_data['X_train'].shape}")
        print(f"验证集形状: {val_data['X_val'].shape}")
        print(f"测试集形状: {test_data['X_test'].shape}")
        print(f"特征数量: {len(metadata['feature_names'])}")
        
        return train_data, val_data, test_data, metadata
        
    except FileNotFoundError as e:
        print(f"❌ 错误: 找不到预处理数据文件: {e}")
        print("请先运行 '06_data_preprocessing_key_genes.py' 生成关键基因的数据。")
        sys.exit(1)

def calculate_metrics(y_true, y_pred, y_proba, le):
    """计算模型评估指标"""
    metrics = {}
    metrics['Accuracy'] = accuracy_score(y_true, y_pred)
    metrics['AUC'] = roc_auc_score(y_true, y_proba, multi_class='ovr')
    metrics['Macro F1'] = f1_score(y_true, y_pred, average='macro')
    metrics['Micro F1'] = f1_score(y_true, y_pred, average='micro')
    metrics['Confusion Matrix'] = confusion_matrix(y_true, y_pred)
    metrics['Classification Report'] = classification_report(
        y_true, y_pred, target_names=le.classes_, output_dict=True
    )
    return metrics

def plot_confusion_matrix(cm, class_names, model_name, save_path=None):
    """绘制混淆矩阵"""
    fig, ax = plt.subplots(figsize=(6, 6))
    sns.heatmap(
        cm,
        annot=True,
        fmt='d',
        cmap='Oranges',
        xticklabels=class_names,
        yticklabels=class_names,
        ax=ax,
        square=True,
        cbar=True,
        cbar_kws={'shrink': 0.6, 'pad': 0.02}
    )
    ax.set_xlabel('Predicted label')
    ax.set_ylabel('True label')
    ax.set_title(model_name, fontsize=12, pad=15)
    plt.tight_layout()

    if save_path:
        base, _ = os.path.splitext(save_path)
        plt.savefig(f"{base}.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{base}.pdf", dpi=300, bbox_inches='tight')
        print(f"混淆矩阵已保存至: {base}.png/pdf")
    plt.show()

def train_gradient_boosting(X_train, y_train):
    """训练Gradient Boosting模型并进行贝叶斯优化"""
    print("\n🚀 开始训练Gradient Boosting模型 (9个关键基因)...")
    
    try:
        model = GradientBoostingClassifier(random_state=42)
        
        params = {
            'n_estimators': Integer(50, 200),
            'max_depth': Integer(3, 10),
            'learning_rate': Real(0.01, 0.3, prior='log-uniform'),
            'min_samples_split': Integer(2, 10),
            'min_samples_leaf': Integer(1, 5),
            'subsample': Real(0.8, 1.0)
        }
        
        print("正在进行贝叶斯超参数优化...")
        bayes_search = BayesSearchCV(
            model, params, n_iter=30, cv=3, scoring='accuracy', random_state=42, n_jobs=4, verbose=1
        )
        
        bayes_search.fit(X_train, y_train)
        best_model = bayes_search.best_estimator_
        
        print(f"✅ 模型训练完成!")
        print(f"最佳参数: {bayes_search.best_params_}")
        print(f"最佳交叉验证得分: {bayes_search.best_score_:.4f}")
        
        return best_model
        
    except Exception as e:
        print(f"❌ 模型训练失败: {str(e)}")
        traceback.print_exc()
        return None

def evaluate_model(model, X_test, y_test, le, output_dir):
    """评估模型性能"""
    print("\n📊 开始模型评估 (9个关键基因)...")
    
    try:
        y_pred = model.predict(X_test)
        y_proba = model.predict_proba(X_test)
        
        metrics = calculate_metrics(y_test, y_pred, y_proba, le)
        
        print(f"\n🎯 测试集性能指标:")
        print(f"  准确率 (Accuracy): {metrics['Accuracy']:.4f}")
        print(f"  多类别AUC: {metrics['AUC']:.4f}")
        print(f"  宏平均F1: {metrics['Macro F1']:.4f}")
        print(f"  微平均F1: {metrics['Micro F1']:.4f}")
        
        cm_path = os.path.join(output_dir, 'confusion_matrix_gb_key_genes')
        plot_confusion_matrix(
            metrics['Confusion Matrix'], 
            le.classes_, 
            'GB Confusion Matrix (9 Key Genes)', 
            save_path=cm_path
        )
        
        return metrics, y_pred
        
    except Exception as e:
        print(f"❌ 模型评估失败: {str(e)}")
        traceback.print_exc()
        return None, None

def save_results(X_test, y_test, y_pred, names_test, feature_names, metrics, model, output_dir):
    """保存模型和评估结果"""
    print("\n💾 正在保存结果 (9个关键基因)...")
    
    try:
        # 1. 保存测试集预测结果
        test_results_df = pd.DataFrame(X_test, columns=feature_names)
        test_results_df.insert(0, 'SampleID', names_test)
        test_results_df['True_label'] = y_test
        test_results_df['Predicted_label'] = y_pred
        
        results_path = os.path.join(output_dir, 'gradient_boosting_test_results_key_genes.csv')
        test_results_df.to_csv(results_path, index=False, encoding='utf-8-sig')
        
        # 2. 保存模型
        model_path = os.path.join(output_dir, 'gradient_boosting_model_key_genes.pkl')
        with open(model_path, 'wb') as f:
            pickle.dump(model, f)
        
        # 3. 保存评估指标
        metrics_path = os.path.join(output_dir, 'gradient_boosting_metrics_key_genes.pkl')
        with open(metrics_path, 'wb') as f:
            pickle.dump(metrics, f)
            
        # 4. 保存详细报告
        report_path = os.path.join(output_dir, 'gradient_boosting_report_key_genes.txt')
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("Gradient Boosting (9个关键基因) 模型评估报告\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"准确率: {metrics['Accuracy']:.4f}\n")
            f.write(f"多类别AUC: {metrics['AUC']:.4f}\n")
            f.write(f"宏平均F1: {metrics['Macro F1']:.4f}\n")
            f.write(f"微平均F1: {metrics['Micro F1']:.4f}\n\n")
            f.write("详细分类报告:\n")
            f.write("-" * 30 + "\n")
            
            # 格式化分类报告
            for class_name, class_metrics in metrics['Classification Report'].items():
                if isinstance(class_metrics, dict):
                    f.write(f"{class_name}:\n")
                    f.write(f"  精确率: {class_metrics.get('precision', 0):.4f}\n")
                    f.write(f"  召回率: {class_metrics.get('recall', 0):.4f}\n")
                    f.write(f"  F1分数: {class_metrics.get('f1-score', 0):.4f}\n")
                    f.write(f"  支持数: {class_metrics.get('support', 0)}\n\n")
            
        print(f"✅ 结果保存完成!")
        print(f"  测试结果: {results_path}")
        print(f"  模型文件: {model_path}")
        print(f"  评估指标: {metrics_path}")
        print(f"  详细报告: {report_path}")

    except Exception as e:
        print(f"❌ 保存结果失败: {str(e)}")
        traceback.print_exc()

def main():
    """主函数"""
    print("\n" + "="*60)
    print("--- 🚀 Gradient Boosting (9个关键基因) ---")
    print("="*60)
    
    # --- 路径设置 ---
    # 脚本所在的目录，即 7-best_model
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # 数据文件和脚本在同一目录下
    data_dir = script_dir
    # 将所有输出保存到一个专用的子目录中
    output_dir = os.path.join(script_dir, "key_genes_model_outputs")
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"\n📁 输入数据目录: {data_dir}")
    print(f"📁 输出结果目录: {output_dir}")

    # 1. 加载数据
    train_data, val_data, test_data, metadata = load_preprocessed_data(data_dir)
    
    X_train, y_train = train_data['X_train'], train_data['y_train']
    X_test, y_test = test_data['X_test'], test_data['y_test']
    names_test = test_data['names_test']
    feature_names = metadata['feature_names']
    le = metadata['label_encoder']
    
    # 2. 训练模型
    model = train_gradient_boosting(X_train, y_train)
    if model is None:
        print("❌ 模型训练失败，程序退出")
        return
    
    # 3. 评估模型
    metrics, y_pred = evaluate_model(model, X_test, y_test, le, output_dir)
    if metrics is None:
        print("❌ 模型评估失败，程序退出")
        return
    
    # 4. 保存结果
    save_results(X_test, y_test, y_pred, names_test, feature_names, metrics, model, output_dir)
    
    print(f"\n🎉 流程完成! 最终测试准确率: {metrics['Accuracy']:.4f}")

if __name__ == "__main__":
    main()
