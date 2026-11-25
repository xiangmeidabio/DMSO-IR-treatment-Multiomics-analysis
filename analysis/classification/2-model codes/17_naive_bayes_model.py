#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Naive Bayes模型训练脚本
从预处理好的数据中加载训练集和测试集，进行朴素贝叶斯模型训练、调参和评估
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pickle
import sys
from sklearn.naive_bayes import GaussianNB
from log_helper import setup_logging

# initialize logging
setup_logging(__file__)
from skopt.space import Real
from skopt import BayesSearchCV
from sklearn.metrics import accuracy_score, roc_auc_score, f1_score, confusion_matrix, classification_report
import traceback

def load_preprocessed_data():
    """加载预处理好的数据"""
    # 获取脚本所在目录
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    results_dir = os.path.join(project_root, "results")
    
    print("正在加载预处理好的数据...")
    
    try:
        # 加载训练集
        with open(os.path.join(results_dir, 'train_data.pkl'), 'rb') as f:
            train_data = pickle.load(f)
        
        # 加载验证集
        with open(os.path.join(results_dir, 'val_data.pkl'), 'rb') as f:
            val_data = pickle.load(f)
            
        # 加载测试集
        with open(os.path.join(results_dir, 'test_data.pkl'), 'rb') as f:
            test_data = pickle.load(f)
            
        # 加载元数据
        with open(os.path.join(results_dir, 'metadata.pkl'), 'rb') as f:
            metadata = pickle.load(f)
        
        print("✅ 数据加载成功!")
        print(f"训练集形状: {train_data['X_train'].shape}")
        print(f"验证集形状: {val_data['X_val'].shape}")
        print(f"测试集形状: {test_data['X_test'].shape}")
        print(f"特征数量: {len(metadata['feature_names'])}")
        print(f"类别分布: {metadata['class_distribution']}")
        
        return train_data, val_data, test_data, metadata
        
    except FileNotFoundError as e:
        print(f"❌ 错误: 找不到预处理数据文件: {e}")
        print("请先运行 data_preprocessing.py 生成预处理数据")
        sys.exit(1)
    except Exception as e:
        print(f"❌ 加载数据时出错: {e}")
        sys.exit(1)

def calculate_metrics(y_true, y_pred, y_proba):
    """计算模型评估指标"""
    metrics = {}
    metrics['Accuracy'] = accuracy_score(y_true, y_pred)
    metrics['AUC'] = roc_auc_score(y_true, y_proba, multi_class='ovr')
    metrics['Macro F1'] = f1_score(y_true, y_pred, average='macro')
    metrics['Micro F1'] = f1_score(y_true, y_pred, average='micro')
    metrics['Confusion Matrix'] = confusion_matrix(y_true, y_pred)
    metrics['Classification Report'] = classification_report(
        y_true, y_pred, target_names=[f'class{int(i)}' for i in np.unique(y_true)], output_dict=True
    )
    return metrics

def plot_confusion_matrix(cm, class_names, model_name, save_path=None):
    """绘制混淆矩阵"""
    plt.figure(figsize=(6, 6))
    ax = sns.heatmap(
        cm,
        annot=True,
        fmt='d',
        cmap='Reds',
        xticklabels=class_names,
        yticklabels=class_names,
        cbar_kws={'shrink': 0.6}
    )
    try:
        cbar = ax.collections[0].colorbar
        cbar.set_label('')
        cbar.ax.tick_params(labelsize=8)
    except Exception:
        pass
    plt.xlabel('Predicted label')
    plt.ylabel('True label')
    plt.title(f'{model_name} Confusion Matrix')
    plt.tight_layout()

    if save_path:
        base, ext = os.path.splitext(save_path)
        png_path = base + '.png'
        pdf_path = base + '.pdf'
        plt.savefig(png_path, dpi=300, bbox_inches='tight')
        plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
        print(f"Confusion matrix saved to: {png_path} and {pdf_path}")

    plt.show()

def train_naive_bayes(X_train, y_train, X_val, y_val):
    """训练Naive Bayes模型并进行贝叶斯优化"""
    print("\n🚀 开始训练Naive Bayes模型...")
    
    try:
        # 合并训练集和验证集用于训练
        X_train_combined = np.vstack([X_train, X_val])
        y_train_combined = np.concatenate([y_train, y_val])
        
        print(f"合并后的训练集大小: {X_train_combined.shape}")
        
        # 定义基础模型
        model = GaussianNB()
        
        # 定义超参数搜索空间
        # Naive Bayes的参数相对较少，主要是方差平滑参数
        params = {
            'var_smoothing': Real(1e-12, 1e-5, prior='log-uniform'),
            'priors': [None]  # 使用默认的先验概率
        }
        
        # 贝叶斯优化
        print("正在进行贝叶斯超参数优化...")
        bayes_search = BayesSearchCV(
            model, params,
            n_iter=30,  # Naive Bayes参数较少，不需要太多搜索次数
            cv=5,       # 5折交叉验证
            scoring='accuracy',
            random_state=42,
            n_jobs=-1,
            verbose=1
        )
        
        bayes_search.fit(X_train_combined, y_train_combined)
        best_model = bayes_search.best_estimator_
        
        print(f"✅ 模型训练完成!")
        print(f"最佳参数: {bayes_search.best_params_}")
        print(f"最佳交叉验证得分: {bayes_search.best_score_:.4f}")
        
        return best_model, bayes_search
        
    except Exception as e:
        print(f"❌ 模型训练失败: {str(e)}")
        traceback.print_exc()
        return None, None

def evaluate_model(model, X_test, y_test, feature_names, class_names, output_dir):
    """评估模型性能"""
    print("\n📊 开始模型评估...")
    
    try:
        # 预测
        y_pred = model.predict(X_test)
        y_proba = model.predict_proba(X_test)
        
        # 计算指标
        metrics = calculate_metrics(y_test, y_pred, y_proba)
        
        # 打印结果
        print(f"\n🎯 测试集性能指标:")
        print(f"  准确率 (Accuracy): {metrics['Accuracy']:.4f}")
        print(f"  多类别AUC: {metrics['AUC']:.4f}")
        print(f"  宏平均F1: {metrics['Macro F1']:.4f}")
        print(f"  微平均F1: {metrics['Micro F1']:.4f}")
        
        # 打印模型信息
        print(f"\n🔍 模型信息:")
        print(f"  方差平滑参数: {model.var_smoothing:.2e}")
        print(f"  特征数量: {len(feature_names)}")
        print(f"  类别数量: {len(model.classes_)}")
        print(f"  类别: {model.classes_}")
        
        # 分析类别先验概率
        if hasattr(model, 'class_prior_') and model.class_prior_ is not None:
            print(f"\n📈 类别先验概率:")
            for i, (class_label, prior) in enumerate(zip(model.classes_, model.class_prior_)):
                print(f"  类别 {class_label}: {prior:.4f}")
        
        # 分析特征方差（选择前10个特征）
        if hasattr(model, 'var_'):
            print(f"\n📊 前10个特征的平均方差:")
            feature_vars = np.mean(model.var_, axis=0)  # 对所有类别求平均
            top_var_indices = np.argsort(feature_vars)[-10:][::-1]
            for i, idx in enumerate(top_var_indices):
                print(f"  {i+1}. {feature_names[idx]}: {feature_vars[idx]:.6f}")
        
        # 绘制混淆矩阵
        cm_path = os.path.join(output_dir, 'confusion_matrix_naive_bayes.png')
        plot_confusion_matrix(
            metrics['Confusion Matrix'], 
            class_names, 
            'Naive Bayes (ANOVA 100基因)',
            save_path=cm_path
        )
        
        return metrics, y_pred, y_proba
        
    except Exception as e:
        print(f"❌ 模型评估失败: {str(e)}")
        traceback.print_exc()
        return None, None, None

def plot_class_probabilities(model, save_path=None):
    """绘制类别先验概率分布"""
    if hasattr(model, 'class_prior_') and model.class_prior_ is not None:
        classes = [f'class{int(i)}' for i in model.classes_]
        probabilities = model.class_prior_

        plt.figure(figsize=(6, 6))
        plt.bar(classes, probabilities, color='lightcoral', edgecolor='darkred')
        plt.xlabel('Class')
        plt.ylabel('Prior probability')
        plt.title('Naive Bayes - Class prior probabilities')
        plt.ylim(0, max(probabilities) * 1.1)

        # add value labels
        for i, prob in enumerate(probabilities):
            plt.text(i, prob + max(probabilities) * 0.01, f'{prob:.3f}',
                     ha='center', va='bottom', fontweight='bold')

        plt.tight_layout()

        if save_path:
            base, ext = os.path.splitext(save_path)
            plt.savefig(base + '.png', dpi=300, bbox_inches='tight')
            plt.savefig(base + '.pdf', dpi=300, bbox_inches='tight')
            print(f"Class prior plot saved to: {base}.png and {base}.pdf")

        plt.show()
    else:
        print("Model uses default priors, skipping class prior plot")

def plot_feature_variance_analysis(model, feature_names, save_path=None, top_n=20):
    """Plot feature variance analysis"""
    if hasattr(model, 'var_'):
        plt.figure(figsize=(6, 6))

        # compute mean variance per feature across classes
        feature_vars = np.mean(model.var_, axis=0)

        # select top N features by variance
        top_var_indices = np.argsort(feature_vars)[-top_n:][::-1]
        top_features = [feature_names[i] for i in top_var_indices]
        top_variances = feature_vars[top_var_indices]

        # horizontal bar plot
        plt.barh(range(len(top_features)), top_variances, color='lightcoral')
        plt.yticks(range(len(top_features)), top_features)
        plt.xlabel('Average variance')
        plt.title(f'Naive Bayes - Top {top_n} high variance features')
        plt.gca().invert_yaxis()
        plt.tight_layout()

        if save_path:
            base, ext = os.path.splitext(save_path)
            plt.savefig(base + '.png', dpi=300, bbox_inches='tight')
            plt.savefig(base + '.pdf', dpi=300, bbox_inches='tight')
            print(f"Feature variance plot saved to: {base}.png and {base}.pdf")

        plt.show()
    else:
        print("No feature variance information available")

def save_results(X_test, y_test, y_pred, names_test, feature_names, metrics, 
                model, output_dir, figure_dir):
    """保存模型结果"""
    print("\n💾 正在保存结果...")
    
    try:
        # 保存测试集预测结果
        test_results_df = pd.DataFrame(X_test, columns=feature_names)
        test_results_df.insert(0, 'SampleID', names_test)
        test_results_df['True_label'] = y_test.values if hasattr(y_test, 'values') else y_test
        test_results_df['Predicted_label'] = y_pred
        
        results_path = os.path.join(output_dir, 'naive_bayes_test_results.csv')
        test_results_df.to_csv(results_path, index=False, encoding='utf-8-sig')
        
        # 保存模型
        model_path = os.path.join(output_dir, 'naive_bayes_model.pkl')
        with open(model_path, 'wb') as f:
            pickle.dump(model, f)
        
        # 保存评估指标
        metrics_path = os.path.join(output_dir, 'naive_bayes_metrics.pkl')
        with open(metrics_path, 'wb') as f:
            pickle.dump(metrics, f)
        
        # 绘制并保存类别概率分布图
        prob_plot_path = os.path.join(figure_dir, 'class_probabilities_naive_bayes.png')
        plot_class_probabilities(model, save_path=prob_plot_path)
        
        # 绘制并保存特征方差分析图
        variance_plot_path = os.path.join(figure_dir, 'feature_variance_naive_bayes.png')
        plot_feature_variance_analysis(model, feature_names, save_path=variance_plot_path)
        
        # 保存详细报告
        report_path = os.path.join(output_dir, 'naive_bayes_report.txt')
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("Naive Bayes（朴素贝叶斯）模型评估报告\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"准确率: {metrics['Accuracy']:.4f}\n")
            f.write(f"多类别AUC: {metrics['AUC']:.4f}\n")
            f.write(f"宏平均F1: {metrics['Macro F1']:.4f}\n")
            f.write(f"微平均F1: {metrics['Micro F1']:.4f}\n\n")
            
            f.write("模型参数:\n")
            f.write("-" * 20 + "\n")
            f.write(f"方差平滑参数: {model.var_smoothing:.2e}\n")
            f.write(f"特征数量: {len(feature_names)}\n")
            f.write(f"类别数量: {len(model.classes_)}\n")
            f.write(f"类别: {list(model.classes_)}\n\n")
            
            f.write("Naive Bayes特点:\n")
            f.write("-" * 20 + "\n")
            f.write("1. 基于贝叶斯定理和特征独立性假设\n")
            f.write("2. 训练速度快，内存需求低\n")
            f.write("3. 对小数据集表现良好\n")
            f.write("4. 对噪声和缺失数据相对鲁棒\n")
            f.write("5. 提供概率预测输出\n")
            f.write("6. 假设特征之间相互独立\n\n")
            
            if hasattr(model, 'class_prior_') and model.class_prior_ is not None:
                f.write("类别先验概率:\n")
                f.write("-" * 20 + "\n")
                for class_label, prior in zip(model.classes_, model.class_prior_):
                    f.write(f"类别 {class_label}: {prior:.4f}\n")
                f.write("\n")
            
            if hasattr(model, 'var_'):
                f.write("前20个高方差特征:\n")
                f.write("-" * 30 + "\n")
                feature_vars = np.mean(model.var_, axis=0)
                top_var_indices = np.argsort(feature_vars)[-20:][::-1]
                for i, idx in enumerate(top_var_indices):
                    f.write(f"{i+1:2d}. {feature_names[idx]:15s}: {feature_vars[idx]:.6f}\n")
                f.write("\n")
            
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
        print(f"  类别概率分布图: {prob_plot_path}")
        print(f"  特征方差分析图: {variance_plot_path}")
        
    except Exception as e:
        print(f"❌ 保存结果失败: {str(e)}")
        traceback.print_exc()

def main():
    """主函数"""
    print("\n" + "="*60)
    print("--- 🚀 Naive Bayes（朴素贝叶斯）模型训练与评估 ---")
    print("         (使用100个ANOVA基因)")
    print("="*60)
    
    # 设置输出目录
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    output_dir = os.path.join(project_root, "results")
    figure_dir = os.path.join(project_root, "figure")
    os.makedirs(figure_dir, exist_ok=True)
    
    # 1. 加载预处理数据
    train_data, val_data, test_data, metadata = load_preprocessed_data()
    
    X_train = train_data['X_train']
    y_train = train_data['y_train']
    X_val = val_data['X_val']
    y_val = val_data['y_val']
    X_test = test_data['X_test']
    y_test = test_data['y_test']
    names_test = test_data['names_test']
    feature_names = metadata['feature_names']
    
    # 准备类别名称
    unique_labels = np.unique(y_train)
    class_names = [f'Class {i}' for i in unique_labels]
    
    # 2. 训练模型
    model, bayes_search = train_naive_bayes(X_train, y_train, X_val, y_val)
    
    if model is None:
        print("❌ 模型训练失败，程序退出")
        return
    
    # 3. 评估模型
    metrics, y_pred, y_proba = evaluate_model(
        model, X_test, y_test, feature_names, class_names, figure_dir
    )
    
    if metrics is None:
        print("❌ 模型评估失败，程序退出")
        return
    
    # 4. 保存结果
    save_results(
        X_test, y_test, y_pred, names_test, 
        feature_names, metrics, model, output_dir, figure_dir
    )
    
    print(f"\n🎉 Naive Bayes模型训练和评估完成!")
    print(f"   最终测试准确率: {metrics['Accuracy']:.4f}")
    print(f"   方差平滑参数: {model.var_smoothing:.2e}")

if __name__ == "__main__":
    main()
