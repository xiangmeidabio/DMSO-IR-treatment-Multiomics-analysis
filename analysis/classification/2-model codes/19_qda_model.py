#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QDA (二次判别分析) 模型训练脚本
从预处理好的数据中加载训练集和测试集，进行模型训练、调参和评估
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pickle
import sys
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from log_helper import setup_logging

# initialize logging
setup_logging(__file__)
from skopt.space import Real, Categorical
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
        y_true, y_pred, target_names=[f'Class {i}' for i in np.unique(y_true)], output_dict=True
    )
    return metrics

def plot_confusion_matrix(cm, class_names, model_name, save_path=None):
    """绘制混淆矩阵"""
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Purples', 
                xticklabels=class_names, yticklabels=class_names,
                cbar_kws={'shrink': 0.6})
    plt.xlabel('Predicted label')
    plt.ylabel('True label')
    plt.title(f'{model_name} Confusion Matrix')
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    base, ext = os.path.splitext(save_path)
    png_path = base + '.png'
    pdf_path = base + '.pdf'
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
    print(f"Confusion matrix saved to: {png_path} and {pdf_path}")
    
    plt.show()

def train_qda(X_train, y_train, X_val, y_val):
    """训练QDA模型并进行贝叶斯优化"""
    print("\n🚀 开始训练QDA模型...")
    
    try:
        # 定义基础模型
        model = QuadraticDiscriminantAnalysis()
        
        # 定义超参数搜索空间
        params = {
            'reg_param': Real(0.0, 1.0),
            'store_covariance': Categorical([True, False]),
            'tol': Real(1e-6, 1e-3, prior='log-uniform')
        }
        
        # 贝叶斯优化
        print("正在进行贝叶斯超参数优化...")
        bayes_search = BayesSearchCV(
            model, params,
            n_iter=30,  # 搜索次数
            cv=3,       # 交叉验证折数
            scoring='accuracy',
            random_state=42,
            n_jobs=-1,  # 使用所有可用CPU核心
            verbose=1
        )
        
        bayes_search.fit(X_train, y_train)
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
        
        # 绘制混淆矩阵
        cm_path = os.path.join(output_dir, 'confusion_matrix_qda.png')
        plot_confusion_matrix(
            metrics['Confusion Matrix'], 
            class_names, 
            'QDA (ANOVA 100基因)',
            save_path=cm_path
        )
        
        return metrics, y_pred, y_proba
        
    except Exception as e:
        print(f"❌ 模型评估失败: {str(e)}")
        traceback.print_exc()
        return None, None, None

def analyze_class_statistics(model, X_train, y_train, feature_names, output_dir):
    """分析QDA模型的类别统计信息"""
    print("\n🔍 分析QDA类别统计信息...")
    
    try:
        # 获取类别先验概率
        class_priors = model.priors_
        classes = model.classes_
        
        # 获取类别均值
        class_means = model.means_
        
        # 创建类别统计DataFrame
        class_stats = []
        for i, class_label in enumerate(classes):
            stats = {
                'Class': f'Class {class_label}',
                'Prior_Probability': class_priors[i],
                'Sample_Count': np.sum(y_train == class_label),
                'Mean_Feature_Value': np.mean(class_means[i])
            }
            class_stats.append(stats)
        
        class_stats_df = pd.DataFrame(class_stats)
        
        # 保存类别统计信息
        stats_path = os.path.join(output_dir, 'qda_class_statistics.csv')
        class_stats_df.to_csv(stats_path, index=False, encoding='utf-8-sig')
        
        # 绘制类别先验概率
        plt.figure(figsize=(10, 6))
        plt.subplot(1, 2, 1)
        plt.bar(class_stats_df['Class'], class_stats_df['Prior_Probability'], color='purple', alpha=0.7)
        plt.title('各类别先验概率')
        plt.xlabel('类别')
        plt.ylabel('先验概率')
        plt.xticks(rotation=45)
        
        # 绘制类别样本数量
        plt.subplot(1, 2, 2)
        plt.bar(class_stats_df['Class'], class_stats_df['Sample_Count'], color='indigo', alpha=0.7)
        plt.title('各类别样本数量')
        plt.xlabel('类别')
        plt.ylabel('样本数量')
        plt.xticks(rotation=45)
        
        plt.tight_layout()
        
        # 保存类别统计图
        stats_plot_path = os.path.join(output_dir, 'qda_class_statistics.png')
        plt.savefig(stats_plot_path, dpi=300, bbox_inches='tight')
        plt.show()
        
        # 分析特征重要性（基于类别均值差异）
        feature_importance = []
        for j, feature in enumerate(feature_names):
            # 计算各类别在该特征上的均值差异
            feature_means = [class_means[i][j] for i in range(len(classes))]
            importance = np.std(feature_means)  # 使用标准差作为重要性指标
            feature_importance.append({
                'feature': feature,
                'importance': importance
            })
        
        feature_importance_df = pd.DataFrame(feature_importance).sort_values('importance', ascending=False)
        
        # 保存特征重要性
        importance_path = os.path.join(output_dir, 'qda_feature_importance.csv')
        feature_importance_df.to_csv(importance_path, index=False, encoding='utf-8-sig')
        
        # 绘制前20个最重要特征
        plt.figure(figsize=(10, 8))
        top_features = feature_importance_df.head(20)
        sns.barplot(data=top_features, x='importance', y='feature', palette='viridis')
        plt.title('QDA模型 - 前20个最重要特征\n(基于类别均值差异)')
        plt.xlabel('重要性 (类别均值标准差)')
        plt.ylabel('基因名称')
        plt.tight_layout()
        
        # 保存特征重要性图
        importance_plot_path = os.path.join(output_dir, 'qda_feature_importance.png')
        plt.savefig(importance_plot_path, dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"✅ 类别统计分析完成!")
        print(f"  类别统计文件: {stats_path}")
        print(f"  类别统计图像: {stats_plot_path}")
        print(f"  特征重要性文件: {importance_path}")
        print(f"  特征重要性图像: {importance_plot_path}")
        
        print(f"\n类别统计摘要:")
        for _, row in class_stats_df.iterrows():
            print(f"  {row['Class']}: 先验概率={row['Prior_Probability']:.4f}, 样本数={row['Sample_Count']}")
        
        print(f"\n前10个最重要特征:")
        for i, (_, row) in enumerate(feature_importance_df.head(10).iterrows(), 1):
            print(f"  {i:2d}. {row['feature']}: {row['importance']:.4f}")
        
        return class_stats_df, feature_importance_df
        
    except Exception as e:
        print(f"❌ 类别统计分析失败: {str(e)}")
        traceback.print_exc()
        return None, None

def save_results(X_test, y_test, y_pred, names_test, feature_names, metrics, model, bayes_search, 
                class_stats_df, feature_importance_df, output_dir):
    """保存模型结果"""
    print("\n💾 正在保存结果...")
    
    try:
        # 保存测试集预测结果
        test_results_df = pd.DataFrame(X_test, columns=feature_names)
        test_results_df.insert(0, 'SampleID', names_test)
        test_results_df['True_label'] = y_test.values
        test_results_df['Predicted_label'] = y_pred
        
        results_path = os.path.join(output_dir, 'qda_test_results.csv')
        test_results_df.to_csv(results_path, index=False, encoding='utf-8-sig')
        
        # 保存模型
        model_path = os.path.join(output_dir, 'qda_model.pkl')
        with open(model_path, 'wb') as f:
            pickle.dump(model, f)
        
        # 保存评估指标
        metrics_path = os.path.join(output_dir, 'qda_metrics.pkl')
        with open(metrics_path, 'wb') as f:
            pickle.dump(metrics, f)
        
        # 保存最佳参数
        params_path = os.path.join(output_dir, 'qda_best_params.pkl')
        with open(params_path, 'wb') as f:
            pickle.dump(bayes_search.best_params_, f)
        
        # 保存详细报告
        report_path = os.path.join(output_dir, 'qda_report.txt')
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("QDA (二次判别分析) 模型评估报告\n")
            f.write("=" * 50 + "\n\n")
            
            f.write("最佳超参数:\n")
            f.write("-" * 20 + "\n")
            for param, value in bayes_search.best_params_.items():
                f.write(f"{param}: {value}\n")
            f.write(f"\n最佳交叉验证得分: {bayes_search.best_score_:.4f}\n\n")
            
            if class_stats_df is not None:
                f.write("类别统计信息:\n")
                f.write("-" * 20 + "\n")
                for _, row in class_stats_df.iterrows():
                    f.write(f"{row['Class']}:\n")
                    f.write(f"  先验概率: {row['Prior_Probability']:.4f}\n")
                    f.write(f"  样本数量: {row['Sample_Count']}\n")
                    f.write(f"  平均特征值: {row['Mean_Feature_Value']:.4f}\n\n")
            
            f.write("测试集性能指标:\n")
            f.write("-" * 20 + "\n")
            f.write(f"准确率: {metrics['Accuracy']:.4f}\n")
            f.write(f"多类别AUC: {metrics['AUC']:.4f}\n")
            f.write(f"宏平均F1: {metrics['Macro F1']:.4f}\n")
            f.write(f"微平均F1: {metrics['Micro F1']:.4f}\n\n")
            
            if feature_importance_df is not None:
                f.write("前20个最重要特征 (基于类别均值差异):\n")
                f.write("-" * 40 + "\n")
                for i, (_, row) in enumerate(feature_importance_df.head(20).iterrows(), 1):
                    f.write(f"{i:2d}. {row['feature']}: {row['importance']:.4f}\n")
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
        print(f"  最佳参数: {params_path}")
        print(f"  详细报告: {report_path}")
        
    except Exception as e:
        print(f"❌ 保存结果失败: {str(e)}")
        traceback.print_exc()

def main():
    """主函数"""
    print("\n" + "="*60)
    print("--- 🚀 QDA (二次判别分析) 模型训练与评估 ---")
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
    model, bayes_search = train_qda(X_train, y_train, X_val, y_val)
    
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
    
    # 4. 分析类别统计信息
    class_stats_df, feature_importance_df = analyze_class_statistics(
        model, X_train, y_train, feature_names, figure_dir
    )
    
    # 5. 保存结果
    save_results(
        X_test, y_test, y_pred, names_test, 
        feature_names, metrics, model, bayes_search, 
        class_stats_df, feature_importance_df, output_dir
    )
    
    print(f"\n🎉 QDA模型训练和评估完成!")
    print(f"   最终测试准确率: {metrics['Accuracy']:.4f}")

if __name__ == "__main__":
    main()
