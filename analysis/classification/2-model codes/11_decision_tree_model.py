#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
决策树模型训练脚本
从预处理好的数据中加载训练集和测试集，进行模型训练、调参和评估
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pickle
import sys
from sklearn.tree import DecisionTreeClassifier
from skopt.space import Integer, Categorical
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
    """Draw confusion matrix: English labels, square, PDF support, small colorbar, single-color cmap."""
    fig, ax = plt.subplots(figsize=(6, 6))
    cmap_choice = 'Greens'
    sns.heatmap(
        cm,
        annot=True,
        fmt='d',
        cmap=cmap_choice,
        xticklabels=[f'Class {i}' for i in range(len(class_names))],
        yticklabels=[f'Class {i}' for i in range(len(class_names))],
        ax=ax,
        square=True,
        cbar=True,
        cbar_kws={'shrink': 0.6, 'pad': 0.02}
    )
    ax.set_xlabel('Predicted label')
    ax.set_ylabel('True label')
    ax.set_title(f'{model_name} (ANOVA 100 genes) Confusion Matrix', fontsize=12, pad=15)
    plt.tight_layout()

    if save_path:
        base, ext = os.path.splitext(save_path)
        if ext.lower() == '.pdf':
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Confusion matrix saved to: {save_path}")
        else:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            pdf_path = base + '.pdf'
            plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
            print(f"Confusion matrix saved to: {save_path} and {pdf_path}")

    plt.show()

def train_decision_tree(X_train, y_train, X_val, y_val):
    """训练决策树模型并进行贝叶斯优化"""
    print("\n🚀 开始训练决策树模型...")
    
    try:
        # 定义基础模型
        model = DecisionTreeClassifier(random_state=42)
        
        # 定义超参数搜索空间
        params = {
            'max_depth': Integer(3, 20),
            'min_samples_split': Integer(2, 20),
            'min_samples_leaf': Integer(1, 10),
            'criterion': Categorical(['gini', 'entropy']),
            'max_features': Categorical(['sqrt', 'log2', None]),
            'splitter': Categorical(['best', 'random'])
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
        cm_path = os.path.join(output_dir, 'confusion_matrix_decision_tree.pdf')
        plot_confusion_matrix(
            metrics['Confusion Matrix'], 
            class_names, 
            'Decision Tree (ANOVA 100 genes)',
            save_path=cm_path
        )
        
        return metrics, y_pred, y_proba
        
    except Exception as e:
        print(f"❌ 模型评估失败: {str(e)}")
        traceback.print_exc()
        return None, None, None

def analyze_feature_importance(model, feature_names, output_dir):
    """分析特征重要性"""
    print("\n🔍 分析特征重要性...")
    
    try:
        # 获取特征重要性
        importances = model.feature_importances_
        
        # 创建特征重要性DataFrame
        feature_importance_df = pd.DataFrame({
            'feature': feature_names,
            'importance': importances
        }).sort_values('importance', ascending=False)
        
        # 保存特征重要性
        importance_path = os.path.join(output_dir, 'decision_tree_feature_importance.csv')
        feature_importance_df.to_csv(importance_path, index=False, encoding='utf-8-sig')
        
        # 绘制前20个最重要特征
        plt.figure(figsize=(10, 8))
        top_features = feature_importance_df.head(20)
        sns.barplot(data=top_features, x='importance', y='feature', palette='viridis')
        plt.title('决策树模型 - 前20个最重要特征')
        plt.xlabel('特征重要性')
        plt.ylabel('基因名称')
        plt.tight_layout()
        
        # 保存特征重要性图
        importance_plot_path = os.path.join(output_dir, 'decision_tree_feature_importance.png')
        plt.savefig(importance_plot_path, dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"✅ 特征重要性分析完成!")
        print(f"  特征重要性文件: {importance_path}")
        print(f"  特征重要性图像: {importance_plot_path}")
        print(f"\n前10个最重要特征:")
        for i, (_, row) in enumerate(feature_importance_df.head(10).iterrows(), 1):
            print(f"  {i:2d}. {row['feature']}: {row['importance']:.4f}")
        
        return feature_importance_df
        
    except Exception as e:
        print(f"❌ 特征重要性分析失败: {str(e)}")
        traceback.print_exc()
        return None

def save_results(X_test, y_test, y_pred, names_test, feature_names, metrics, model, bayes_search, feature_importance_df, output_dir):
    """保存模型结果"""
    print("\n💾 正在保存结果...")
    
    try:
        # 保存测试集预测结果
        test_results_df = pd.DataFrame(X_test, columns=feature_names)
        test_results_df.insert(0, 'SampleID', names_test)
        test_results_df['True_label'] = y_test.values
        test_results_df['Predicted_label'] = y_pred
        
        results_path = os.path.join(output_dir, 'decision_tree_test_results.csv')
        test_results_df.to_csv(results_path, index=False, encoding='utf-8-sig')
        
        # 保存模型
        model_path = os.path.join(output_dir, 'decision_tree_model.pkl')
        with open(model_path, 'wb') as f:
            pickle.dump(model, f)
        
        # 保存评估指标
        metrics_path = os.path.join(output_dir, 'decision_tree_metrics.pkl')
        with open(metrics_path, 'wb') as f:
            pickle.dump(metrics, f)
        
        # 保存最佳参数
        params_path = os.path.join(output_dir, 'decision_tree_best_params.pkl')
        with open(params_path, 'wb') as f:
            pickle.dump(bayes_search.best_params_, f)
        
        # 保存详细报告
        report_path = os.path.join(output_dir, 'decision_tree_report.txt')
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("决策树模型评估报告\n")
            f.write("=" * 50 + "\n\n")
            
            f.write("最佳超参数:\n")
            f.write("-" * 20 + "\n")
            for param, value in bayes_search.best_params_.items():
                f.write(f"{param}: {value}\n")
            f.write(f"\n最佳交叉验证得分: {bayes_search.best_score_:.4f}\n\n")
            
            f.write("测试集性能指标:\n")
            f.write("-" * 20 + "\n")
            f.write(f"准确率: {metrics['Accuracy']:.4f}\n")
            f.write(f"多类别AUC: {metrics['AUC']:.4f}\n")
            f.write(f"宏平均F1: {metrics['Macro F1']:.4f}\n")
            f.write(f"微平均F1: {metrics['Micro F1']:.4f}\n\n")
            
            if feature_importance_df is not None:
                f.write("前20个最重要特征:\n")
                f.write("-" * 30 + "\n")
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
    print("--- 🚀 决策树模型训练与评估 ---")
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
    model, bayes_search = train_decision_tree(X_train, y_train, X_val, y_val)
    
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
    
    # 4. 分析特征重要性
    feature_importance_df = analyze_feature_importance(model, feature_names, figure_dir)
    
    # 5. 保存结果
    save_results(
        X_test, y_test, y_pred, names_test, 
        feature_names, metrics, model, bayes_search, feature_importance_df, output_dir
    )
    
    print(f"\n🎉 决策树模型训练和评估完成!")
    print(f"   最终测试准确率: {metrics['Accuracy']:.4f}")

if __name__ == "__main__":
    main()
