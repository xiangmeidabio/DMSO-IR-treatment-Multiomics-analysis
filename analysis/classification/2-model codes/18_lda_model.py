#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LDA（线性判别分析）模型训练脚本
从预处理好的数据中加载训练集和测试集，进行线性判别分析模型训练、调参和评估
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pickle
import sys
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
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
        cmap='cividis',
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

def train_lda(X_train, y_train, X_val, y_val):
    """训练LDA模型并进行贝叶斯优化"""
    print("\n🚀 开始训练LDA模型...")
    
    try:
        # 合并训练集和验证集用于训练
        X_train_combined = np.vstack([X_train, X_val])
        y_train_combined = np.concatenate([y_train, y_val])
        
        print(f"合并后的训练集大小: {X_train_combined.shape}")
        
        # 定义基础模型
        model = LinearDiscriminantAnalysis()
        
        # 定义超参数搜索空间
        params = {
            'solver': Categorical(['lsqr', 'eigen']),  # 添加eigen求解器
            'shrinkage': Real(0.0, 1.0),  # 收缩参数
            'tol': Real(1e-6, 1e-3, prior='log-uniform')  # 容忍度
        }
        
        # 贝叶斯优化
        print("正在进行贝叶斯超参数优化...")
        bayes_search = BayesSearchCV(
            model, params,
            n_iter=40,  # 增加搜索次数
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
        print(f"  求解器: {model.solver}")
        print(f"  收缩参数: {model.shrinkage}")
        print(f"  容忍度: {model.tol}")
        print(f"  类别数量: {len(model.classes_)}")
        print(f"  特征数量: {len(feature_names)}")
        
        # 分析判别函数和降维
        if hasattr(model, 'scalings_'):
            n_components = model.scalings_.shape[1]
            print(f"  判别成分数量: {n_components}")
            print(f"  降维后维度: {n_components}")
            
            # 计算解释方差比例
            if hasattr(model, 'explained_variance_ratio_'):
                print(f"  前3个成分解释方差比例: {model.explained_variance_ratio_[:3]}")
                print(f"  累计解释方差比例: {np.cumsum(model.explained_variance_ratio_)[:3]}")
        
        # 分析类别中心
        if hasattr(model, 'means_'):
            print(f"\n📈 类别中心信息:")
            for i, class_label in enumerate(model.classes_):
                mean_values = model.means_[i]
                print(f"  类别 {class_label}: 平均值范围 [{mean_values.min():.4f}, {mean_values.max():.4f}]")
        
        # 绘制混淆矩阵
        cm_path = os.path.join(output_dir, 'confusion_matrix_lda.png')
        plot_confusion_matrix(
            metrics['Confusion Matrix'], 
            class_names, 
            'LDA (ANOVA 100基因)',
            save_path=cm_path
        )
        
        return metrics, y_pred, y_proba
        
    except Exception as e:
        print(f"❌ 模型评估失败: {str(e)}")
        traceback.print_exc()
        return None, None, None

def plot_lda_components(model, X_test, y_test, save_path=None):
    """绘制LDA降维后的数据分布"""
    if hasattr(model, 'transform'):
        try:
            # 将测试数据投影到LDA空间
            X_lda = model.transform(X_test)
            
            # 如果有多个成分，绘制前两个成分
            if X_lda.shape[1] >= 2:
                plt.figure(figsize=(10, 8))
                
                # 为每个类别使用不同颜色
                unique_classes = np.unique(y_test)
                colors = plt.cm.Set1(np.linspace(0, 1, len(unique_classes)))
                
                for i, class_label in enumerate(unique_classes):
                    mask = y_test == class_label
                    plt.scatter(X_lda[mask, 0], X_lda[mask, 1], 
                              c=[colors[i]], label=f'Class {class_label}', 
                              alpha=0.7, s=50)
                
                plt.xlabel('第一判别成分')
                plt.ylabel('第二判别成分')
                plt.title('LDA降维后的数据分布')
                plt.legend()
                plt.grid(True, alpha=0.3)
                plt.tight_layout()
                
                if save_path:
                    base, ext = os.path.splitext(save_path)
                    plt.savefig(base + '.png', dpi=300, bbox_inches='tight')
                    plt.savefig(base + '.pdf', dpi=300, bbox_inches='tight')
                    print(f"LDA components plot saved to: {base}.png and {base}.pdf")
                
                plt.show()
            
            elif X_lda.shape[1] == 1:
                # 如果只有一个成分，绘制直方图
                plt.figure(figsize=(10, 6))
                
                unique_classes = np.unique(y_test)
                for class_label in unique_classes:
                    mask = y_test == class_label
                    plt.hist(X_lda[mask, 0], alpha=0.7, label=f'Class {class_label}', bins=20)
                
                plt.xlabel('判别成分值')
                plt.ylabel('频数')
                plt.title('LDA降维后的数据分布（单成分）')
                plt.legend()
                plt.grid(True, alpha=0.3)
                plt.tight_layout()
                
                if save_path:
                    base, ext = os.path.splitext(save_path)
                    plt.savefig(base + '.png', dpi=300, bbox_inches='tight')
                    plt.savefig(base + '.pdf', dpi=300, bbox_inches='tight')
                    print(f"LDA components plot saved to: {base}.png and {base}.pdf")
                
                plt.show()
                
        except Exception as e:
            print(f"绘制LDA成分图时出错: {e}")

def plot_explained_variance(model, save_path=None):
    """绘制解释方差比例图"""
    if hasattr(model, 'explained_variance_ratio_'):
        plt.figure(figsize=(10, 6))
        
        n_components = len(model.explained_variance_ratio_)
        components = range(1, n_components + 1)
        
        # 绘制解释方差比例
        plt.subplot(1, 2, 1)
        plt.bar(components, model.explained_variance_ratio_, color='skyblue', edgecolor='navy')
        plt.xlabel('判别成分')
        plt.ylabel('解释方差比例')
        plt.title('各成分解释方差比例')
        plt.xticks(components)
        
        # 绘制累计解释方差比例
        plt.subplot(1, 2, 2)
        cumsum_ratio = np.cumsum(model.explained_variance_ratio_)
        plt.plot(components, cumsum_ratio, 'o-', color='red', linewidth=2, markersize=8)
        plt.xlabel('判别成分')
        plt.ylabel('累计解释方差比例')
        plt.title('累计解释方差比例')
        plt.xticks(components)
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            base, ext = os.path.splitext(save_path)
            plt.savefig(base + '.png', dpi=300, bbox_inches='tight')
            plt.savefig(base + '.pdf', dpi=300, bbox_inches='tight')
            print(f"Explained variance plot saved to: {base}.png and {base}.pdf")
        
        plt.show()
    else:
        print("模型没有解释方差信息")

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
        
        results_path = os.path.join(output_dir, 'lda_test_results.csv')
        test_results_df.to_csv(results_path, index=False, encoding='utf-8-sig')
        
        # 保存模型
        model_path = os.path.join(output_dir, 'lda_model.pkl')
        with open(model_path, 'wb') as f:
            pickle.dump(model, f)
        
        # 保存评估指标
        metrics_path = os.path.join(output_dir, 'lda_metrics.pkl')
        with open(metrics_path, 'wb') as f:
            pickle.dump(metrics, f)
        
        # 绘制并保存LDA成分分布图
        components_plot_path = os.path.join(figure_dir, 'lda_components.png')
        plot_lda_components(model, X_test, y_test, save_path=components_plot_path)
        
        # 绘制并保存解释方差图
        variance_plot_path = os.path.join(figure_dir, 'lda_explained_variance.png')
        plot_explained_variance(model, save_path=variance_plot_path)
        
        # 保存详细报告
        report_path = os.path.join(output_dir, 'lda_report.txt')
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("LDA（线性判别分析）模型评估报告\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"准确率: {metrics['Accuracy']:.4f}\n")
            f.write(f"多类别AUC: {metrics['AUC']:.4f}\n")
            f.write(f"宏平均F1: {metrics['Macro F1']:.4f}\n")
            f.write(f"微平均F1: {metrics['Micro F1']:.4f}\n\n")
            
            f.write("模型参数:\n")
            f.write("-" * 20 + "\n")
            f.write(f"求解器: {model.solver}\n")
            f.write(f"收缩参数: {model.shrinkage}\n")
            f.write(f"容忍度: {model.tol}\n")
            f.write(f"类别数量: {len(model.classes_)}\n")
            f.write(f"特征数量: {len(feature_names)}\n")
            
            if hasattr(model, 'scalings_'):
                f.write(f"判别成分数量: {model.scalings_.shape[1]}\n")
            f.write("\n")
            
            f.write("LDA特点:\n")
            f.write("-" * 20 + "\n")
            f.write("1. 基于贝叶斯决策理论的线性分类器\n")
            f.write("2. 寻找最佳的线性判别函数\n")
            f.write("3. 同时进行分类和降维\n")
            f.write("4. 假设各类别协方差矩阵相同\n")
            f.write("5. 最大化类间距离，最小化类内距离\n")
            f.write("6. 提供概率预测输出\n\n")
            
            if hasattr(model, 'explained_variance_ratio_'):
                f.write("解释方差比例:\n")
                f.write("-" * 20 + "\n")
                for i, ratio in enumerate(model.explained_variance_ratio_):
                    f.write(f"成分 {i+1}: {ratio:.4f}\n")
                f.write(f"累计解释方差: {np.sum(model.explained_variance_ratio_):.4f}\n\n")
            
            if hasattr(model, 'means_'):
                f.write("类别中心信息:\n")
                f.write("-" * 20 + "\n")
                for i, class_label in enumerate(model.classes_):
                    mean_values = model.means_[i]
                    f.write(f"类别 {class_label}: 平均值范围 [{mean_values.min():.4f}, {mean_values.max():.4f}]\n")
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
        print(f"  LDA成分分布图: {components_plot_path}")
        print(f"  解释方差图: {variance_plot_path}")
        
    except Exception as e:
        print(f"❌ 保存结果失败: {str(e)}")
        traceback.print_exc()

def main():
    """主函数"""
    print("\n" + "="*60)
    print("--- 🚀 LDA（线性判别分析）模型训练与评估 ---")
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
    class_names = [f'class{int(i)}' for i in unique_labels]
    
    # 2. 训练模型
    model, bayes_search = train_lda(X_train, y_train, X_val, y_val)
    
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
    
    print(f"\n🎉 LDA模型训练和评估完成!")
    print(f"   最终测试准确率: {metrics['Accuracy']:.4f}")
    print(f"   求解器: {model.solver}")
    if hasattr(model, 'scalings_'):
        print(f"   判别成分数量: {model.scalings_.shape[1]}")

if __name__ == "__main__":
    main()
