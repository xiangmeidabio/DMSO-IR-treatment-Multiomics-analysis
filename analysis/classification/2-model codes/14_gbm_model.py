#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GBM（梯度提升机）模型训练脚本
从预处理好的数据中加载训练集和测试集，进行GBM模型训练、调参和评估
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pickle
import sys
from sklearn.ensemble import GradientBoostingClassifier
from skopt import BayesSearchCV
from skopt.space import Integer, Real, Categorical
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
        cmap='Oranges',
        xticklabels=class_names,
        yticklabels=class_names,
        cbar_kws={'shrink': 0.6}
    )
    # remove colorbar label if present and reduce tick size
    try:
        cbar = ax.collections[0].colorbar
        cbar.set_label('')
        cbar.ax.tick_params(labelsize=8)
    except Exception:
        pass
    plt.xlabel('Predicted label', fontsize=11)
    plt.ylabel('True label', fontsize=11)
    plt.title(f'{model_name} Confusion Matrix', fontsize=12, pad=15)
    plt.tight_layout()
    
    if save_path:
        base, ext = os.path.splitext(save_path)
        png_path = base + '.png'
        pdf_path = base + '.pdf'
        plt.savefig(png_path, dpi=300, bbox_inches='tight')
        plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
        print(f"Confusion matrix saved to: {png_path} and {pdf_path}")
    
    plt.show()

def train_gbm(X_train, y_train, X_val, y_val):
    """训练GBM模型并进行贝叶斯优化"""
    print("\n🚀 开始训练GBM模型...")
    print("=" * 50)
    
    try:
        # 合并训练集和验证集用于训练
        X_train_combined = np.vstack([X_train, X_val])
        y_train_combined = np.concatenate([y_train, y_val])
        
        print(f"合并后的训练集大小: {X_train_combined.shape}")
        
        # 定义基础模型
        model = GradientBoostingClassifier(
            random_state=42,
            loss='log_loss'  # 多分类损失函数
        )
        
        # 定义超参数搜索空间
        param_space = {
            'n_estimators': Integer(100, 300),
            'max_depth': Integer(3, 12),
            'learning_rate': Real(0.01, 0.3, prior='log-uniform'),
            'min_samples_split': Integer(2, 20),
            'min_samples_leaf': Integer(1, 10),
            'subsample': Real(0.7, 1.0),
            'max_features': Categorical(['sqrt', 'log2', None]),
            'min_impurity_decrease': Real(0.0, 0.1),
            'validation_fraction': Real(0.1, 0.2)  # 早停验证比例
        }
        
        # 贝叶斯优化
        print(f"🔍 开始贝叶斯超参数搜索...")
        bayes_search = BayesSearchCV(
            estimator=model,
            search_spaces=param_space,
            n_iter=50,  # 增加搜索次数
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
        
        # 特征重要性分析
        feature_importance = pd.DataFrame({
            'feature': feature_names,
            'importance': model.feature_importances_
        }).sort_values('importance', ascending=False)
        
        print(f"\n🔍 前10个最重要的特征:")
        for i, (_, row) in enumerate(feature_importance.head(10).iterrows()):
            print(f"  {i+1}. {row['feature']}: {row['importance']:.4f}")
        
        # 保存特征重要性
        importance_path = os.path.join(output_dir, 'gbm_feature_importance.csv')
        feature_importance.to_csv(importance_path, index=False, encoding='utf-8-sig')
        print(f"\n特征重要性已保存到: {importance_path}")
        
        # 打印模型信息
        print(f"\n🌳 模型信息:")
        print(f"  提升树数量: {model.n_estimators}")
        print(f"  最大深度: {model.max_depth}")
        print(f"  学习率: {model.learning_rate:.4f}")
        print(f"  子采样比例: {model.subsample:.4f}")
        print(f"  最小分割样本数: {model.min_samples_split}")
        print(f"  最小叶子样本数: {model.min_samples_leaf}")
        print(f"  特征选择策略: {model.max_features}")
        
        # 绘制混淆矩阵
        cm_path = os.path.join(output_dir, 'confusion_matrix_gbm.png')
        plot_confusion_matrix(
            metrics['Confusion Matrix'], 
            class_names, 
            'GBM (ANOVA 100 genes)',
            save_path=cm_path
        )
        
        return metrics, y_pred, y_proba, feature_importance
        
    except Exception as e:
        print(f"❌ 模型评估失败: {str(e)}")
        traceback.print_exc()
        return None, None, None, None

def plot_feature_importance(feature_importance, save_path=None, top_n=20):
    """绘制特征重要性图"""
    plt.figure(figsize=(6, 6))
    # 选择前N个最重要的特征
    top_features = feature_importance.head(top_n)

    # 水平条形图，单色（Oranges）
    plt.barh(range(len(top_features)), top_features['importance'], color='orange')
    plt.yticks(range(len(top_features)), top_features['feature'])
    plt.xlabel('Feature importance')
    plt.title(f'GBM - Top {top_n} important features')
    plt.gca().invert_yaxis()
    plt.tight_layout()

    if save_path:
        base, ext = os.path.splitext(save_path)
        png_path = base + '.png'
        pdf_path = base + '.pdf'
        plt.savefig(png_path, dpi=300, bbox_inches='tight')
        plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
        print(f"Feature importance plot saved to: {png_path} and {pdf_path}")

    plt.show()

def plot_training_progress(model, save_path=None):
    """绘制训练过程中的损失变化"""
    if hasattr(model, 'train_score_'):
        plt.figure(figsize=(10, 6))
        
        # 绘制训练损失
        plt.plot(model.train_score_, label='训练损失', color='blue', linewidth=2)
        
        # 如果有验证损失，也绘制出来
        if hasattr(model, 'validation_scores_') and model.validation_scores_ is not None:
            plt.plot(model.validation_scores_, label='验证损失', color='red', linewidth=2)
        
        plt.xlabel('提升轮次')
        plt.ylabel('损失值')
        plt.title('GBM训练过程 - 损失变化')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"训练过程图已保存到: {save_path}")
        
        plt.show()
    else:
        print("模型没有训练损失记录，跳过训练过程图绘制")

def save_results(X_test, y_test, y_pred, names_test, feature_names, metrics, 
                model, feature_importance, output_dir, figure_dir):
    """保存模型结果"""
    print("\n💾 正在保存结果...")
    
    try:
        # 保存测试集预测结果
        test_results_df = pd.DataFrame(X_test, columns=feature_names)
        test_results_df.insert(0, 'SampleID', names_test)
        test_results_df['True_label'] = y_test.values if hasattr(y_test, 'values') else y_test
        test_results_df['Predicted_label'] = y_pred
        
        results_path = os.path.join(output_dir, 'gbm_test_results.csv')
        test_results_df.to_csv(results_path, index=False, encoding='utf-8-sig')
        
        # 保存模型
        model_path = os.path.join(output_dir, 'gbm_model.pkl')
        with open(model_path, 'wb') as f:
            pickle.dump(model, f)
        
        # 保存评估指标
        metrics_path = os.path.join(output_dir, 'gbm_metrics.pkl')
        with open(metrics_path, 'wb') as f:
            pickle.dump(metrics, f)
        
        # 绘制并保存特征重要性图
        importance_plot_path = os.path.join(figure_dir, 'feature_importance_gbm.png')
        plot_feature_importance(feature_importance, save_path=importance_plot_path)
        
        # 绘制并保存训练过程图
        training_plot_path = os.path.join(figure_dir, 'training_progress_gbm.png')
        plot_training_progress(model, save_path=training_plot_path)
        
        # 保存详细报告
        report_path = os.path.join(output_dir, 'gbm_report.txt')
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("GBM（梯度提升机）模型评估报告\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"准确率: {metrics['Accuracy']:.4f}\n")
            f.write(f"多类别AUC: {metrics['AUC']:.4f}\n")
            f.write(f"宏平均F1: {metrics['Macro F1']:.4f}\n")
            f.write(f"微平均F1: {metrics['Micro F1']:.4f}\n\n")
            
            f.write("模型参数:\n")
            f.write("-" * 20 + "\n")
            for param, value in model.get_params().items():
                f.write(f"{param}: {value}\n")
            f.write("\n")
            
            f.write("GBM特点:\n")
            f.write("-" * 20 + "\n")
            f.write("1. 序列化训练，每棵树纠正前面树的错误\n")
            f.write("2. 使用梯度下降优化损失函数\n")
            f.write("3. 支持多种损失函数和正则化\n")
            f.write("4. 对特征缩放不敏感\n")
            f.write("5. 提供特征重要性分析\n\n")
            
            f.write("前20个最重要特征:\n")
            f.write("-" * 30 + "\n")
            for i, (_, row) in enumerate(feature_importance.head(20).iterrows()):
                f.write(f"{i+1:2d}. {row['feature']:15s}: {row['importance']:.6f}\n")
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
        print(f"  特征重要性图: {importance_plot_path}")
        print(f"  训练过程图: {training_plot_path}")
        
    except Exception as e:
        print(f"❌ 保存结果失败: {str(e)}")
        traceback.print_exc()

def main():
    """主函数"""
    print("\n" + "="*60)
    print("--- 🚀 GBM（梯度提升机）模型训练与评估 ---")
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
    model, bayes_search = train_gbm(X_train, y_train, X_val, y_val)
    
    if model is None:
        print("❌ 模型训练失败，程序退出")
        return
    
    # 3. 评估模型
    metrics, y_pred, y_proba, feature_importance = evaluate_model(
        model, X_test, y_test, feature_names, class_names, figure_dir
    )
    
    if metrics is None:
        print("❌ 模型评估失败，程序退出")
        return
    
    # 4. 保存结果
    save_results(
        X_test, y_test, y_pred, names_test, 
        feature_names, metrics, model, feature_importance, output_dir, figure_dir
    )
    
    print(f"\n🎉 GBM模型训练和评估完成!")
    print(f"   最终测试准确率: {metrics['Accuracy']:.4f}")
    print(f"   模型包含 {model.n_estimators} 棵提升树")
    print(f"   学习率: {model.learning_rate:.4f}")

if __name__ == "__main__":
    main()
