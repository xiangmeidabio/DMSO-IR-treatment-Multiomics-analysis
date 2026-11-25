#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AdaBoost模型训练脚本
从预处理好的数据中加载训练集和测试集，进行模型训练、调参和评估
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pickle
import sys
from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from skopt.space import Integer, Real, Categorical
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
    # make square figure and use a single-color colormap (blue)
    plt.figure(figsize=(6, 6))
    ax = sns.heatmap(
        cm,
        annot=True,
        fmt='d',
        cmap='Blues',
        xticklabels=class_names,
        yticklabels=class_names,
        cbar_kws={'shrink': 0.6}  # smaller colorbar
    )
    # remove colorbar label and set tick size
    try:
        cbar = ax.collections[0].colorbar
        cbar.set_label('')
        cbar.ax.tick_params(labelsize=8)
    except Exception:
        pass
    # English labels
    plt.xlabel('Predicted label')
    plt.ylabel('True label')
    plt.title(f'{model_name} Confusion Matrix')
    plt.tight_layout()
    
    if save_path:
        # save both png and pdf as square
        base, ext = os.path.splitext(save_path)
        png_path = base + '.png'
        pdf_path = base + '.pdf'
        plt.savefig(png_path, dpi=300, bbox_inches='tight')
        plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
        print(f"Confusion matrix saved to: {png_path} and {pdf_path}")
    
    plt.show()

def train_adaboost(X_train, y_train, X_val, y_val):
    """训练AdaBoost模型并进行贝叶斯优化"""
    print("\n🚀 开始训练AdaBoost模型...")
    
    try:
        # 定义基础模型 - 使用决策树桩（max_depth=1）作为弱学习器
        base_estimator = DecisionTreeClassifier(max_depth=1, random_state=42)
        
        # 定义AdaBoost模型
        model = AdaBoostClassifier(
            estimator=base_estimator,
            random_state=42,
            algorithm='SAMME'  # 适用于多类别分类
        )
        
        # 定义超参数搜索空间
        params = {
            'n_estimators': Integer(50, 200),
            'learning_rate': Real(0.01, 2.0, prior='log-uniform'),
            'algorithm': Categorical(['SAMME'])  # 只保留SAMME，SAMME.R已被弃用
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
        cm_path = os.path.join(output_dir, 'confusion_matrix_adaboost.png')
        plot_confusion_matrix(
            metrics['Confusion Matrix'], 
            class_names, 
            'AdaBoost (ANOVA 100 genes)',
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
        importance_path = os.path.join(output_dir, 'adaboost_feature_importance.csv')
        feature_importance_df.to_csv(importance_path, index=False, encoding='utf-8-sig')
        # plot top 20 features (single color, square) and save PNG+PDF
        plt.figure(figsize=(6, 6))
        top_features = feature_importance_df.head(20)
        sns.barplot(data=top_features, x='importance', y='feature', color='blue')
        plt.title('AdaBoost - Top 20 important features')
        plt.xlabel('Feature importance')
        plt.ylabel('Gene')
        plt.tight_layout()

        importance_base = os.path.join(output_dir, 'adaboost_feature_importance')
        importance_png = importance_base + '.png'
        importance_pdf = importance_base + '.pdf'
        plt.savefig(importance_png, dpi=300, bbox_inches='tight')
        plt.savefig(importance_pdf, dpi=300, bbox_inches='tight')
        plt.show()

        print(f"✅ 特征重要性分析完成!")
        print(f"  特征重要性文件: {importance_path}")
        print(f"  特征重要性图像: {importance_png} and {importance_pdf}")
        print(f"\n前10个最重要特征:")
        for i, (_, row) in enumerate(feature_importance_df.head(10).iterrows(), 1):
            print(f"  {i:2d}. {row['feature']}: {row['importance']:.4f}")

        return feature_importance_df

    except Exception as e:
        print(f"❌ 特征重要性分析失败: {str(e)}")
        traceback.print_exc()
        return None

def analyze_learning_curve(model, output_dir):
    """分析学习曲线（弱学习器的权重）"""
    print("\n📈 分析AdaBoost学习曲线...")
    
    try:
        # 获取弱学习器权重
        estimator_weights = model.estimator_weights_
        estimator_errors = model.estimator_errors_
        
        # plot learning curve (square) and save PNG+PDF
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 8))

        # weak learner weights
        ax1.plot(range(1, len(estimator_weights) + 1), estimator_weights, 'b-', marker='o')
        ax1.set_xlabel('Weak learner index')
        ax1.set_ylabel('Weight')
        ax1.set_title('AdaBoost Weak Learner Weights')
        ax1.grid(True, alpha=0.3)

        # weak learner errors
        ax2.plot(range(1, len(estimator_errors) + 1), estimator_errors, 'r-', marker='s')
        ax2.set_xlabel('Weak learner index')
        ax2.set_ylabel('Error')
        ax2.set_title('AdaBoost Weak Learner Errors')
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()

        learning_base = os.path.join(output_dir, 'adaboost_learning_curve')
        learning_png = learning_base + '.png'
        learning_pdf = learning_base + '.pdf'
        plt.savefig(learning_png, dpi=300, bbox_inches='tight')
        plt.savefig(learning_pdf, dpi=300, bbox_inches='tight')
        plt.show()

        # 保存学习曲线数据
        learning_data = pd.DataFrame({
            'estimator_index': range(1, len(estimator_weights) + 1),
            'weight': estimator_weights,
            'error': estimator_errors
        })

        learning_data_path = os.path.join(output_dir, 'adaboost_learning_data.csv')
        learning_data.to_csv(learning_data_path, index=False, encoding='utf-8-sig')

        print(f"✅ 学习曲线分析完成!")
        print(f"  学习曲线图像: {learning_png} and {learning_pdf}")
        print(f"  学习曲线数据: {learning_data_path}")
        print(f"  弱学习器总数: {len(estimator_weights)}")
        print(f"  平均权重: {np.mean(estimator_weights):.4f}")
        print(f"  平均错误率: {np.mean(estimator_errors):.4f}")

        return learning_data

    except Exception as e:
        print(f"❌ 学习曲线分析失败: {str(e)}")
        traceback.print_exc()
        return None

def save_results(X_test, y_test, y_pred, names_test, feature_names, metrics, model, bayes_search, 
                feature_importance_df, learning_data, output_dir):
    """保存模型结果"""
    print("\n💾 正在保存结果...")
    
    try:
        # 保存测试集预测结果
        test_results_df = pd.DataFrame(X_test, columns=feature_names)
        test_results_df.insert(0, 'SampleID', names_test)
        test_results_df['True_label'] = y_test.values
        test_results_df['Predicted_label'] = y_pred
        
        results_path = os.path.join(output_dir, 'adaboost_test_results.csv')
        test_results_df.to_csv(results_path, index=False, encoding='utf-8-sig')
        
        # 保存模型
        model_path = os.path.join(output_dir, 'adaboost_model.pkl')
        with open(model_path, 'wb') as f:
            pickle.dump(model, f)
        
        # 保存评估指标
        metrics_path = os.path.join(output_dir, 'adaboost_metrics.pkl')
        with open(metrics_path, 'wb') as f:
            pickle.dump(metrics, f)
        
        # 保存最佳参数
        params_path = os.path.join(output_dir, 'adaboost_best_params.pkl')
        with open(params_path, 'wb') as f:
            pickle.dump(bayes_search.best_params_, f)
        
        # 保存详细报告
        report_path = os.path.join(output_dir, 'adaboost_report.txt')
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("AdaBoost模型评估报告\n")
            f.write("=" * 50 + "\n\n")
            
            f.write("最佳超参数:\n")
            f.write("-" * 20 + "\n")
            for param, value in bayes_search.best_params_.items():
                f.write(f"{param}: {value}\n")
            f.write(f"\n最佳交叉验证得分: {bayes_search.best_score_:.4f}\n\n")
            
            f.write("模型信息:\n")
            f.write("-" * 20 + "\n")
            f.write(f"弱学习器数量: {len(model.estimator_weights_)}\n")
            f.write(f"平均弱学习器权重: {np.mean(model.estimator_weights_):.4f}\n")
            f.write(f"平均弱学习器错误率: {np.mean(model.estimator_errors_):.4f}\n\n")
            
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
        print(f"  格式说明: 第一列=样本名称, 中间列=特征, 倒数第二列=真实标签, 最后一列=预测标签")
        
    except Exception as e:
        print(f"❌ 保存结果失败: {str(e)}")
        traceback.print_exc()

def main():
    """主函数"""
    print("\n" + "="*60)
    print("--- 🚀 AdaBoost模型训练与评估 ---")
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
    model, bayes_search = train_adaboost(X_train, y_train, X_val, y_val)
    
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
    
    # 5. 分析学习曲线
    learning_data = analyze_learning_curve(model, figure_dir)
    
    # 6. 保存结果
    save_results(
        X_test, y_test, y_pred, names_test, 
        feature_names, metrics, model, bayes_search, 
        feature_importance_df, learning_data, output_dir
    )
    
    print(f"\n🎉 AdaBoost模型训练和评估完成!")
    print(f"   最终测试准确率: {metrics['Accuracy']:.4f}")

if __name__ == "__main__":
    main()
