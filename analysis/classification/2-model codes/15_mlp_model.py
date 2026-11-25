#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MLP (多层感知机) 模型训练脚本
从预处理好的数据中加载训练集和测试集，进行模型训练、调参和评估
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pickle
import sys
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import accuracy_score, roc_auc_score, f1_score, confusion_matrix, classification_report
import traceback
import matplotlib

matplotlib.rcParams['font.sans-serif'] = ['SimHei']  # 或 'Microsoft YaHei', 'Arial Unicode MS' 等
matplotlib.rcParams['axes.unicode_minus'] = False


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
        cmap='Greens',
        xticklabels=class_names,
        yticklabels=class_names,
        cbar_kws={'shrink': 0.6}
    )
    plt.xlabel('Predicted label', fontsize=11)
    plt.ylabel('True label', fontsize=11)
    # ensure title is consistent
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

def train_mlp(X_train, y_train, X_val, y_val):
    """训练MLP模型并进行网格搜索优化"""
    print("\n🚀 开始训练MLP模型...")
    
    try:
        # 定义基础模型
        model = MLPClassifier(
            random_state=42,
            early_stopping=True,
            validation_fraction=0.1,
            max_iter=1000
        )
        
        # 定义超参数搜索空间
        param_grid = {
            'hidden_layer_sizes': [(50,), (100,), (150,), (50, 50), (100, 50), (100, 100)],
            'activation': ['relu', 'tanh', 'logistic'],
            'solver': ['adam', 'sgd', 'lbfgs'],
            'alpha': [0.0001, 0.001, 0.01, 0.1],
            'learning_rate': ['constant', 'adaptive'],
        }
        
        # 计算参数组合总数
        param_count = (len(param_grid['hidden_layer_sizes']) * 
                      len(param_grid['activation']) * 
                      len(param_grid['solver']) * 
                      len(param_grid['alpha']) * 
                      len(param_grid['learning_rate']))
        
        print(f"🔍 开始网格搜索超参数优化 (总共 {param_count} 个参数组合, 5折交叉验证)...")
        
        # 网格搜索
        grid_search = GridSearchCV(
            estimator=model,
            param_grid=param_grid,
            cv=5,
            scoring='accuracy',
            n_jobs=-1,
            verbose=1
        )
        
        grid_search.fit(X_train, y_train)
        best_model = grid_search.best_estimator_
        
        print(f"✅ 模型训练完成!")
        print(f"最佳参数: {grid_search.best_params_}")
        print(f"最佳交叉验证得分: {grid_search.best_score_:.4f}")
        
        return best_model, grid_search
        
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
        cm_path = os.path.join(output_dir, 'confusion_matrix_mlp.png')
        plot_confusion_matrix(
            metrics['Confusion Matrix'], 
            class_names, 
            'MLP (ANOVA 100基因)',
            save_path=cm_path
        )
        
        return metrics, y_pred, y_proba
        
    except Exception as e:
        print(f"❌ 模型评估失败: {str(e)}")
        traceback.print_exc()
        return None, None, None

def analyze_training_history(model, output_dir):
    """分析训练历史（损失曲线）"""
    print("\n📈 分析MLP训练历史...")
    
    try:
        # 检查是否有损失曲线数据
        if hasattr(model, 'loss_curve_'):
            loss_curve = model.loss_curve_
            
            # 绘制损失曲线
            plt.figure(figsize=(10, 6))
            plt.plot(range(1, len(loss_curve) + 1), loss_curve, 'b-', linewidth=2, marker='o', markersize=4)
            plt.xlabel('迭代次数')
            plt.ylabel('损失值')
            plt.title('MLP训练损失曲线')
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            
            # 保存损失曲线图
            loss_base = os.path.join(output_dir, 'mlp_loss_curve')
            plt.savefig(loss_base + '.png', dpi=300, bbox_inches='tight')
            plt.savefig(loss_base + '.pdf', dpi=300, bbox_inches='tight')
            plt.show()
            
            # 保存损失曲线数据
            loss_data = pd.DataFrame({
                'iteration': range(1, len(loss_curve) + 1),
                'loss': loss_curve
            })
            
            loss_data_path = os.path.join(output_dir, 'mlp_loss_data.csv')
            loss_data.to_csv(loss_data_path, index=False, encoding='utf-8-sig')
            
            print(f"✅ 训练历史分析完成!")
            print(f"  损失 curve plots: {loss_base}.png and {loss_base}.pdf")
            print(f"  损失曲线数据: {loss_data_path}")
            print(f"  总迭代次数: {len(loss_curve)}")
            print(f"  最终损失值: {loss_curve[-1]:.6f}")
            print(f"  最小损失值: {min(loss_curve):.6f}")
            
            return loss_data
        else:
            print("⚠️ 模型没有损失曲线数据")
            return None
        
    except Exception as e:
        print(f"❌ 训练历史分析失败: {str(e)}")
        traceback.print_exc()
        return None

def analyze_network_structure(model):
    """分析网络结构"""
    print("\n🔍 分析MLP网络结构...")
    
    try:
        print(f"网络结构信息:")
        print(f"  隐藏层结构: {model.hidden_layer_sizes}")
        print(f"  激活函数: {model.activation}")
        print(f"  求解器: {model.solver}")
        print(f"  学习率: {model.learning_rate}")
        print(f"  正则化参数 (alpha): {model.alpha}")
        print(f"  总迭代次数: {model.n_iter_}")
        
        # 计算网络参数数量
        if hasattr(model, 'coefs_'):
            total_params = 0
            layer_info = []
            
            for i, coef in enumerate(model.coefs_):
                layer_params = coef.size + model.intercepts_[i].size
                total_params += layer_params
                layer_info.append({
                    'layer': f'Layer {i+1}',
                    'input_size': coef.shape[0],
                    'output_size': coef.shape[1],
                    'parameters': layer_params
                })
            
            print(f"  网络总参数数量: {total_params:,}")
            print(f"  各层参数详情:")
            for info in layer_info:
                print(f"    {info['layer']}: {info['input_size']} -> {info['output_size']}, 参数数: {info['parameters']:,}")
        
        return {
            'hidden_layer_sizes': model.hidden_layer_sizes,
            'activation': model.activation,
            'solver': model.solver,
            'learning_rate': model.learning_rate,
            'alpha': model.alpha,
            'n_iter': model.n_iter_,
            'total_params': total_params if 'total_params' in locals() else 0
        }
        
    except Exception as e:
        print(f"❌ 网络结构分析失败: {str(e)}")
        traceback.print_exc()
        return None

def save_results(X_test, y_test, y_pred, names_test, feature_names, metrics, model, grid_search, 
                loss_data, network_info, output_dir):
    """保存模型结果"""
    print("\n💾 正在保存结果...")
    
    try:
        # 保存测试集预测结果
        test_results_df = pd.DataFrame(X_test, columns=feature_names)
        test_results_df.insert(0, 'SampleID', names_test)
        test_results_df['True_label'] = y_test.values
        test_results_df['Predicted_label'] = y_pred
        
        results_path = os.path.join(output_dir, 'mlp_test_results.csv')
        test_results_df.to_csv(results_path, index=False, encoding='utf-8-sig')
        
        # 保存模型
        model_path = os.path.join(output_dir, 'mlp_model.pkl')
        with open(model_path, 'wb') as f:
            pickle.dump(model, f)
        
        # 保存评估指标
        metrics_path = os.path.join(output_dir, 'mlp_metrics.pkl')
        with open(metrics_path, 'wb') as f:
            pickle.dump(metrics, f)
        
        # 保存最佳参数
        params_path = os.path.join(output_dir, 'mlp_best_params.pkl')
        with open(params_path, 'wb') as f:
            pickle.dump(grid_search.best_params_, f)
        
        # 保存详细报告
        report_path = os.path.join(output_dir, 'mlp_report.txt')
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("MLP (多层感知机) 模型评估报告\n")
            f.write("=" * 50 + "\n\n")
            
            f.write("最佳超参数:\n")
            f.write("-" * 20 + "\n")
            for param, value in grid_search.best_params_.items():
                f.write(f"{param}: {value}\n")
            f.write(f"\n最佳交叉验证得分: {grid_search.best_score_:.4f}\n\n")
            
            if network_info:
                f.write("网络结构信息:\n")
                f.write("-" * 20 + "\n")
                f.write(f"隐藏层结构: {network_info.get('hidden_layer_sizes', 'N/A')}\n")
                f.write(f"激活函数: {network_info.get('activation', 'N/A')}\n")
                f.write(f"求解器: {network_info.get('solver', 'N/A')}\n")
                f.write(f"学习率策略: {network_info.get('learning_rate', 'N/A')}\n")
                f.write(f"正则化参数: {network_info.get('alpha', 'N/A')}\n")
                f.write(f"训练迭代次数: {network_info.get('n_iter', 'N/A')}\n")
                f.write(f"网络总参数数量: {network_info.get('total_params', 'N/A'):,}\n\n")
            
            f.write("测试集性能指标:\n")
            f.write("-" * 20 + "\n")
            f.write(f"准确率: {metrics['Accuracy']:.4f}\n")
            f.write(f"多类别AUC: {metrics['AUC']:.4f}\n")
            f.write(f"宏平均F1: {metrics['Macro F1']:.4f}\n")
            f.write(f"微平均F1: {metrics['Micro F1']:.4f}\n\n")
            
            if loss_data is not None:
                f.write("训练历史信息:\n")
                f.write("-" * 20 + "\n")
                f.write(f"总迭代次数: {len(loss_data)}\n")
                f.write(f"最终损失值: {loss_data['loss'].iloc[-1]:.6f}\n")
                f.write(f"最小损失值: {loss_data['loss'].min():.6f}\n\n")
            
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
    print("--- 🚀 MLP (多层感知机) 模型训练与评估 ---")
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
    model, grid_search = train_mlp(X_train, y_train, X_val, y_val)
    
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
    
    # 4. 分析训练历史
    loss_data = analyze_training_history(model, figure_dir)
    
    # 5. 分析网络结构
    network_info = analyze_network_structure(model)
    
    # 6. 保存结果
    save_results(
        X_test, y_test, y_pred, names_test, 
        feature_names, metrics, model, grid_search, 
        loss_data, network_info, output_dir
    )
    
    print(f"\n🎉 MLP模型训练和评估完成!")
    print(f"   最终测试准确率: {metrics['Accuracy']:.4f}")

if __name__ == "__main__":
    main()
