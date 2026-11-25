#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ResNet-50模型训练脚本
从预处理好的数据中加载训练集和测试集，进行ResNet-50深度学习模型训练和评估
"""

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pickle
import sys
from torchvision import models
from sklearn.metrics import accuracy_score, roc_auc_score, f1_score, confusion_matrix, classification_report
import traceback
from log_helper import setup_logging

# initialize logging
setup_logging(__file__)

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
    sns.heatmap(
        cm, 
        annot=True, 
        fmt='d', 
        cmap='Reds',
        xticklabels=class_names, 
        yticklabels=class_names,
    cbar_kws={'shrink': 0.6}
    )
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

def create_resnet50_model(input_features, num_classes=4, device='cpu'):
    """创建适用于表格数据的ResNet-50模型"""
    print("正在创建ResNet-50模型...")
    
    # 加载预训练的ResNet-50模型
    model = models.resnet50(weights=models.ResNet50_Weights.DEFAULT)
    
    # 修改第一层卷积层以适应1D特征数据
    # 将特征数据reshape为 (batch, 1, 1, features) 的形式
    model.conv1 = nn.Conv2d(
        in_channels=1,
        out_channels=64,
        kernel_size=(1, input_features),
        stride=1,
        padding=0,
        bias=False
    )
    
    # 修改最后的全连接层
    num_ftrs = model.fc.in_features
    model.fc = nn.Linear(num_ftrs, num_classes)
    
    # 移动到指定设备
    model = model.to(device)
    
    print(f"✅ ResNet-50模型创建完成")
    print(f"  输入特征数: {input_features}")
    print(f"  输出类别数: {num_classes}")
    print(f"  模型参数数量: {sum(p.numel() for p in model.parameters()):,}")
    print(f"  可训练参数数量: {sum(p.numel() for p in model.parameters() if p.requires_grad):,}")
    
    return model

def train_resnet50(X_train, y_train, X_val, y_val, feature_names, device='auto', 
                  epochs=30, batch_size=32, learning_rate=1e-4):
    """训练ResNet-50模型"""
    
    # 设备选择
    if device == 'auto':
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    print(f"\n🚀 开始训练ResNet-50模型...")
    print(f"使用设备: {device}")
    print(f"训练样本数: {len(X_train)}, 验证样本数: {len(X_val)}")
    
    try:
        # 将字符串标签转换为数字
        unique_labels = sorted(np.unique(np.concatenate([y_train, y_val])))
        label_to_idx = {label: idx for idx, label in enumerate(unique_labels)}
        
        y_train_numeric = np.array([label_to_idx[label] for label in y_train])
        y_val_numeric = np.array([label_to_idx[label] for label in y_val])
        
        # 转换为张量并reshape为适合ResNet的格式
        X_train_tensor = torch.FloatTensor(X_train).unsqueeze(1).unsqueeze(2)  # (N, 1, 1, features)
        X_val_tensor = torch.FloatTensor(X_val).unsqueeze(1).unsqueeze(2)
        
        y_train_tensor = torch.LongTensor(y_train_numeric)
        y_val_tensor = torch.LongTensor(y_val_numeric)
        
        print(f"输入数据形状: {X_train_tensor.shape} (N, C, H, W)")
        
        # 创建数据加载器
        train_dataset = TensorDataset(X_train_tensor, y_train_tensor)
        val_dataset = TensorDataset(X_val_tensor, y_val_tensor)
        
        train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
        val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
        
        # 创建模型
        num_classes = len(unique_labels)
        input_features = X_train.shape[1]
        model = create_resnet50_model(input_features, num_classes, device)
        
        # 损失函数和优化器
        criterion = nn.CrossEntropyLoss().to(device)
        optimizer = optim.AdamW(
            model.parameters(),
            lr=learning_rate,
            weight_decay=1e-5
        )
        scheduler = optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, mode='max', factor=0.5, patience=3, verbose=True
        )
        
        # 训练参数
        best_val_auc = 0.0
        best_model_state = None
        patience = 5
        no_improve_cnt = 0
        
        train_losses = []
        val_accuracies = []
        val_aucs = []
        
        print("\n开始训练...")
        for epoch in range(epochs):
            # 训练阶段
            model.train()
            train_loss = 0.0
            
            for batch_x, batch_y in train_loader:
                batch_x, batch_y = batch_x.to(device), batch_y.to(device)
                
                outputs = model(batch_x)
                loss = criterion(outputs, batch_y)
                
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
                
                train_loss += loss.item() * batch_x.size(0)
            
            avg_train_loss = train_loss / len(train_loader.dataset)
            train_losses.append(avg_train_loss)
            
            # 验证阶段
            model.eval()
            val_preds = []
            val_probs = []
            val_true = []
            
            with torch.no_grad():
                for batch_x, batch_y in val_loader:
                    batch_x, batch_y = batch_x.to(device), batch_y.to(device)
                    outputs = model(batch_x)
                    probs = torch.softmax(outputs, dim=1).cpu().numpy()
                    preds = torch.argmax(outputs, dim=1).cpu().numpy()
                    
                    val_probs.extend(probs)
                    val_preds.extend(preds)
                    val_true.extend(batch_y.cpu().numpy())
            
            val_metrics = calculate_metrics(val_true, val_preds, val_probs)
            val_auc = val_metrics['AUC']
            val_acc = val_metrics['Accuracy']
            
            val_accuracies.append(val_acc)
            val_aucs.append(val_auc)
            
            # 学习率调度
            scheduler.step(val_auc)
            
            # 早停检查
            if val_auc > best_val_auc:
                best_val_auc = val_auc
                best_model_state = model.state_dict()
                no_improve_cnt = 0
            else:
                no_improve_cnt += 1
                if no_improve_cnt >= patience:
                    print(f"  早停触发: 验证AUC连续{patience}轮未改善")
                    break
            
            # 打印训练进度
            if (epoch + 1) % 2 == 0:
                print(f"  轮次 {epoch+1}/{epochs} | 训练损失: {avg_train_loss:.4f} | "
                      f"验证准确率: {val_acc:.4f} | 验证AUC: {val_auc:.4f}")
        
        # 加载最佳模型
        model.load_state_dict(best_model_state)
        
        print(f"✅ ResNet-50训练完成!")
        print(f"最佳验证AUC: {best_val_auc:.4f}")
        
        return model, unique_labels, label_to_idx, train_losses, val_accuracies, val_aucs
        
    except Exception as e:
        print(f"❌ ResNet-50训练失败: {str(e)}")
        traceback.print_exc()
        return None, None, None, None, None, None

def evaluate_model(model, X_test, y_test, feature_names, class_names, 
                  unique_labels, label_to_idx, device, output_dir):
    """评估模型性能"""
    print("\n📊 开始模型评估...")
    
    try:
        # 将测试集标签转换为数字
        y_test_numeric = np.array([label_to_idx[label] for label in y_test])
        
        # 转换为张量
        X_test_tensor = torch.FloatTensor(X_test).unsqueeze(1).unsqueeze(2)
        y_test_tensor = torch.LongTensor(y_test_numeric)
        
        # 创建测试数据加载器
        test_dataset = TensorDataset(X_test_tensor, y_test_tensor)
        test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False)
        
        # 测试阶段
        model.eval()
        test_preds = []
        test_probs = []
        test_true = []
        
        with torch.no_grad():
            for batch_x, batch_y in test_loader:
                batch_x, batch_y = batch_x.to(device), batch_y.to(device)
                outputs = model(batch_x)
                probs = torch.softmax(outputs, dim=1).cpu().numpy()
                preds = torch.argmax(outputs, dim=1).cpu().numpy()
                
                test_probs.extend(probs)
                test_preds.extend(preds)
                test_true.extend(batch_y.cpu().numpy())
        
        # 将数字预测结果转换回原始标签
        idx_to_label = {idx: label for label, idx in label_to_idx.items()}
        test_preds_original = [idx_to_label[pred] for pred in test_preds]
        test_true_original = [idx_to_label[true] for true in test_true]
        
        test_metrics = calculate_metrics(test_true, test_preds, test_probs)
        
        print(f"\n🎯 测试集性能指标:")
        print(f"  准确率: {test_metrics['Accuracy']:.4f}")
        print(f"  多类别AUC: {test_metrics['AUC']:.4f}")
        print(f"  宏平均F1: {test_metrics['Macro F1']:.4f}")
        print(f"  微平均F1: {test_metrics['Micro F1']:.4f}")
        
        # 绘制混淆矩阵
        cm_path = os.path.join(output_dir, 'confusion_matrix_resnet50.png')
        plot_confusion_matrix(
            test_metrics['Confusion Matrix'], 
            class_names, 
            'ResNet-50 (ANOVA 100基因)',
            save_path=cm_path
        )
        
        return test_metrics, test_preds_original, test_true_original
        
    except Exception as e:
        print(f"❌ 模型评估失败: {str(e)}")
        traceback.print_exc()
        return None, None, None

def plot_training_history(train_losses, val_accuracies, val_aucs, save_path=None):
    """绘制训练历史"""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # 训练损失
    axes[0].plot(train_losses, 'b-', linewidth=2)
    axes[0].set_title('训练损失')
    axes[0].set_xlabel('轮次')
    axes[0].set_ylabel('损失值')
    axes[0].grid(True, alpha=0.3)
    
    # 验证准确率
    axes[1].plot(val_accuracies, 'g-', linewidth=2)
    axes[1].set_title('验证准确率')
    axes[1].set_xlabel('轮次')
    axes[1].set_ylabel('准确率')
    axes[1].grid(True, alpha=0.3)
    
    # 验证AUC
    axes[2].plot(val_aucs, 'r-', linewidth=2)
    axes[2].set_title('验证AUC')
    axes[2].set_xlabel('轮次')
    axes[2].set_ylabel('AUC值')
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"训练历史图已保存到: {save_path}")
    
    plt.show()

def save_results(model, X_test, y_test_true, y_test_pred, names_test, feature_names, 
                metrics, train_losses, val_accuracies, val_aucs, output_dir, figure_dir):
    """保存模型结果"""
    print("\n💾 正在保存结果...")
    
    try:
        # 保存测试集预测结果
        test_results_df = pd.DataFrame(X_test, columns=feature_names)
        test_results_df.insert(0, 'SampleID', names_test)
        test_results_df['True_label'] = y_test_true
        test_results_df['Predicted_label'] = y_test_pred
        
        results_path = os.path.join(output_dir, 'resnet50_test_results.csv')
        test_results_df.to_csv(results_path, index=False, encoding='utf-8-sig')
        
        # 保存模型
        model_path = os.path.join(output_dir, 'resnet50_model.pth')
        torch.save(model.state_dict(), model_path)
        
        # 保存完整模型（包含架构）
        full_model_path = os.path.join(output_dir, 'resnet50_full_model.pth')
        torch.save(model, full_model_path)
        
        # 保存评估指标
        metrics_path = os.path.join(output_dir, 'resnet50_metrics.pkl')
        with open(metrics_path, 'wb') as f:
            pickle.dump(metrics, f)
        
        # 绘制并保存训练历史图
        history_plot_path = os.path.join(figure_dir, 'resnet50_training_history.png')
        plot_training_history(train_losses, val_accuracies, val_aucs, save_path=history_plot_path)
        
        # 保存详细报告
        report_path = os.path.join(output_dir, 'resnet50_report.txt')
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("ResNet-50模型评估报告\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"准确率: {metrics['Accuracy']:.4f}\n")
            f.write(f"多类别AUC: {metrics['AUC']:.4f}\n")
            f.write(f"宏平均F1: {metrics['Macro F1']:.4f}\n")
            f.write(f"微平均F1: {metrics['Micro F1']:.4f}\n\n")
            
            f.write("模型架构:\n")
            f.write("-" * 20 + "\n")
            f.write("基础架构: ResNet-50\n")
            f.write("预训练权重: ImageNet\n")
            f.write("输入适配: 1D特征 -> 2D卷积\n")
            f.write(f"输入特征数: {len(feature_names)}\n")
            f.write(f"输出类别数: {len(np.unique(y_test_true))}\n")
            f.write(f"总参数数量: {sum(p.numel() for p in model.parameters()):,}\n")
            f.write(f"可训练参数: {sum(p.numel() for p in model.parameters() if p.requires_grad):,}\n\n")
            
            f.write("ResNet-50特点:\n")
            f.write("-" * 20 + "\n")
            f.write("1. 深度残差网络，解决梯度消失问题\n")
            f.write("2. 使用跳跃连接（skip connections）\n")
            f.write("3. 预训练权重提供强大的特征提取能力\n")
            f.write("4. 适配表格数据的1D->2D转换\n")
            f.write("5. 端到端训练，自动特征学习\n")
            f.write("6. 早停和学习率调度防止过拟合\n\n")
            
            f.write("训练配置:\n")
            f.write("-" * 20 + "\n")
            f.write("优化器: AdamW\n")
            f.write("损失函数: CrossEntropyLoss\n")
            f.write("学习率调度: ReduceLROnPlateau\n")
            f.write("早停策略: 验证AUC连续5轮未改善\n")
            f.write(f"最终训练轮次: {len(train_losses)}\n")
            f.write(f"最佳验证AUC: {max(val_aucs):.4f}\n\n")
            
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
        print(f"  模型权重: {model_path}")
        print(f"  完整模型: {full_model_path}")
        print(f"  评估指标: {metrics_path}")
        print(f"  详细报告: {report_path}")
        print(f"  训练历史图: {history_plot_path}")
        
    except Exception as e:
        print(f"❌ 保存结果失败: {str(e)}")
        traceback.print_exc()

def main():
    """主函数"""
    print("\n" + "="*60)
    print("--- 🚀 ResNet-50模型训练与评估 ---")
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
    model, unique_labels, label_to_idx, train_losses, val_accuracies, val_aucs = train_resnet50(
        X_train, y_train, X_val, y_val, feature_names
    )
    
    if model is None:
        print("❌ 模型训练失败，程序退出")
        return
    
    # 3. 评估模型
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    metrics, y_pred, y_true = evaluate_model(
        model, X_test, y_test, feature_names, class_names,
        unique_labels, label_to_idx, device, figure_dir
    )
    
    if metrics is None:
        print("❌ 模型评估失败，程序退出")
        return
    
    # 4. 保存结果
    save_results(
        model, X_test, y_true, y_pred, names_test,
        feature_names, metrics, train_losses, val_accuracies, val_aucs,
        output_dir, figure_dir
    )
    
    print(f"\n🎉 ResNet-50模型训练和评估完成!")
    print(f"   最终测试准确率: {metrics['Accuracy']:.4f}")
    print(f"   模型参数数量: {sum(p.numel() for p in model.parameters()):,}")

if __name__ == "__main__":
    main()
