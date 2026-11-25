#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mamba模型训练脚本
基于TabMamba架构的表格数据分类模型
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
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report, roc_auc_score, f1_score
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
    # square figure, English labels, and smaller colorbar
    fig, ax = plt.subplots(figsize=(6, 6))
    # use a different colormap than the logistic model and labeled classes
    cmap_choice = 'Reds'  # single-color red palette
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
    ax.set_xlabel('Predicted label', fontsize=11)
    ax.set_ylabel('True label', fontsize=11)
    ax.set_title(f'TabMamba (ANOVA 100 genes) Confusion Matrix', fontsize=12, pad=15)
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

class SimplifiedMambaLayer(nn.Module):
    """简化的Mamba层实现"""
    def __init__(self, d_model, d_state=16, dropout=0.1):
        super().__init__()
        self.d_model = d_model
        self.d_state = d_state
        
        self.input_proj = nn.Linear(d_model, d_model)
        self.A = nn.Parameter(torch.ones(1, d_state) * -1.0)
        self.B = nn.Linear(d_model, d_state)
        self.C = nn.Linear(d_state, d_model)
        self.D = nn.Parameter(torch.zeros(1, d_model))
        self.silu = nn.SiLU()
        self.dropout = nn.Dropout(dropout)
        
    def forward(self, x):
        x_proj = self.input_proj(x)
        x_act = self.silu(x_proj)
        
        state = self.B(x_act)
        out = self.C(state) + self.D * x_act
        out = self.dropout(out)
        
        return out

class TabMamba(nn.Module):
    """适用于表格数据的Mamba模型"""
    def __init__(self, input_dim=100, d_model=64, n_mamba_layers=1, dropout=0.2, num_classes=4):
        super().__init__()
        self.input_dim = input_dim
        self.d_model = d_model
        
        # 特征投影层
        self.feat_proj = nn.Sequential(
            nn.Linear(input_dim, d_model),
            nn.BatchNorm1d(d_model),
            nn.SiLU(),
            nn.Dropout(dropout)
        )
        
        # Mamba层
        self.mamba_layers = nn.ModuleList([
            SimplifiedMambaLayer(
                d_model=d_model,
                d_state=16,
                dropout=dropout
            ) for _ in range(n_mamba_layers)
        ])
        
        # 输出头
        self.output_head = nn.Sequential(
            nn.Linear(d_model, d_model // 2),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(d_model // 2, num_classes)
        )
        
    def forward(self, x):
        x = self.feat_proj(x)
        for mamba_layer in self.mamba_layers:
            x = mamba_layer(x)
        x = self.output_head(x)
        return x

def train_mamba_model(X_train, y_train, X_val, y_val, X_test, y_test, 
                     feature_names, device='auto', epochs=50, batch_size=32):
    """训练Mamba模型"""
    
    # 设备选择
    if device == 'auto':
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    print(f"\n🚀 开始训练Mamba模型...")
    print(f"使用设备: {device}")
    print(f"训练样本数: {len(X_train)}, 测试样本数: {len(X_test)}")
    
    try:
        # 转换为张量
        X_train_tensor = torch.FloatTensor(X_train)
        X_val_tensor = torch.FloatTensor(X_val)
        X_test_tensor = torch.FloatTensor(X_test)
        
        # 将字符串标签转换为数字
        unique_labels = sorted(np.unique(np.concatenate([y_train, y_val, y_test])))
        label_to_idx = {label: idx for idx, label in enumerate(unique_labels)}
        
        y_train_numeric = np.array([label_to_idx[label] for label in y_train])
        y_val_numeric = np.array([label_to_idx[label] for label in y_val])
        y_test_numeric = np.array([label_to_idx[label] for label in y_test])
        
        y_train_tensor = torch.LongTensor(y_train_numeric)
        y_val_tensor = torch.LongTensor(y_val_numeric)
        y_test_tensor = torch.LongTensor(y_test_numeric)
        
        # 创建数据加载器
        train_dataset = TensorDataset(X_train_tensor, y_train_tensor)
        val_dataset = TensorDataset(X_val_tensor, y_val_tensor)
        test_dataset = TensorDataset(X_test_tensor, y_test_tensor)
        
        train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
        val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
        test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)
        
        # 创建模型
        num_classes = len(unique_labels)
        model = TabMamba(
            input_dim=X_train_tensor.size(1),
            d_model=64,
            n_mamba_layers=1,
            dropout=0.2,
            num_classes=num_classes
        ).to(device)
        
        # 损失函数和优化器
        criterion = nn.CrossEntropyLoss().to(device)
        optimizer = optim.AdamW(
            model.parameters(),
            lr=5e-4,
            weight_decay=1e-4
        )
        scheduler = optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, mode='max', factor=0.5, patience=3, verbose=True
        )
        
        # 训练参数
        best_val_auc = 0.0
        best_model_state = None
        patience = 5
        no_improve_cnt = 0
        
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
        
        # 将数字标签转换回原始标签
        test_preds_original = [unique_labels[pred] for pred in test_preds]
        test_true_original = [unique_labels[true] for true in test_true]
        
        test_metrics = calculate_metrics(test_true, test_preds, test_probs)
        
        print(f"\n✅ Mamba训练完成 - 关键测试集指标:")
        print(f"  准确率: {test_metrics['Accuracy']:.4f}")
        print(f"  多类别AUC: {test_metrics['AUC']:.4f}")
        print(f"  宏平均F1: {test_metrics['Macro F1']:.4f}")
        
        return model, test_metrics, test_preds_original, test_true_original, unique_labels
        
    except Exception as e:
        print(f"❌ Mamba训练失败: {str(e)}")
        traceback.print_exc()
        return None, None, None, None, None

def save_results(model, metrics, X_test, y_test_true, y_test_pred, names_test, 
                feature_names, class_names, output_dir, figure_dir):
    """保存模型结果"""
    print("\n💾 正在保存结果...")
    
    try:
        # 保存测试集预测结果
        test_results_df = pd.DataFrame(X_test, columns=feature_names)
        test_results_df.insert(0, 'SampleID', names_test)
        test_results_df['True_label'] = y_test_true
        test_results_df['Predicted_label'] = y_test_pred
        
        results_path = os.path.join(output_dir, 'mamba_test_results.csv')
        test_results_df.to_csv(results_path, index=False, encoding='utf-8-sig')
        
        # 保存模型
        model_path = os.path.join(output_dir, 'mamba_model.pth')
        torch.save(model.state_dict(), model_path)
        
        # 保存完整模型（包含架构）
        full_model_path = os.path.join(output_dir, 'mamba_full_model.pth')
        torch.save(model, full_model_path)
        
        # 保存评估指标
        metrics_path = os.path.join(output_dir, 'mamba_metrics.pkl')
        with open(metrics_path, 'wb') as f:
            pickle.dump(metrics, f)
        
        # 绘制并保存混淆矩阵
        cm_path = os.path.join(figure_dir, 'confusion_matrix_mamba.pdf')
        plot_confusion_matrix(
            metrics['Confusion Matrix'],
            class_names,
            'TabMamba',
            save_path=cm_path
        )
        
        # 保存详细报告
        report_path = os.path.join(output_dir, 'mamba_report.txt')
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("TabMamba模型评估报告\n")
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
        print(f"  模型权重: {model_path}")
        print(f"  完整模型: {full_model_path}")
        print(f"  评估指标: {metrics_path}")
        print(f"  详细报告: {report_path}")
        print(f"  混淆矩阵: {cm_path}")
        
    except Exception as e:
        print(f"❌ 保存结果失败: {str(e)}")
        traceback.print_exc()

def main():
    """主函数"""
    print("\n" + "="*60)
    print("--- 🚀 TabMamba模型训练与评估 ---")
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
    
    # 2. 训练模型
    model, metrics, y_pred, y_true, class_names = train_mamba_model(
        X_train, y_train, X_val, y_val, X_test, y_test, feature_names
    )
    
    if model is None:
        print("❌ 模型训练失败，程序退出")
        return
    
    # 3. 保存结果
    save_results(
        model, metrics, X_test, y_true, y_pred, names_test,
        feature_names, class_names, output_dir, figure_dir
    )
    
    print(f"\n🎉 TabMamba模型训练和评估完成!")
    print(f"   最终测试准确率: {metrics['Accuracy']:.4f}")

if __name__ == "__main__":
    main()
