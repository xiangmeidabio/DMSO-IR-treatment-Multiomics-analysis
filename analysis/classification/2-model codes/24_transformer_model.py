#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Transformer (TabTransformer) 模型训练脚本
使用专门为表格数据设计的Transformer架构进行基因表达数据分类
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
import math
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
    sns.heatmap(cm, annot=True, fmt='d', cmap='Greens', 
                xticklabels=class_names, yticklabels=class_names,
                cbar_kws={'shrink': 0.6})
    plt.xlabel('Predicted label', fontsize=11)
    plt.ylabel('True label', fontsize=11)
    plt.title(f'{model_name} Confusion Matrix', fontsize=12, pad=15)
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

class PositionalEncoding(nn.Module):
    """位置编码模块，为序列数据添加位置信息"""
    
    def __init__(self, d_model, max_len=100, dropout=0.1):
        super().__init__()
        self.dropout = nn.Dropout(p=dropout)
        
        # 创建位置编码矩阵
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * (-math.log(10000.0) / d_model))
        
        # 使用正弦和余弦函数生成位置编码
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0)
        
        # 注册为缓冲区，不参与梯度更新
        self.register_buffer('pe', pe)

    def forward(self, x):
        # 添加位置编码并应用dropout
        x = x + self.pe[:, :x.size(1), :]
        return self.dropout(x)

class TabTransformerFull(nn.Module):
    """专门为表格数据设计的Transformer模型"""
    
    def __init__(self, input_dim, num_classes=4, d_model=64, nhead=4, num_layers=2, 
                 dim_feedforward=128, dropout=0.2):
        super().__init__()
        self.input_dim = input_dim
        self.d_model = d_model
        
        # 输入投影层：将原始特征映射到模型维度
        self.input_projection = nn.Sequential(
            nn.Linear(input_dim, d_model),
            nn.BatchNorm1d(d_model),
            nn.ReLU(),
            nn.Dropout(dropout)
        )
        
        # 位置编码
        self.pos_encoder = PositionalEncoding(d_model, max_len=1, dropout=dropout)
        
        # Transformer编码器层
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=nhead,
            dim_feedforward=dim_feedforward,
            dropout=dropout,
            batch_first=True  # 使用batch_first格式
        )
        self.transformer_encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)
        
        # 分类头
        self.fc = nn.Sequential(
            nn.Linear(d_model, d_model//2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(d_model//2, num_classes)
        )

    def forward(self, x):
        batch_size = x.size(0)
        
        # 输入投影
        x = self.input_projection(x)
        
        # 添加序列维度 (batch_size, seq_len=1, d_model)
        x = x.unsqueeze(1)
        
        # 位置编码
        x = self.pos_encoder(x)
        
        # Transformer编码
        x = self.transformer_encoder(x)
        
        # 移除序列维度
        x = x.squeeze(1)
        
        # 分类
        x = self.fc(x)
        return x

def prepare_data_loaders(train_data, val_data, test_data, batch_size=32):
    """准备PyTorch数据加载器"""
    print("正在准备PyTorch数据加载器...")
    
    try:
        # 转换为numpy数组
        X_train_np = train_data['X_train']
        X_val_np = val_data['X_val']
        X_test_np = test_data['X_test']
        
        y_train_np = train_data['y_train'].astype(int)
        y_val_np = val_data['y_val'].astype(int)
        y_test_np = test_data['y_test'].astype(int)
        
        # 转换为PyTorch张量
        X_train_tensor = torch.FloatTensor(X_train_np)
        X_val_tensor = torch.FloatTensor(X_val_np)
        X_test_tensor = torch.FloatTensor(X_test_np)
        
        y_train_tensor = torch.LongTensor(y_train_np)
        y_val_tensor = torch.LongTensor(y_val_np)
        y_test_tensor = torch.LongTensor(y_test_np)
        
        print(f"数据张量形状: {X_train_tensor.shape}")
        
        # 创建数据集
        train_dataset = TensorDataset(X_train_tensor, y_train_tensor)
        val_dataset = TensorDataset(X_val_tensor, y_val_tensor)
        test_dataset = TensorDataset(X_test_tensor, y_test_tensor)
        
        # 创建数据加载器
        train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
        val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
        test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)
        
        print(f"✅ 数据加载器准备完成!")
        print(f"  批次大小: {batch_size}")
        print(f"  训练批次数: {len(train_loader)}")
        print(f"  验证批次数: {len(val_loader)}")
        print(f"  测试批次数: {len(test_loader)}")
        
        return train_loader, val_loader, test_loader, X_train_tensor.shape[1]
        
    except Exception as e:
        print(f"❌ 数据加载器准备失败: {str(e)}")
        traceback.print_exc()
        return None, None, None, None

def create_transformer_model(input_dim, num_classes=4):
    """创建TabTransformer模型"""
    print("\n🏗️ 创建TabTransformer模型...")
    
    try:
        model = TabTransformerFull(
            input_dim=input_dim,
            num_classes=num_classes,
            d_model=64,           # 模型维度
            nhead=4,              # 注意力头数
            num_layers=2,         # Transformer层数
            dim_feedforward=128,  # 前馈网络维度
            dropout=0.2           # Dropout率
        )
        
        print(f"✅ TabTransformer模型创建完成!")
        print(f"  输入维度: {input_dim}")
        print(f"  模型维度: 64")
        print(f"  注意力头数: 4")
        print(f"  Transformer层数: 2")
        print(f"  前馈网络维度: 128")
        print(f"  输出类别数: {num_classes}")
        
        return model
        
    except Exception as e:
        print(f"❌ 模型创建失败: {str(e)}")
        traceback.print_exc()
        return None

def train_model(model, train_loader, val_loader, device, epochs=30, patience=5):
    """训练Transformer模型"""
    print("\n🚀 开始训练TabTransformer模型...")
    
    try:
        # 损失函数和优化器
        criterion = nn.CrossEntropyLoss().to(device)
        optimizer = optim.AdamW(
            model.parameters(),
            lr=3e-4,              # 较高的学习率，适合Transformer
            weight_decay=1e-4     # 较强的权重衰减
        )
        scheduler = optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, mode='max', factor=0.5, patience=3, verbose=True
        )
        
        # 训练记录
        train_losses = []
        val_accuracies = []
        val_aucs = []
        
        best_val_auc = 0.0
        best_model_state = None
        no_improve_cnt = 0
        
        print(f"训练配置:")
        print(f"  训练轮数: {epochs}")
        print(f"  早停耐心: {patience}")
        print(f"  学习率: 3e-4 (适合Transformer)")
        print(f"  权重衰减: 1e-4")
        print(f"  设备: {device}")
        
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
                
                # 梯度裁剪，防止梯度爆炸
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                
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
                best_model_state = model.state_dict().copy()
                no_improve_cnt = 0
            else:
                no_improve_cnt += 1
                if no_improve_cnt >= patience:
                    print(f"  早停触发: 验证AUC连续{patience}轮无改善")
                    break
            
            # 打印进度
            if (epoch + 1) % 2 == 0 or epoch == 0:
                print(f"  轮次 {epoch+1}/{epochs} | 训练损失: {avg_train_loss:.4f} | "
                      f"验证准确率: {val_acc:.4f} | 验证AUC: {val_auc:.4f}")
        
        # 加载最佳模型
        model.load_state_dict(best_model_state)
        
        print(f"✅ 模型训练完成!")
        print(f"  最佳验证AUC: {best_val_auc:.4f}")
        print(f"  实际训练轮数: {epoch + 1}")
        
        # 返回训练历史
        training_history = {
            'train_losses': train_losses,
            'val_accuracies': val_accuracies,
            'val_aucs': val_aucs,
            'best_val_auc': best_val_auc,
            'epochs_trained': epoch + 1
        }
        
        return model, training_history
        
    except Exception as e:
        print(f"❌ 模型训练失败: {str(e)}")
        traceback.print_exc()
        return None, None

def evaluate_model(model, test_loader, device, feature_names, class_names, output_dir):
    """评估模型性能"""
    print("\n📊 开始模型评估...")
    
    try:
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
        
        # 计算指标
        metrics = calculate_metrics(test_true, test_preds, test_probs)
        
        # 打印结果
        print(f"\n🎯 测试集性能指标:")
        print(f"  准确率 (Accuracy): {metrics['Accuracy']:.4f}")
        print(f"  多类别AUC: {metrics['AUC']:.4f}")
        print(f"  宏平均F1: {metrics['Macro F1']:.4f}")
        print(f"  微平均F1: {metrics['Micro F1']:.4f}")
        
        # 绘制混淆矩阵
        cm_path = os.path.join(output_dir, 'confusion_matrix_transformer.png')
        plot_confusion_matrix(
            metrics['Confusion Matrix'], 
            class_names, 
            'TabTransformer (ANOVA 100 genes)',
            save_path=cm_path
        )
        
        return metrics, test_preds, test_probs
        
    except Exception as e:
        print(f"❌ 模型评估失败: {str(e)}")
        traceback.print_exc()
        return None, None, None

def plot_training_history(training_history, output_dir):
    """绘制训练历史曲线"""
    print("\n📈 绘制训练历史曲线...")
    
    try:
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
        
        epochs = range(1, len(training_history['train_losses']) + 1)
        
        # 训练损失
        ax1.plot(epochs, training_history['train_losses'], 'b-', linewidth=2, marker='o', markersize=4)
        ax1.set_xlabel('训练轮次')
        ax1.set_ylabel('训练损失')
        ax1.set_title('TabTransformer 训练损失曲线')
        ax1.grid(True, alpha=0.3)
        
        # 验证准确率
        ax2.plot(epochs, training_history['val_accuracies'], 'g-', linewidth=2, marker='s', markersize=4)
        ax2.set_xlabel('训练轮次')
        ax2.set_ylabel('验证准确率')
        ax2.set_title('TabTransformer 验证准确率曲线')
        ax2.grid(True, alpha=0.3)
        
        # 验证AUC
        ax3.plot(epochs, training_history['val_aucs'], 'r-', linewidth=2, marker='^', markersize=4)
        ax3.set_xlabel('训练轮次')
        ax3.set_ylabel('验证AUC')
        ax3.set_title('TabTransformer 验证AUC曲线')
        ax3.grid(True, alpha=0.3)
        
        # 综合对比
        ax4.plot(epochs, training_history['val_accuracies'], 'g-', linewidth=2, label='验证准确率')
        ax4.plot(epochs, training_history['val_aucs'], 'r-', linewidth=2, label='验证AUC')
        ax4.set_xlabel('训练轮次')
        ax4.set_ylabel('指标值')
        ax4.set_title('TabTransformer 验证指标对比')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # 保存训练历史图
        history_plot_path = os.path.join(output_dir, 'transformer_training_history.png')
        plt.savefig(history_plot_path, dpi=300, bbox_inches='tight')
        plt.show()
        
        # 保存训练历史数据
        history_data = pd.DataFrame({
            'epoch': epochs,
            'train_loss': training_history['train_losses'],
            'val_accuracy': training_history['val_accuracies'],
            'val_auc': training_history['val_aucs']
        })
        
        history_data_path = os.path.join(output_dir, 'transformer_training_history.csv')
        history_data.to_csv(history_data_path, index=False, encoding='utf-8-sig')
        
        print(f"✅ 训练历史分析完成!")
        print(f"  训练历史图像: {history_plot_path}")
        print(f"  训练历史数据: {history_data_path}")
        
        return history_data
        
    except Exception as e:
        print(f"❌ 训练历史分析失败: {str(e)}")
        traceback.print_exc()
        return None

def analyze_attention_weights(model, test_loader, device, output_dir):
    """分析注意力权重（如果可能的话）"""
    print("\n🔍 分析Transformer注意力机制...")
    
    try:
        model.eval()
        
        # 获取一个批次的数据进行分析
        with torch.no_grad():
            for batch_x, batch_y in test_loader:
                batch_x = batch_x.to(device)
                
                # 由于我们的模型结构相对简单，这里主要分析模型的内部表示
                x = model.input_projection(batch_x)
                x = x.unsqueeze(1)
                x = model.pos_encoder(x)
                
                # 通过Transformer编码器
                encoded = model.transformer_encoder(x)
                
                print(f"✅ 注意力机制分析:")
                print(f"  输入维度: {batch_x.shape}")
                print(f"  投影后维度: {x.shape}")
                print(f"  编码后维度: {encoded.shape}")
                print(f"  注意力头数: 4")
                print(f"  Transformer层数: 2")
                
                break  # 只分析第一个批次
        
        return True
        
    except Exception as e:
        print(f"❌ 注意力分析失败: {str(e)}")
        traceback.print_exc()
        return False

def save_results(X_test, y_test, y_pred, names_test, feature_names, metrics, model, 
                training_history, history_data, output_dir):
    """保存模型结果"""
    print("\n💾 正在保存结果...")
    
    try:
        # 保存测试集预测结果
        test_results_df = pd.DataFrame(X_test, columns=feature_names)
        test_results_df.insert(0, 'SampleID', names_test)
        test_results_df['True_label'] = y_test
        test_results_df['Predicted_label'] = y_pred
        
        results_path = os.path.join(output_dir, 'transformer_test_results.csv')
        test_results_df.to_csv(results_path, index=False, encoding='utf-8-sig')
        
        # 保存模型
        model_path = os.path.join(output_dir, 'transformer_model.pth')
        torch.save(model.state_dict(), model_path)
        
        # 保存完整模型（用于推理）
        full_model_path = os.path.join(output_dir, 'transformer_full_model.pth')
        torch.save(model, full_model_path)
        
        # 保存评估指标
        metrics_path = os.path.join(output_dir, 'transformer_metrics.pkl')
        with open(metrics_path, 'wb') as f:
            pickle.dump(metrics, f)
        
        # 保存训练历史
        history_path = os.path.join(output_dir, 'transformer_training_history.pkl')
        with open(history_path, 'wb') as f:
            pickle.dump(training_history, f)
        
        # 保存详细报告
        report_path = os.path.join(output_dir, 'transformer_report.txt')
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("TabTransformer 模型评估报告\n")
            f.write("=" * 50 + "\n\n")
            
            f.write("模型架构信息:\n")
            f.write("-" * 20 + "\n")
            f.write("基础架构: TabTransformer (专为表格数据设计)\n")
            f.write("输入投影: Linear + BatchNorm + ReLU + Dropout\n")
            f.write("位置编码: 正弦余弦位置编码\n")
            f.write("Transformer编码器: 2层, 4个注意力头\n")
            f.write("模型维度: 64\n")
            f.write("前馈网络维度: 128\n")
            f.write("分类头: 两层全连接网络\n\n")
            
            f.write("训练配置:\n")
            f.write("-" * 20 + "\n")
            f.write("优化器: AdamW (lr=3e-4, weight_decay=1e-4)\n")
            f.write("学习率调度: ReduceLROnPlateau\n")
            f.write("梯度裁剪: max_norm=1.0\n")
            f.write("早停耐心: 5轮\n")
            f.write("批次大小: 32\n")
            f.write(f"实际训练轮数: {training_history['epochs_trained']}\n")
            f.write(f"最佳验证AUC: {training_history['best_val_auc']:.4f}\n\n")
            
            f.write("测试集性能指标:\n")
            f.write("-" * 20 + "\n")
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
        print(f"  训练历史: {history_path}")
        print(f"  详细报告: {report_path}")
        
    except Exception as e:
        print(f"❌ 保存结果失败: {str(e)}")
        traceback.print_exc()

def main():
    """主函数"""
    print("\n" + "="*60)
    print("--- 🚀 TabTransformer 模型训练与评估 ---")
    print("         (使用100个ANOVA基因)")
    print("="*60)
    
    # 检查CUDA可用性
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"使用设备: {device}")
    if torch.cuda.is_available():
        print(f"GPU型号: {torch.cuda.get_device_name(0)}")
    
    # 设置输出目录
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    output_dir = os.path.join(project_root, "results")
    figure_dir = os.path.join(project_root, "figure")
    os.makedirs(figure_dir, exist_ok=True)
    
    # 1. 加载预处理数据
    train_data, val_data, test_data, metadata = load_preprocessed_data()
    
    feature_names = metadata['feature_names']
    names_test = test_data['names_test']
    
    # 准备类别名称
    unique_labels = np.unique(train_data['y_train'])
    class_names = [f'Class {i}' for i in unique_labels]
    
    # 2. 准备数据加载器
    train_loader, val_loader, test_loader, input_dim = prepare_data_loaders(
        train_data, val_data, test_data, batch_size=32
    )
    
    if train_loader is None:
        print("❌ 数据加载器准备失败，程序退出")
        return
    
    # 3. 创建模型
    model = create_transformer_model(input_dim, num_classes=4)
    if model is None:
        print("❌ 模型创建失败，程序退出")
        return
    
    model = model.to(device)
    
    # 计算模型参数
    total_params = sum(p.numel() for p in model.parameters())
    trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"模型参数数量: {total_params:,}")
    print(f"可训练参数数量: {trainable_params:,}")
    
    # 4. 训练模型
    model, training_history = train_model(
        model, train_loader, val_loader, device, epochs=30, patience=5
    )
    
    if model is None:
        print("❌ 模型训练失败，程序退出")
        return
    
    # 5. 评估模型
    metrics, y_pred, y_proba = evaluate_model(
        model, test_loader, device, feature_names, class_names, figure_dir
    )
    
    if metrics is None:
        print("❌ 模型评估失败，程序退出")
        return
    
    # 6. 分析注意力机制
    analyze_attention_weights(model, test_loader, device, figure_dir)
    
    # 7. 绘制训练历史
    history_data = plot_training_history(training_history, figure_dir)
    
    # 8. 保存结果
    save_results(
        test_data['X_test'], test_data['y_test'], y_pred, names_test,
        feature_names, metrics, model, training_history, history_data, output_dir
    )
    
    print(f"\n🎉 TabTransformer模型训练和评估完成!")
    print(f"   最终测试准确率: {metrics['Accuracy']:.4f}")
    print(f"   最终测试AUC: {metrics['AUC']:.4f}")

if __name__ == "__main__":
    main()
