import warnings
import pickle
import os
import sys
import shap
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shap import kmeans
from IPython.display import display

warnings.filterwarnings('ignore')


# ==================== 新增：创建 SHAP 输出目录 ====================
# 直接保存在当前脚本目录（即4-SHAP/）
script_dir = os.path.dirname(os.path.abspath(__file__))
shap_output_dir = script_dir
print(f"📁 SHAP 输出目录: {shap_output_dir}")
# ===============================================================

# ==================== 保存图像函数（保存到当前目录） ====================
def save_figure(fig, filename, dpi=300):
    """保存图像为 PNG 和 PDF 到当前目录（4-SHAP/）"""
    png_path = os.path.join(shap_output_dir, f"{filename}.png")
    pdf_path = os.path.join(shap_output_dir, f"{filename}.pdf")
    fig.savefig(png_path, dpi=dpi, bbox_inches='tight')
    fig.savefig(pdf_path, dpi=dpi, bbox_inches='tight')
    print(f"✅ 已保存图像: {filename}{{.png, .pdf}}")
# ===================================================================

# ==================== 加载模型和数据（从 results/） ====================
def load_model_and_data():
    project_root = os.path.dirname(script_dir)
    results_dir = os.path.join(project_root, "results")
    
    print("正在加载模型和数据...")
    try:
        model_path = os.path.join(results_dir, 'gradient_boosting_model.pkl')
        with open(model_path, 'rb') as f:
            model = pickle.load(f)
        print(f"✅ 模型加载成功: {model_path}")

        with open(os.path.join(results_dir, 'train_data.pkl'), 'rb') as f:
            train_data = pickle.load(f)
        X_train_scaled = train_data['X_train']

        with open(os.path.join(results_dir, 'val_data.pkl'), 'rb') as f:
            val_data = pickle.load(f)
        X_val_scaled = val_data['X_val']
        y_val = val_data['y_val']
        
        with open(os.path.join(results_dir, 'metadata.pkl'), 'rb') as f:
            metadata = pickle.load(f)
        feature_names = metadata['feature_names']
        
        print("✅ 数据加载成功!")
        print(f"训练集形状: {X_train_scaled.shape}")
        print(f"验证集形状: {X_val_scaled.shape}")
        print(f"特征数量: {len(feature_names)}")
        
        return model, X_train_scaled, X_val_scaled, y_val, feature_names
        
    except FileNotFoundError as e:
        print(f"❌ 错误: 找不到文件: {e}")
        print("请确保已运行 '05_data_preprocessing.py' 和 '06_gradient_boosting_model.py'")
        sys.exit(1)
    except Exception as e:
        print(f"❌ 加载时出错: {e}")
        sys.exit(1)
# ===================================================================

# ==================== 主程序 ====================
print("\n  Gradient Boosting SHAP 可视化分析")
print("=" * 60)

gb_model, X_train_scaled, X_val_scaled, y_val, feature_names = load_model_and_data()

print("  创建KernelExplainer的背景数据集...")
background_data = kmeans(X_train_scaled, 50)

print("  初始化SHAP KernelExplainer (这可能需要一些时间)...")
explainer = shap.KernelExplainer(gb_model.predict_proba, background_data)

sample_size = min(100, len(X_val_scaled))
X_val_sample = pd.DataFrame(X_val_scaled[:sample_size], columns=feature_names)
y_val_sample = y_val.iloc[:sample_size]

print(f"  计算SHAP值 (使用{sample_size}个样本)...")
shap_values = explainer.shap_values(X_val_sample)

print(f"✅ SHAP值计算完成，为 {len(shap_values)} 个类别生成了SHAP值。")

# ==================== 保存 SHAP 结果为 .pkl ====================
print("\n💾 正在保存 SHAP 结果到 'shap_results.pkl'...")
shap_results = {
    'shap_values': shap_values,
    'X_val_sample': X_val_sample,
    'y_val_sample': y_val_sample,
    'feature_names': feature_names,
    'expected_values': explainer.expected_value,
    # 注意：explainer 本身可能无法可靠 pickle，故不保存
}
shap_pkl_path = os.path.join(shap_output_dir, "shap_results.pkl")
with open(shap_pkl_path, 'wb') as f:
    pickle.dump(shap_results, f)
print(f"✅ SHAP 结果已保存至: {shap_pkl_path}")
# ===============================================================

print("  开始生成可视化图表...")

# ==================== 1. SHAP Summary Plot (Bar) ====================
print("\n  1. 生成SHAP Summary Plot (Bar)...")
plt.figure()
shap.summary_plot(shap_values, X_val_sample, plot_type="bar", show=False)
plt.title('Gradient Boosting SHAP Summary Plot (All Classes)', fontsize=16, fontweight='bold', pad=20)
plt.tight_layout()
save_figure(plt.gcf(), "shap_summary_bar")
plt.show()

# ==================== 2. SHAP Beeswarm Plot (All Classes) ====================
print("\n  2. 生成 SHAP Beeswarm Plot (All Classes)...")
plt.figure()
shap.summary_plot(shap_values, X_val_sample, show=False)
plt.title('Gradient Boosting SHAP Beeswarm Plot (All Classes)', fontsize=16, fontweight='bold', pad=20)
plt.tight_layout()
save_figure(plt.gcf(), "shap_beeswarm_all")
plt.show()

# ==================== 3. 按类别分析 Beeswarm Plot ====================
print("\n  3. 按类别生成 SHAP Beeswarm Plots...")
class_names = [f'Class {i}' for i in range(len(shap_values))]
for i, class_name in enumerate(class_names):
    plt.figure()
    shap.summary_plot(shap_values[i], X_val_sample, show=False)
    plt.title(f'SHAP Beeswarm Plot for {class_name}', fontsize=16, fontweight='bold', pad=20)
    plt.tight_layout()
    save_figure(plt.gcf(), f"shap_beeswarm_{class_name.replace(' ', '_')}")
    plt.show()

# ==================== 4. SHAP Waterfall Plot ====================
print("\n  4. 生成 SHAP Waterfall Plot...")
sample_idx = 0
for i, class_name in enumerate(class_names):
    shap_explanation_for_class = shap.Explanation(
        values=shap_values[i][sample_idx, :],
        base_values=explainer.expected_value[i],
        data=X_val_sample.iloc[sample_idx],
        feature_names=feature_names
    )
    plt.figure()
    shap.plots.waterfall(shap_explanation_for_class, max_display=20, show=False)
    plt.title(f'Waterfall for Sample {sample_idx} - Prediction for {class_name}', fontweight='bold')
    plt.tight_layout()
    save_figure(plt.gcf(), f"shap_waterfall_sample{sample_idx}_class{i}")
    plt.show()

# ==================== 5. SHAP Dependence Plot ====================
print("\n  5. 生成SHAP Dependence Plots...")
global_shap_importance = np.abs(np.array(shap_values)).mean(axis=(0, 1))
top_features_indices = np.argsort(-global_shap_importance)[:4]
class_to_plot = 1 
print(f"  (为 Class {class_to_plot} 绘制依赖图)")

fig, axes = plt.subplots(2, 2, figsize=(15, 12))
axes = axes.flatten()

for i, feature_idx in enumerate(top_features_indices):
    feature_name = feature_names[feature_idx]
    plt.sca(axes[i])
    shap.dependence_plot(
        feature_name,
        shap_values[class_to_plot],
        X_val_sample,
        ax=axes[i],
        show=False
    )
    axes[i].set_title(f'SHAP Dependence: {feature_name} (for Class {class_to_plot})', fontweight='bold', fontsize=12)

plt.tight_layout()
save_figure(fig, "shap_dependence_top4")
plt.show()

# ==================== 6. SHAP Force Plot ====================
print("\n  6. 生成SHAP Force Plot...")
shap.initjs()
for i, class_name in enumerate(class_names):
    print(f"  Force plot for sample 0, class {class_name}")
    display(shap.force_plot(explainer.expected_value[i], shap_values[i][0,:], X_val_sample.iloc[0,:], feature_names))

# ==================== 7. SHAP Heatmap ====================
print("\n  7. 生成SHAP Heatmap...")
plt.figure(figsize=(15, 8))
shap_explanation_heatmap = shap.Explanation(
    values=shap_values[class_to_plot],
    base_values=explainer.expected_value[class_to_plot],
    data=X_val_sample,
    feature_names=feature_names
)
shap.plots.heatmap(shap_explanation_heatmap, max_display=20, show=False)
plt.title(f'SHAP Heatmap for Class {class_to_plot}', fontsize=16, fontweight='bold')
plt.tight_layout()
save_figure(plt.gcf(), "shap_heatmap_class1")
plt.show()

print("\n" + "=" * 60)
print("  Gradient Boosting SHAP 可视化分析完成!")
print("=" * 60)