#!/usr/bin/env python3
"""
绘制AF3-Homo与其他9种蛋白质预测方法的lDDT分析图
基于plot_lddt_comparison.py的风格，包含柱状图和箱线图，使用英文标签，保存为PDF格式
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# 设置字体和样式
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['axes.edgecolor'] = 'black'

def load_and_process_data():
    """加载并处理AF3-Homo和其他方法的数据"""
    print("Loading AF3-Homo and other methods data...")
    
    # 加载AF3-Homo结果
    af3_homo_file = "/Users/david/Desktop/h5o/af3_homo_analysis_20250627_180331/analysis_results.csv"
    af3_homo_df = pd.read_csv(af3_homo_file)
    
    # 加载其他方法结果
    other_methods_file = "/Users/david/Desktop/h5o/comprehensive_homo_analysis_fixed_20250627_134353/analysis_results.csv"
    other_methods_df = pd.read_csv(other_methods_file)
    
    print(f"AF3-Homo raw data: {len(af3_homo_df)} records")
    print(f"Other methods raw data: {len(other_methods_df)} records")
    
    # 过滤成功的结果和有效的分数
    af3_homo_df = af3_homo_df[
        (af3_homo_df['status'] == 'success') & 
        (af3_homo_df['lddt_score'] >= 0) & 
        (af3_homo_df['dockq_score'] >= 0)
    ].copy()
    
    other_methods_df = other_methods_df[
        (other_methods_df['status'] == 'success') & 
        (other_methods_df['lddt_score'] >= 0) & 
        (other_methods_df['dockq_score'] >= 0)
    ].copy()
    
    # 重命名AF3-Homo的方法名
    af3_homo_df['method_new'] = 'AlphaFold3'
    
    # 处理ColabFold分类，基于文件名
    def classify_method(row):
        method = row['method']
        filename = row['filename']
        
        if method == 'colabfold':
            if 'unrelaxed' in filename.lower():
                return 'ColabFold_Unrelaxed'
            elif 'relaxed' in filename.lower():
                return 'ColabFold_Relaxed'
            else:
                return 'ColabFold_Relaxed'
        else:
            return method
    
    # 应用分类函数到其他方法
    other_methods_df['method_new'] = other_methods_df.apply(classify_method, axis=1)
    
    # 合并数据
    combined_df = pd.concat([af3_homo_df, other_methods_df], ignore_index=True)
    
    print(f"AF3-Homo valid data: {len(af3_homo_df)} records")
    print(f"Other methods valid data: {len(other_methods_df)} records")
    print(f"Combined valid data: {len(combined_df)} records")
    
    return combined_df

def create_lddt_plots(df, output_dir):
    """创建lDDT分析图"""
    
    # 方法名映射（用于显示）
    method_mapping = {
        'AlphaFold3': 'AlphaFold3',
        'boltz2': 'Boltz',
        'chai1_esm': 'Chai1-ESM',
        'chai1_msa': 'Chai1-MSA',
        'ColabFold_Relaxed': 'ColabFold-Relaxed',
        'ColabFold_Unrelaxed': 'ColabFold-Unrelaxed',
        'protenix_3cycle_esm': 'Protenix-3cycle-ESM',
        'protenix_3cycle_msa': 'Protenix-3cycle-MSA',
        'protenix_esm': 'Protenix-ESM',
        'protenix_msa': 'Protenix-MSA'
    }
    
    # 按照指定顺序排列方法
    method_order = [
        'AlphaFold3',
        'boltz2', 
        'chai1_esm',
        'chai1_msa',
        'ColabFold_Unrelaxed',
        'ColabFold_Relaxed',
        'protenix_3cycle_esm',
        'protenix_3cycle_msa',
        'protenix_esm',
        'protenix_msa'
    ]
    
    # 获取实际存在的方法
    available_methods = set(df['method_new'].unique())
    methods = [method for method in method_order if method in available_methods]
    print(f"Found {len(methods)} methods: {methods}")
    
    # 计算每种方法的统计数据
    method_stats = []
    for method in methods:
        method_data = df[df['method_new'] == method]['lddt_score']
        stats = {
            'method': method,
            'display_name': method_mapping.get(method, method),
            'mean': method_data.mean(),
            'std': method_data.std(),
            'median': method_data.median(),
            'count': len(method_data)
        }
        method_stats.append(stats)
    
    stats_df = pd.DataFrame(method_stats)
    
    # 颜色方案 - 使用专业的色彩搭配
    colors = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D', '#592E83', 
              '#1A5490', '#E07A5F', '#3D5A80', '#98C1D9', '#EE6C4D']
    
    # 创建图表
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # 1. 柱状图 - 显示平均lDDT分数
    x_pos = np.arange(len(stats_df))
    bars = ax1.bar(x_pos, stats_df['mean'], 
                   color=[colors[i % len(colors)] for i in range(len(stats_df))], 
                   alpha=0.8, edgecolor='black', linewidth=0.5)
    
    # 添加误差棒
    ax1.errorbar(x_pos, stats_df['mean'], yerr=stats_df['std'], 
                fmt='none', color='black', capsize=3, capthick=1)
    
    # 添加数值标签
    for i, (bar, mean_val, std_val) in enumerate(zip(bars, stats_df['mean'], stats_df['std'])):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + std_val + 0.01,
                f'{mean_val:.3f}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    ax1.set_xlabel('Prediction Methods', fontweight='bold')
    ax1.set_ylabel('Mean lDDT Score', fontweight='bold')
    ax1.set_title('Average lDDT Scores Across Ten Methods', fontweight='bold', fontsize=12)
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(stats_df['display_name'], rotation=45, ha='right', fontsize=9)
    ax1.set_ylim(0, 1)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # 2. 箱线图
    lddt_data = [df[df['method_new'] == method]['lddt_score'].values for method in methods]
    
    bp = ax2.boxplot(lddt_data, labels=stats_df['display_name'], patch_artist=True,
                     boxprops=dict(linewidth=1.2),
                     whiskerprops=dict(linewidth=1.2),
                     capprops=dict(linewidth=1.2),
                     medianprops=dict(linewidth=1.5, color='black'))
    
    # 为每个箱子设置颜色
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)
    
    ax2.set_xlabel('Prediction Methods', fontweight='bold')
    ax2.set_ylabel('lDDT Score', fontweight='bold')
    ax2.set_title('lDDT Score Distribution by Method', fontweight='bold', fontsize=12)
    ax2.tick_params(axis='x', rotation=45, labelsize=9)
    ax2.set_ylim(0, 1)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    plt.tight_layout()
    
    # 保存为PDF
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = os.path.join(output_dir, f'AF3_Homo_lDDT_analysis_{timestamp}.pdf')
    plt.savefig(output_file, format='pdf', dpi=300, bbox_inches='tight')
    print(f"AF3-Homo lDDT analysis plot saved to: {output_file}")
    
    # 也保存为PNG以便查看
    png_file = os.path.join(output_dir, f'AF3_Homo_lDDT_analysis_{timestamp}.png')
    plt.savefig(png_file, format='png', dpi=300, bbox_inches='tight')
    print(f"AF3-Homo lDDT analysis plot also saved as PNG: {png_file}")
    
    plt.show()
    
    # 输出详细统计
    print("\n=== lDDT Statistics ===")
    for _, row in stats_df.iterrows():
        print(f"{row['display_name']:20s}: Mean={row['mean']:.3f}, "
              f"Std={row['std']:.3f}, "
              f"Median={row['median']:.3f}, "
              f"Count={row['count']}")
    
    return stats_df

def main():
    """主函数"""
    # 创建输出目录
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = f"af3_homo_lddt_analysis_{timestamp}"
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"AF3-Homo lDDT Analysis for Ten Protein Prediction Methods")
    print(f"Output directory: {output_dir}")
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # 加载和处理数据
    df = load_and_process_data()
    
    # 创建lDDT分析图
    stats_df = create_lddt_plots(df, output_dir)
    
    # 保存统计数据
    stats_file = os.path.join(output_dir, f'AF3_Homo_lDDT_statistics_{timestamp}.csv')
    stats_df.to_csv(stats_file, index=False)
    print(f"Statistics saved to: {stats_file}")
    
    print(f"\nAnalysis completed! Results saved in: {output_dir}")

if __name__ == "__main__":
    main() 