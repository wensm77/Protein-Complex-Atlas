#!/usr/bin/env python3
"""
绘制AF3-Homo与其他9种蛋白质预测方法的DockQ分析图
基于plot_dockq_analysis_副本.py的风格，包含均值柱状图和质量等级分类柱状图，使用英文标签，保存为PDF格式
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

def classify_dockq_quality(score):
    """根据DockQ分数分类质量等级"""
    if score <= 0.23:
        return 'Incorrect'
    elif score < 0.49:
        return 'Acceptable'
    elif score < 0.8:
        return 'Medium'
    else:
        return 'High'

def create_dockq_plots(df, output_dir):
    """创建DockQ分析图"""
    
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
    
    # 按照参考图片的顺序排列方法
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
    
    # 为每个预测添加质量分类
    df['dockq_quality'] = df['dockq_score'].apply(classify_dockq_quality)
    
    # 计算每种方法的统计数据
    method_stats = []
    quality_stats = []
    
    for method in methods:
        method_data = df[df['method_new'] == method]
        dockq_scores = method_data['dockq_score']
        
        # 基本统计
        stats = {
            'method': method,
            'display_name': method_mapping.get(method, method),
            'mean': dockq_scores.mean(),
            'std': dockq_scores.std(),
            'median': dockq_scores.median(),
            'count': len(dockq_scores)
        }
        method_stats.append(stats)
        
        # 质量等级统计
        total_count = len(method_data)
        quality_counts = method_data['dockq_quality'].value_counts()
        
        quality_stat = {
            'method': method,
            'display_name': method_mapping.get(method, method),
            'total': total_count,
            'incorrect_count': quality_counts.get('Incorrect', 0),
            'acceptable_count': quality_counts.get('Acceptable', 0),
            'medium_count': quality_counts.get('Medium', 0),
            'high_count': quality_counts.get('High', 0),
            'incorrect_pct': quality_counts.get('Incorrect', 0) / total_count * 100,
            'acceptable_pct': quality_counts.get('Acceptable', 0) / total_count * 100,
            'medium_pct': quality_counts.get('Medium', 0) / total_count * 100,
            'high_pct': quality_counts.get('High', 0) / total_count * 100
        }
        quality_stats.append(quality_stat)
    
    stats_df = pd.DataFrame(method_stats)
    quality_df = pd.DataFrame(quality_stats)
    
    # 颜色方案
    colors = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D', '#592E83', 
              '#1A5490', '#E07A5F', '#3D5A80', '#98C1D9', '#EE6C4D']
    
    # 质量等级颜色
    quality_colors = {
        'High': '#2E8B57',      # 深绿色
        'Medium': '#FFD700',    # 金色
        'Acceptable': '#FF8C00', # 橙色
        'Incorrect': '#DC143C'   # 深红色
    }
    
    # 创建图表
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7))
    
    # 1. 均值柱状图
    x_pos = np.arange(len(stats_df))
    bars = ax1.bar(x_pos, stats_df['mean'], 
                   color=[colors[i % len(colors)] for i in range(len(stats_df))], 
                   alpha=0.8, edgecolor='black', linewidth=0.5)
    
    # 添加数值标签
    for i, (bar, mean_val) in enumerate(zip(bars, stats_df['mean'])):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 0.02,
                f'{mean_val:.3f}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    ax1.set_xlabel('Prediction Methods', fontweight='bold')
    ax1.set_ylabel('Mean DockQ Score', fontweight='bold')
    ax1.set_title('Average DockQ Scores Across Ten Methods', fontweight='bold', fontsize=12)
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(stats_df['display_name'], rotation=45, ha='right', fontsize=9)
    ax1.set_ylim(0, 0.8)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # 2. 质量等级分类堆叠柱状图
    x_pos2 = np.arange(len(quality_df))
    width = 0.8
    
    # 堆叠柱状图
    bottom_incorrect = np.zeros(len(quality_df))
    bottom_acceptable = quality_df['incorrect_pct'].values
    bottom_medium = bottom_acceptable + quality_df['acceptable_pct'].values
    bottom_high = bottom_medium + quality_df['medium_pct'].values
    
    bars_incorrect = ax2.bar(x_pos2, quality_df['incorrect_pct'], width, 
                            label='Incorrect (≤0.23)', color=quality_colors['Incorrect'], alpha=0.8)
    bars_acceptable = ax2.bar(x_pos2, quality_df['acceptable_pct'], width, 
                             bottom=bottom_acceptable, label='Acceptable (0.23-0.49)', 
                             color=quality_colors['Acceptable'], alpha=0.8)
    bars_medium = ax2.bar(x_pos2, quality_df['medium_pct'], width, 
                         bottom=bottom_medium, label='Medium (0.49-0.8)', 
                         color=quality_colors['Medium'], alpha=0.8)
    bars_high = ax2.bar(x_pos2, quality_df['high_pct'], width, 
                       bottom=bottom_high, label='High (≥0.8)', 
                       color=quality_colors['High'], alpha=0.8)
    
    # 添加百分比标签
    for i, row in quality_df.iterrows():
        # 只显示比例较大的标签，避免拥挤
        if row['high_pct'] > 5:
            ax2.text(i, bottom_high[i] + row['high_pct']/2, f"{row['high_pct']:.1f}%", 
                    ha='center', va='center', fontsize=8, fontweight='bold')
        if row['medium_pct'] > 5:
            ax2.text(i, bottom_medium[i] + row['medium_pct']/2, f"{row['medium_pct']:.1f}%", 
                    ha='center', va='center', fontsize=8, fontweight='bold')
        if row['acceptable_pct'] > 5:
            ax2.text(i, bottom_acceptable[i] + row['acceptable_pct']/2, f"{row['acceptable_pct']:.1f}%", 
                    ha='center', va='center', fontsize=8, fontweight='bold')
        if row['incorrect_pct'] > 5:
            ax2.text(i, row['incorrect_pct']/2, f"{row['incorrect_pct']:.1f}%", 
                    ha='center', va='center', fontsize=8, fontweight='bold')
    
    ax2.set_xlabel('Prediction Methods', fontweight='bold')
    ax2.set_ylabel('Percentage (%)', fontweight='bold')
    ax2.set_title('DockQ Quality Distribution by Method', fontweight='bold', fontsize=12)
    ax2.set_xticks(x_pos2)
    ax2.set_xticklabels(quality_df['display_name'], rotation=45, ha='right', fontsize=9)
    ax2.set_ylim(0, 100)
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    plt.tight_layout()
    
    # 保存为PDF
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = os.path.join(output_dir, f'AF3_Homo_DockQ_analysis_{timestamp}.pdf')
    plt.savefig(output_file, format='pdf', dpi=300, bbox_inches='tight')
    print(f"AF3-Homo DockQ analysis plot saved to: {output_file}")
    
    # 也保存为PNG以便查看
    png_file = os.path.join(output_dir, f'AF3_Homo_DockQ_analysis_{timestamp}.png')
    plt.savefig(png_file, format='png', dpi=300, bbox_inches='tight')
    print(f"AF3-Homo DockQ analysis plot also saved as PNG: {png_file}")
    
    plt.show()
    
    # 输出详细统计
    print("\n=== DockQ Statistics ===")
    for _, row in stats_df.iterrows():
        print(f"{row['display_name']:20s}: Mean={row['mean']:.3f}, "
              f"Std={row['std']:.3f}, "
              f"Median={row['median']:.3f}, "
              f"Count={row['count']}")
    
    print("\n=== DockQ Quality Distribution ===")
    for _, row in quality_df.iterrows():
        print(f"{row['display_name']:20s}: "
              f"High={row['high_pct']:.1f}%, "
              f"Medium={row['medium_pct']:.1f}%, "
              f"Acceptable={row['acceptable_pct']:.1f}%, "
              f"Incorrect={row['incorrect_pct']:.1f}%")
    
    return stats_df, quality_df

def main():
    """主函数"""
    # 创建输出目录
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = f"af3_homo_dockq_analysis_{timestamp}"
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"AF3-Homo DockQ Analysis for Ten Protein Prediction Methods")
    print(f"Output directory: {output_dir}")
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # 加载和处理数据
    df = load_and_process_data()
    
    # 创建DockQ分析图
    stats_df, quality_df = create_dockq_plots(df, output_dir)
    
    # 保存统计数据
    stats_file = os.path.join(output_dir, f'AF3_Homo_DockQ_statistics_{timestamp}.csv')
    stats_df.to_csv(stats_file, index=False)
    
    quality_file = os.path.join(output_dir, f'AF3_Homo_DockQ_quality_distribution_{timestamp}.csv')
    quality_df.to_csv(quality_file, index=False)
    
    print(f"Statistics saved to: {stats_file}")
    print(f"Quality distribution saved to: {quality_file}")
    
    print(f"\nAnalysis completed! Results saved in: {output_dir}")

if __name__ == "__main__":
    main() 