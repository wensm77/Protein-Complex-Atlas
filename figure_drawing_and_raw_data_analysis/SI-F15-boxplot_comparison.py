import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# 设置matplotlib字体
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Helvetica']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = 16  # 增大默认字体大小

# 文件路径
file1 = '/Users/david/Desktop/0531/fig4-low/merged_data.csv'
file2 = '/Users/david/Desktop/0531/fig4-low/115w_scores_plus_LIS_0409_human_fasta6_both_matched.csv'

print("Reading CSV files...")

# 读取两个CSV文件
print("Reading first file: merged_data.csv")
df1 = pd.read_csv(file1)
print(f"First file: {len(df1)} rows, {len(df1.columns)} columns")

print("Reading second file: 115w_scores_plus_LIS_0409_human_fasta6_both_matched.csv")
df2 = pd.read_csv(file2)
print(f"Second file: {len(df2)} rows, {len(df2.columns)} columns")

# 计算len1+len2
print("\nCalculating total length...")
df1['total_length'] = df1['len1'] + df1['len2']
df2['total_length'] = df2['len1'] + df2['len2']

# 准备绘图数据
print("Preparing plot data...")
# 为plddt_x创建数据
plddt_data = []
plddt_data.extend([(val, 'merged_data') for val in df1['plddt_x'].dropna()])
plddt_data.extend([(val, 'human_fasta6') for val in df2['plddt_x'].dropna()])
plddt_df = pd.DataFrame(plddt_data, columns=['plddt_x', 'dataset'])

# 为total_length创建数据
length_data = []
length_data.extend([(val, 'merged_data') for val in df1['total_length'].dropna()])
length_data.extend([(val, 'human_fasta6') for val in df2['total_length'].dropna()])
length_df = pd.DataFrame(length_data, columns=['total_length', 'dataset'])

print(f"plddt_x data points: {len(plddt_df)}")
print(f"total_length data points: {len(length_df)}")

# 计算统计数据
print("Calculating statistics...")
stats_data = []

# plddt_x statistics
for dataset_name, df in [('merged_data', df1), ('human_fasta6', df2)]:
    plddt_stats = df['plddt_x'].describe(percentiles=[0.25, 0.5, 0.75])
    length_stats = df['total_length'].describe(percentiles=[0.25, 0.5, 0.75])
    
    stats_data.append({
        'Dataset': dataset_name,
        'Metric': 'plddt_x',
        'Count': plddt_stats['count'],
        'Mean': plddt_stats['mean'],
        'Std': plddt_stats['std'],
        'Min': plddt_stats['min'],
        'Q1': plddt_stats['25%'],
        'Median': plddt_stats['50%'],
        'Q3': plddt_stats['75%'],
        'Max': plddt_stats['max']
    })
    
    stats_data.append({
        'Dataset': dataset_name,
        'Metric': 'total_length',
        'Count': length_stats['count'],
        'Mean': length_stats['mean'],
        'Std': length_stats['std'],
        'Min': length_stats['min'],
        'Q1': length_stats['25%'],
        'Median': length_stats['50%'],
        'Q3': length_stats['75%'],
        'Max': length_stats['max']
    })

# 创建统计数据DataFrame并保存到CSV
stats_df = pd.DataFrame(stats_data)
stats_csv_file = '/Users/david/Desktop/0531/fig4-low/boxplot_statistics.csv'
stats_df.to_csv(stats_csv_file, index=False)
print(f"Statistics saved to: {stats_csv_file}")

# 设置图形样式
sns.set_style("white")
plt.style.use('default')

# 创建子图，调整为更窄的尺寸
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 8))

# 绘制plddt_x的箱线图
print("Creating plddt_x boxplot...")
box1 = sns.boxplot(data=plddt_df, x='dataset', y='plddt_x', ax=ax1, 
                   palette=['#1f77b4', '#ff7f0e'], 
                   showfliers=False,  # 不显示异常值
                   linewidth=2.5)
ax1.set_title('plddt_x Distribution Comparison', fontsize=24, fontweight='bold', fontfamily='Arial')
ax1.set_xlabel('Dataset', fontsize=20, fontfamily='Arial')
ax1.set_ylabel('plddt_x', fontsize=20, fontfamily='Arial')
ax1.tick_params(axis='both', labelsize=16)  # 刻度标签字体大小
ax1.grid(False)
ax1.set_facecolor('white')

# 绘制total_length的箱线图
print("Creating total_length boxplot...")
box2 = sns.boxplot(data=length_df, x='dataset', y='total_length', ax=ax2,
                   palette=['#1f77b4', '#ff7f0e'],
                   showfliers=False,  # 不显示异常值
                   linewidth=2.5)
ax2.set_title('Total Length Distribution Comparison', fontsize=24, fontweight='bold', fontfamily='Arial')
ax2.set_xlabel('Dataset', fontsize=20, fontfamily='Arial')
ax2.set_ylabel('Total Length (len1+len2)', fontsize=20, fontfamily='Arial')
ax2.set_ylim(0, 2000)  # 设置y轴范围为0-2000
ax2.tick_params(axis='both', labelsize=16)  # 刻度标签字体大小
ax2.grid(False)
ax2.set_facecolor('white')

# 设置图形背景为白色，去除灰色背景
fig.patch.set_facecolor('white')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

# 调整布局
plt.tight_layout()

# 保存为PDF
output_file = '/Users/david/Desktop/0531/fig4-low/boxplot_comparison.pdf'
print(f"\nSaving plot to: {output_file}")
plt.savefig(output_file, format='pdf', dpi=300, bbox_inches='tight', 
           facecolor='white', edgecolor='none')

print("Plot saved successfully!")

# 显示统计信息
print(f"\nStatistics Summary:")
print(f"merged_data.csv:")
print(f"  - plddt_x: Mean={df1['plddt_x'].mean():.2f}, Median={df1['plddt_x'].median():.2f}")
print(f"  - total_length: Mean={df1['total_length'].mean():.2f}, Median={df1['total_length'].median():.2f}")

print(f"\nhuman_fasta6.csv:")
print(f"  - plddt_x: Mean={df2['plddt_x'].mean():.2f}, Median={df2['plddt_x'].median():.2f}")
print(f"  - total_length: Mean={df2['total_length'].mean():.2f}, Median={df2['total_length'].median():.2f}")

# 显示四分位数信息
print(f"\nBoxplot Statistics:")
print(f"\nplddt_x - merged_data:")
q1, median, q3 = df1['plddt_x'].quantile([0.25, 0.5, 0.75])
print(f"  Q1: {q1:.2f}, Median: {median:.2f}, Q3: {q3:.2f}")

print(f"\nplddt_x - human_fasta6:")
q1, median, q3 = df2['plddt_x'].quantile([0.25, 0.5, 0.75])
print(f"  Q1: {q1:.2f}, Median: {median:.2f}, Q3: {q3:.2f}")

print(f"\ntotal_length - merged_data:")
q1, median, q3 = df1['total_length'].quantile([0.25, 0.5, 0.75])
print(f"  Q1: {q1:.2f}, Median: {median:.2f}, Q3: {q3:.2f}")

print(f"\ntotal_length - human_fasta6:")
q1, median, q3 = df2['total_length'].quantile([0.25, 0.5, 0.75])
print(f"  Q1: {q1:.2f}, Median: {median:.2f}, Q3: {q3:.2f}")

print("\nBoxplot generation completed!")
print(f"Files generated:")
print(f"  - Plot: {output_file}")
print(f"  - Statistics: {stats_csv_file}") 