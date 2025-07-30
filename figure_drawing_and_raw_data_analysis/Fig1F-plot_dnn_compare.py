import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

folder = os.path.dirname(os.path.abspath(__file__))
input_path = os.path.join(folder, 'heterodimer_batch123_sn_diff.csv')
output_pdf = os.path.join(folder, 'DNN_compare_hist.pdf')
output_bins_csv = os.path.join(folder, 'DNN_compare_hist_bins.csv')

# 读取数据
df = pd.read_csv(input_path)

# 两组筛选
cond = (df['LIS'] >= 0.203) & (df['LIA'] >= 3432)
group1 = df[cond]['DNN']
group2 = df[~cond]['DNN']

# 只统计1到10的DNN
bins = np.arange(1, 12) - 0.5  # 1到10的整数bin
labels = np.arange(1, 11)

hist1, _ = np.histogram(group1, bins=bins)
hist2, _ = np.histogram(group2, bins=bins)

x = np.arange(1, 11)
x_csv = x - 1  # 横坐标减去1，变为0到9
bar_width = 0.4  # 设置为0.4

plt.figure(figsize=(7, 6))  # 图宽减小
# 论文绿 #2ca02c 用于满足条件，论文蓝 #1f77b4 用于Other
bars1 = plt.bar(x - bar_width/2, hist1, width=bar_width, color='#2ca02c', label='LIS≥0.203 & LIA≥3432', zorder=3, linewidth=1.2, edgecolor='black')
bars2 = plt.bar(x + bar_width/2, hist2, width=bar_width, color='#1f77b4', label='Other', zorder=2, alpha=0.7, linewidth=1.2, edgecolor='black')

# 所有柱子描边统一粗细
# （已在bar绘制时设置linewidth=1.2）

plt.xlabel('DNN', fontname='Arial', fontsize=14)
plt.ylabel('Count', fontname='Arial', fontsize=14)
plt.title('DNN Distribution Comparison', fontname='Arial', fontsize=16)
plt.xticks(x)
plt.legend()
plt.grid(False)
plt.gca().set_facecolor('white')
plt.gcf().patch.set_facecolor('white')
plt.tight_layout()
plt.savefig(output_pdf, bbox_inches='tight', dpi=300, transparent=True)
plt.close()
print(f'已输出: {output_pdf}')

# 保存柱子数值到csv
bins_df = pd.DataFrame({'DNN': x_csv, 'LIS>=0.203_LIA>=3432': hist1, 'Other': hist2})
bins_df.to_csv(output_bins_csv, index=False)
print(f'已输出: {output_bins_csv}') 