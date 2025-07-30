import pandas as pd
import matplotlib.pyplot as plt
import os

folder = os.path.dirname(os.path.abspath(__file__))
input_path = os.path.join(folder, 'homodimer_seq1_matched.csv')
output_pdf = os.path.join(folder, 'homodimer_LIS_hist.pdf')
output_bins_csv = os.path.join(folder, 'homodimer_LIS_hist_bins.csv')

# 读取数据
df = pd.read_csv(input_path, usecols=['LIS'])

# 去除缺失值
lis = df['LIS'].dropna()

plt.figure(figsize=(8, 6))
bins = 15
counts, bins_edges, _ = plt.hist(lis, bins=bins, color='#348ABD', edgecolor='black')
plt.xlabel('LIS', fontname='Arial', fontsize=14)
plt.ylabel('Count', fontname='Arial', fontsize=14)
plt.title('LIS Distribution (Homodimer Matched)', fontname='Arial', fontsize=16)
plt.axvline(x=0.203, color='gray', linestyle='--', linewidth=1.5)
plt.grid(False)
plt.gca().set_facecolor('white')
plt.gcf().patch.set_facecolor('white')
plt.savefig(output_pdf, bbox_inches='tight', dpi=300, transparent=True)
plt.close()
print(f'已输出: {output_pdf}')

# 输出每个柱子的区间和数量到csv
bins_left = bins_edges[:-1]
bins_right = bins_edges[1:]
bins_df = pd.DataFrame({'bin_left': bins_left, 'bin_right': bins_right, 'count': counts.astype(int)})
bins_df.to_csv(output_bins_csv, index=False)
print(f'已输出: {output_bins_csv}')

# 统计LIS>=0.203的数量
count_ge_0203 = (lis >= 0.203).sum()
print(f'LIS >= 0.203 的数量: {count_ge_0203}') 