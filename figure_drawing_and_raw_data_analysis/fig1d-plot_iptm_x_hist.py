import pandas as pd
import matplotlib.pyplot as plt
import os

folder = os.path.dirname(os.path.abspath(__file__))
input_path = os.path.join(folder, 'heterodimer_batch123_matched_unique.csv')
output_pdf = os.path.join(folder, 'iptm_x_hist.pdf')
output_bins_csv = os.path.join(folder, 'iptm_x_hist_bins.csv')

# 读取数据
try:
    df = pd.read_csv(input_path, usecols=['iptm_x'])
except Exception:
    df = pd.read_csv(input_path)

# 去除缺失值
iptm_x = df['iptm_x'].dropna()

plt.figure(figsize=(8, 6))
bins = 15
counts, bins_edges, _ = plt.hist(iptm_x, bins=bins, color='c', edgecolor='black')
plt.xlabel('iptm_x', fontname='Arial', fontsize=14)
plt.ylabel('Count', fontname='Arial', fontsize=14)
plt.title('iptm_x Distribution', fontname='Arial', fontsize=16)
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