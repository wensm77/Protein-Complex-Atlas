import pandas as pd
import matplotlib.pyplot as plt
import os

folder = os.path.dirname(os.path.abspath(__file__))
input_path = os.path.join(folder, 'heterodimer_batch123_matched_unique.csv')
output_pdf = os.path.join(folder, 'LIS_hist_broken_y.pdf')
output_pdf_normal = os.path.join(folder, 'LIS_hist_normal.pdf')
output_bins_csv = os.path.join(folder, 'LIS_hist_bins.csv')

# 读取数据
try:
    df = pd.read_csv(input_path, usecols=['LIS'])
except Exception:
    df = pd.read_csv(input_path)

# 去除缺失值
lis = df['LIS'].dropna()

# 不分段直方图
plt.figure(figsize=(8, 6))
bins = 15
counts, bins_edges, _ = plt.hist(lis, bins=bins, color='#348ABD', edgecolor='black')
plt.xlabel('LIS', fontname='Arial', fontsize=14)
plt.ylabel('Count', fontname='Arial', fontsize=14)
plt.title('LIS Distribution', fontname='Arial', fontsize=16)
plt.axvline(x=0.203, color='red', linestyle='--', linewidth=1.5)
plt.grid(False)
plt.gca().set_facecolor('white')
plt.gcf().patch.set_facecolor('white')
plt.savefig(output_pdf_normal, bbox_inches='tight', dpi=300, transparent=True)
plt.close()
print(f'已输出: {output_pdf_normal}')

# 输出每个柱子的区间和数量到csv
bins_left = bins_edges[:-1]
bins_right = bins_edges[1:]
bins_df = pd.DataFrame({'bin_left': bins_left, 'bin_right': bins_right, 'count': counts.astype(int)})
bins_df.to_csv(output_bins_csv, index=False)
print(f'已输出: {output_bins_csv}')

# 分段断轴直方图
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 6), gridspec_kw={'height_ratios': [1, 2]})
ax2.hist(lis, bins=bins, color='#348ABD', edgecolor='black')
ax1.hist(lis, bins=bins, color='#348ABD', edgecolor='black')

# 设置y轴范围
ax2.set_ylim(0, 12000)
ax1.set_ylim(90000, 120000)

# 辅助线
ax1.axvline(x=0.203, color='red', linestyle='--', linewidth=1.5)
ax2.axvline(x=0.203, color='red', linestyle='--', linewidth=1.5)

# 隐藏上子图的x轴
ax1.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax1.tick_params(labeltop=False, top=False)
ax2.xaxis.tick_bottom()

# 断轴锯齿
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot([-0.01, +0.01], [-0.02, +0.02], **kwargs)
ax1.plot([0.99, 1.01], [-0.02, +0.02], **kwargs)
kwargs.update(transform=ax2.transAxes)
ax2.plot([-0.01, +0.01], [0.98, 1.02], **kwargs)
ax2.plot([0.99, 1.01], [0.98, 1.02], **kwargs)

# 字体和标签
ax2.set_xlabel('LIS', fontname='Arial', fontsize=14)
ax2.set_ylabel('Count', fontname='Arial', fontsize=14)
ax1.set_ylabel('Count', fontname='Arial', fontsize=14)
fig.suptitle('LIS Distribution (Broken Y)', fontname='Arial', fontsize=16)

for ax in [ax1, ax2]:
    ax.set_facecolor('white')
fig.patch.set_facecolor('white')

plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig(output_pdf, bbox_inches='tight', dpi=300, transparent=True)
plt.close()
print(f'已输出: {output_pdf}') 