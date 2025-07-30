import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from matplotlib import font_manager

# 数据
species = ['Prokaryotic species', 'Human', 'Human-virus', 'Arabidopsis', 'Mouse']
counts = [626867, 17959, 4920, 5665, 13581]

# 按照数量从大到小排序
data_pairs = list(zip(species, counts))
data_pairs.sort(key=lambda x: x[1], reverse=True)
species_sorted, counts_sorted = zip(*data_pairs)

# 设置字体为 Arial
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial']
# 设置PDF字体嵌入为可编辑格式
plt.rcParams['pdf.fonttype'] = 42  # TrueType字体，可编辑
plt.rcParams['ps.fonttype'] = 42   # PostScript字体，可编辑

# 创建图形和轴，调整尺寸比例
fig, ax = plt.subplots(figsize=(12, 8))

# 创建对数柱状图，使用更细的柱子
bars = ax.bar(species_sorted, counts_sorted, color='#4A5568', width=0.4, alpha=0.9)

# 设置y轴为对数刻度，调整范围使柱子更短
ax.set_yscale('log')
ax.set_ylim(1000, 10000000)  # 扩大范围，使柱子相对更短

# 设置y轴刻度，使用标准数字表示法
ax.set_yticks([1000, 10000, 100000, 1000000, 10000000])
ax.set_yticklabels(['10³', '10⁴', '10⁵', '10⁶', '10⁷'])

# 设置标题和标签
ax.set_ylabel('Counts', fontsize=16, fontweight='bold', color='black')
ax.set_xlabel('')  # 不显示x轴标签

# 设置x轴和y轴标签，增大字体
plt.xticks(rotation=0, ha='center', fontsize=14, color='black')
plt.yticks(fontsize=14, color='black')

# 在每个柱子上显示数值，调整位置
for bar, count in zip(bars, counts_sorted):
    height = bar.get_height()
    # 将数字放在柱子上方
    y_pos = height * 1.15
    ax.text(bar.get_x() + bar.get_width()/2., y_pos,
            f'{count:,}',
            ha='center', va='bottom', fontsize=12, fontweight='bold', color='black')

# 移除上边框和右边框，保持左边框和下边框为黑色
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_color('black')
ax.spines['left'].set_linewidth(1.5)
ax.spines['bottom'].set_color('black')
ax.spines['bottom'].set_linewidth(1.5)

# 移除网格线
ax.grid(False)

# 设置坐标轴颜色和粗细
ax.tick_params(colors='black', width=1.5)

# 调整布局，增加边距
plt.tight_layout(pad=2.0)

# 保存为可编辑的PDF
plt.savefig('species_count_log_chart_editable.pdf', 
            format='pdf', 
            dpi=300, 
            bbox_inches='tight',
            transparent=False,
            facecolor='white')

# 显示图表
plt.show()

print("可编辑PDF图表已保存为 'species_count_log_chart_editable.pdf'") 