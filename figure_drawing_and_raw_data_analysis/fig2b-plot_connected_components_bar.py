import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# 读取数据
input_file = 'verifytable/heterodimer_batch123_matched_LIS_LIA_filtered.csv'
df = pd.read_csv(input_file)

# 统计每个菌的连通子图（元素>=3）数量
grouped = df.groupby('source')
result = {}
for source, group in grouped:
    G = nx.Graph()
    for n1, n2 in zip(group['SN1'], group['SN2']):
        if pd.notna(n1) and pd.notna(n2):
            G.add_edge(str(n1), str(n2))
    count = 0
    for comp in nx.connected_components(G):
        if len(comp) >= 3:
            count += 1
    result[source] = count

# 转为DataFrame并排序
res_df = pd.DataFrame(list(result.items()), columns=['菌名', '连通子图数(元素≥3)'])
res_df = res_df.sort_values('连通子图数(元素≥3)', ascending=False)

# 绘制柱状图
plt.figure(figsize=(12,6))
plt.bar(res_df['菌名'], res_df['连通子图数(元素≥3)'])
plt.xlabel('菌名')
plt.ylabel('连通子图数(元素≥3)')
plt.title('每个菌中元素≥3的连通子图个数统计')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('connected_components_bar.png', dpi=200)
plt.show()

print('统计和绘图已完成，图片已保存为 connected_components_bar.png') 