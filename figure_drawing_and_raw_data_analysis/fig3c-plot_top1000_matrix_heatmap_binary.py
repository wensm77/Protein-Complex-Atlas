import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from matplotlib.colors import ListedColormap

def create_binary_matrix_heatmap():
    """创建二元颜色的聚类矩阵热图（蓝色-白色）"""
    
    # 文件路径
    csv_file = '115w_scores_plus_LIS_0409_human_fasta6_both_matched_LIS_LIA_filtered.csv'
    
    print("正在读取数据...")
    df = pd.read_csv(csv_file)
    print(f"总记录数: {len(df):,}")
    
    # 过滤掉空值
    df_clean = df.dropna(subset=['seq1_name', 'seq2_name'])
    print(f"有效连接记录数: {len(df_clean):,}")
    
    # 统计连接度数
    all_proteins = list(df_clean['seq1_name']) + list(df_clean['seq2_name'])
    protein_counts = Counter(all_proteins)
    
    # 选择Top 1000蛋白质
    top_proteins = [protein for protein, count in protein_counts.most_common(1000)]
    print(f"选择Top 1000蛋白质进行矩阵可视化")
    
    # 创建邻接矩阵
    protein_to_idx = {protein: i for i, protein in enumerate(top_proteins)}
    matrix = np.zeros((1000, 1000))
    
    print("正在构建邻接矩阵...")
    # 填充矩阵
    for idx, row in df_clean.iterrows():
        if idx % 5000 == 0:
            print(f"  处理进度: {idx:,}/{len(df_clean):,} ({idx/len(df_clean)*100:.1f}%)")
        
        seq1 = row['seq1_name']
        seq2 = row['seq2_name']
        
        if seq1 in protein_to_idx and seq2 in protein_to_idx:
            i = protein_to_idx[seq1]
            j = protein_to_idx[seq2]
            matrix[i, j] = 1
            matrix[j, i] = 1  # 对称矩阵
    
    # 计算连接统计
    total_connections = np.sum(matrix) / 2  # 除以2因为是对称矩阵
    max_possible = 1000 * 999 / 2
    density = total_connections / max_possible * 100
    
    print(f"\n矩阵统计:")
    print(f"  矩阵大小: 1000 × 1000")
    print(f"  实际连接数: {int(total_connections):,}")
    print(f"  连接密度: {density:.4f}%")
    
    # 创建二元颜色映射：白色(0) -> 蓝色(1)
    colors = ['white', '#1f77b4']  # 白色和蓝色
    cmap = ListedColormap(colors)
    
    # 创建聚类热图，显示聚类树，聚类树高度较矮
    g = sns.clustermap(matrix,
                       method='ward',
                       metric='euclidean',
                       cmap=cmap,
                       linewidths=0,  # 不显示网格线
                       figsize=(20, 20),
                       xticklabels=False,  # 不显示x轴标签
                       yticklabels=False,  # 不显示y轴标签
                       cbar=False,  # 不显示颜色条
                       dendrogram_ratio=(0.05, 0.05))  # 聚类树高度较矮
    # 强制去除colorbar和所有数字
    if hasattr(g, 'cax') and g.cax is not None:
        g.cax.set_visible(False)
    for ax in [g.ax_heatmap, g.ax_row_dendrogram, g.ax_col_dendrogram]:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_title("")
    
    # 保存高分辨率图片
    output_file = 'top2000_binary_clustered_matrix.pdf'
    plt.savefig(output_file, dpi=300, bbox_inches='tight', pad_inches=0)
    plt.close()
    
    print(f"二元聚类热图已保存: {output_file}")
    
    # 创建简化版本（降采样）
    create_binary_downsampled_heatmap(matrix, cmap)
    
    # 创建无聚类版本（按原始顺序）
    create_binary_original_heatmap(matrix, cmap)

def create_binary_downsampled_heatmap(matrix, cmap):
    """创建降采样版本的二元热图"""
    print("\n正在生成降采样版本热图...")
    
    # 每10个蛋白质取1个，得到200×200矩阵
    step = 10
    downsampled_matrix = matrix[::step, ::step]
    
    plt.figure(figsize=(12, 12))
    
    # 创建热图，完全无标签
    ax = sns.heatmap(downsampled_matrix,
                     cmap=cmap,
                     square=True,
                     linewidths=0,
                     xticklabels=False,
                     yticklabels=False,
                     cbar=False)
    
    # 去掉所有标签和标题
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title('')
    ax.set_xticks([])
    ax.set_yticks([])
    
    plt.tight_layout()
    plt.savefig('top2000_binary_downsampled_matrix.pdf', dpi=300, bbox_inches='tight', pad_inches=0)
    plt.close()
    
    print(f"降采样二元热图已保存: top2000_binary_downsampled_matrix.pdf")

def create_binary_original_heatmap(matrix, cmap):
    """创建无聚类的原始顺序二元热图"""
    print("\n正在生成原始顺序热图...")
    
    plt.figure(figsize=(20, 20))
    
    # 创建热图，完全无标签，不进行聚类
    ax = sns.heatmap(matrix,
                     cmap=cmap,
                     square=True,
                     linewidths=0,
                     xticklabels=False,
                     yticklabels=False,
                     cbar=False)
    
    # 去掉所有标签和标题
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title('')
    ax.set_xticks([])
    ax.set_yticks([])
    
    plt.tight_layout()
    plt.savefig('top2000_binary_original_matrix.pdf', dpi=300, bbox_inches='tight', pad_inches=0)
    plt.close()
    
    print(f"原始顺序二元热图已保存: top2000_binary_original_matrix.pdf")

if __name__ == "__main__":
    create_binary_matrix_heatmap() 