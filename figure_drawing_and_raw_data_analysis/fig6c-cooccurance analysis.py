#!/usr/bin/env python3
import os
import sys
import time

# === 配置区 ===
WORK_DIR      = "/Users/david/Desktop/m8"
UNIQUE_FILE   = os.path.join(WORK_DIR, "unique_first_column.txt")
OUTPUT_DIR    = "/Users/david/Desktop/m8/m8split_output"
OUTPUT_TSV    = os.path.join(WORK_DIR, "colocalization_scores3.tsv")
MIN_SET_SIZE  = 0    # 最小集合大小阈值
# =============

def load_genes(unique_file):
    with open(unique_file, 'r') as f:
        return [line.strip() for line in f if line.strip()]

def build_gene_data(genes, output_dir):
    gene_data = {}
    for gene in genes:
        files = [fn for fn in os.listdir(output_dir) if fn.startswith(f"{gene}_")]
        if not files:
            continue
        s = set()
        for fn in files:
            path = os.path.join(output_dir, fn)
            with open(path) as f:
                for line in f:
                    fields = line.strip().split()
                    if len(fields) < 2:
                        continue
                    genome = fields[1].split("_GCF")[0]
                    s.add(genome)
        gene_data[gene] = {
            "unique_set": s,
            "set_size":  len(s)
        }
    return gene_data

def filter_genes_by_size(genes, gene_data, min_size):
    kept = [g for g in genes if g in gene_data and gene_data[g]["set_size"] >= min_size]
    print(f"过滤后，满足集合大小 ≥ {min_size} 的基因有 {len(kept)} 个（共 {len(genes)} 个）")
    return kept

def compute_coloc(filtered_genes, gene_data, output_tsv):
    n = len(filtered_genes)
    total_pairs = n * (n - 1) // 2
    print(f"开始计算 {total_pairs} 对基因的共定位分数（交集/并集）。")

    pair_count = 0
    interval = max(total_pairs // 100, 1)
    start = time.time()

    with open(output_tsv, 'w') as out:
        out.write("gene1\tgene2\tscore\n")
        for i in range(n):
            g1 = filtered_genes[i]
            set1 = gene_data[g1]["unique_set"]
            for j in range(i+1, n):
                g2 = filtered_genes[j]
                set2 = gene_data[g2]["unique_set"]
                inter = len(set1 & set2)
                union = len(set1 | set2)  # 并集大小
                score = inter / union if union else 0.0
                out.write(f"{g1}\t{g2}\t{score:.4f}\n")
                pair_count += 1
                if pair_count % interval == 0 or pair_count == total_pairs:
                    pct = pair_count / total_pairs * 100
                    elapsed = time.time() - start
                    print(f"  已处理 {pair_count}/{total_pairs} 对 ({pct:.2f}%)，用时 {elapsed:.1f}s")

    print("计算完成，结果保存在：", output_tsv)

if __name__ == "__main__":
    genes = load_genes(UNIQUE_FILE)
    gene_data = build_gene_data(genes, OUTPUT_DIR)
    filtered = filter_genes_by_size(genes, gene_data, MIN_SET_SIZE)
    compute_coloc(filtered, gene_data, OUTPUT_TSV)