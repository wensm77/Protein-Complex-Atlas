import csv

def read_cluster_results(tsv_file):
    cluster_dict = {}
    
    with open(tsv_file, 'r') as file:
        for line in file:
            cluster, sequence = line.strip().split('\t')
            
          
            if cluster not in cluster_dict:
                cluster_dict[cluster] = []
            
          
            cluster_dict[cluster].append(sequence)
    
    return cluster_dict

def write_clusters_to_csv(cluster_dict, output_csv):
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['CHAINID', 'CLUSTER'])
        
        for cluster, sequences in sorted(cluster_dict.items()):
            for sequence in sequences:
                writer.writerow([sequence, cluster])

def main():
    tsv_file = '.../MPNN_data/foldseek/clu_cluster.tsv'  # 替换为你的TSV文件路径
    output_csv = '.../MPNN_data/clusters_with_representatives.csv'  # 输出CSV文件名

 
    cluster_dict = read_cluster_results(tsv_file)

    write_clusters_to_csv(cluster_dict, output_csv)
    print(f"Clusters with representatives written to {output_csv}")

if __name__ == "__main__":
    main()
