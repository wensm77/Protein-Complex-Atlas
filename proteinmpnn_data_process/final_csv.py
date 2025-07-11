import os
import pandas as pd

def get_chain_ids_from_directory(directory, suffixes):
    """从指定目录中获取符合后缀条件的文件ID名（包含后缀）"""
    chain_ids = []
    
    for filename in os.listdir(directory):
        for suffix in suffixes:
            if filename.endswith(suffix):
                chain_id = filename[:-len('.pt')]  # 去除 .pt 后缀
                chain_ids.append(chain_id)
                break
    
    return chain_ids

def main():
 
    directory = '.../MPNN_data/pt'
    updated_list_csv = '.../MPNN_data/List_new.csv'
    output_csv = '.../MPNN_data/list_new.csv'
    
 
    suffixes = ['_A.pt', '_B.pt']
    chain_ids = get_chain_ids_from_directory(directory, suffixes)
    
    # 读取Updated_List.csv文件
    updated_df = pd.read_csv(updated_list_csv)
    
    # 创建一个新的DataFrame
    new_df = pd.DataFrame(columns=['CHAINID', 'DEPOSITION', 'RESOLUTION', 'HASH', 'CLUSTER', 'SEQUENCE'])
    
 
    for chain_id in chain_ids:
        base_id = chain_id.split('_')[0]   
        matching_row = updated_df[updated_df['CHAINID'] == base_id]
        
        if not matching_row.empty:
            deposition = '2017/2/27'
            resolution = 2.0
            hash_value = matching_row['HASH'].values[0]
            cluster = matching_row['CLUSTER'].values[0]
            sequence = matching_row['SEQUENCE'].values[0]
            
       
            new_df.loc[len(new_df)] = [chain_id, deposition, resolution, hash_value, cluster, sequence]
        else:
            print(f"Matching row not found for {base_id} in Updated_List.csv")
 
    new_df.to_csv(output_csv, index=False)
    print(f"New CSV file with required columns written to {output_csv}")

if __name__ == "__main__":
    main()
