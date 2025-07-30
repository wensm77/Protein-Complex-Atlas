import os
import sys
import json
import glob
import argparse
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from scipy.spatial.distance import pdist, squareform

###############################################################################
# Helper: Parse the PDB once and return the structure
###############################################################################
def get_structure(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("cf", pdb_file)
    return structure

###############################################################################
# 1) Parse PDB to identify chain lengths (using the structure)
###############################################################################
def parse_pdb_chains_from_structure(structure):
    """
    Given a parsed structure, returns:
      - chain_order: list of chain IDs (e.g. ['A','B','C',...])
      - chain_lengths: list of residue counts for each chain.
    (No filtering for water residues.)
    """
    chain_order = []
    chain_lengths = []
    for chain in structure.get_chains():
        chain_order.append(chain.id)
        # Count every residue (do not skip water)
        count = sum(1 for residue in chain)
        chain_lengths.append(count)
    return chain_order, chain_lengths

###############################################################################
# 2) Read ColabFold JSON
###############################################################################
def read_colabfold_json(json_file):
    """
    Reads a typical ColabFold JSON that has keys:
      - 'pae'   => NxN array
      - 'ptm'   => float
      - 'iptm'  => float
      - 'plddt' => array of length N
    """
    with open(json_file, 'r') as f:
        data = json.load(f)

    pae = np.array(data.get('pae', []), dtype=float)
    pae = np.nan_to_num(pae)
    ptm = data.get('ptm', 0.0)
    iptm = data.get('iptm', 0.0)
    plddt = np.array(data.get('plddt', []), dtype=float)

    return {
        "pae": pae,
        "ptm": ptm,
        "iptm": iptm,
        "plddt": plddt,
    }

###############################################################################
# 3) Transform PAE => LIS
###############################################################################
def transform_pae_matrix(pae_matrix, pae_cutoff):
    """
    For each PAE value < cutoff, compute LIS = 1 - (PAE/cutoff); else 0.
    """
    transformed = np.zeros_like(pae_matrix)
    mask = (pae_matrix < pae_cutoff)
    transformed[mask] = 1.0 - (pae_matrix[mask] / pae_cutoff)
    return transformed

###############################################################################
# 4) Calculate contact map from structure
###############################################################################
def calculate_contact_map_from_structure(structure, distance_threshold=8.0):
    coords = []
    for chain in structure.get_chains():
        for residue in chain:
            if residue.has_id('CB'):
                coords.append(residue['CB'].coord)
            elif residue.has_id('CA'):
                coords.append(residue['CA'].coord)
    
    if len(coords) == 0:
        return np.zeros((0, 0), dtype=int)
    
    coords = np.array(coords)
    distmat = squareform(pdist(coords))
    contact_map = (distmat < distance_threshold).astype(int)
    return contact_map

###############################################################################
# 5) Calculate mean LIS (or cLIS) for each chain pair sub-block
###############################################################################
def calculate_mean_lis(transformed_map, subunit_number):
    """
    For each submatrix (chain_i x chain_j) in transformed_map,
    average the values > 0.
    """
    cum_lengths = np.cumsum(subunit_number)
    n_chains = len(subunit_number)
    mean_matrix = np.zeros((n_chains, n_chains))
    starts = np.concatenate(([0], cum_lengths[:-1]))
    for i in range(n_chains):
        for j in range(n_chains):
            si, ei = starts[i], cum_lengths[i]
            sj, ej = starts[j], cum_lengths[j]
            block = transformed_map[si:ei, sj:ej]
            positive = block[block > 0]
            mean_matrix[i, j] = positive.mean() if len(positive) > 0 else 0.0
    return mean_matrix

###############################################################################
# 6) Main Analysis Function: Build and return DataFrame
###############################################################################
def colabfold_lis_analysis_to_df(pdb_files, json_files, pae_cutoff=12.0, distance_cutoff=8.0):
    """
    For each (pdb, json) pair, compute metrics for every chain pair:
      - LIA: count of positions with transformed PAE > 0.
      - LIR: count of unique residue indices (rows + columns) in the submatrix.
      - cLIA and cLIR: same as above but for contact-based (cLIA) map.
      - LIS: mean of non-zero values in submatrix of transformed PAE.
      - cLIS: mean of submatrix from contact-filtered map.
      - iLIS = sqrt(LIS * cLIS)
      - iptm, ptm, and average plddt from JSON.
      - Also extracts 'rank' and 'model' from the pdb file name.
      - File name is extracted as the part before "_relaxed_".
      - For each chain pair, returns LIR_indices_A, LIR_indices_B, cLIR_indices_A, and cLIR_indices_B.
    """
    all_rows = []
    for pdb_file, json_file in zip(pdb_files, json_files):
        folder_name = os.path.basename(os.path.dirname(pdb_file))
        base_name = os.path.basename(pdb_file)

        # Extract file_name as the part before "_relaxed_"
        file_name = base_name.split("_relaxed_")[0] if "_relaxed_" in base_name else base_name.split(".")[0]

        # Extract rank
        rank = None
        if "_rank_" in base_name:
            try:
                rank = int(base_name.split("_rank_")[1].split("_")[0])
            except ValueError:
                rank = None

        # Extract model number
        model_val = None
        if "model_" in base_name:
            try:
                model_val = base_name.split("model_")[1].split("_")[0]
            except IndexError:
                model_val = None
        
        print(f"处理文件: {folder_name}")
        
        # Parse the PDB structure once
        try:
            structure = get_structure(pdb_file)
            chain_order, chain_lengths = parse_pdb_chains_from_structure(structure)
            total_res = sum(chain_lengths)
        except Exception as e:
            print(f"[错误] 无法解析PDB文件 {pdb_file}: {e}")
            continue
        
        # Read JSON for PAE, ptm, iptm, and plddt
        try:
            data = read_colabfold_json(json_file)
            pae_matrix = data["pae"]
            ptm = data["ptm"]
            iptm = data["iptm"]
            plddt = data["plddt"]
            confidence = 0.8 * iptm + 0.2 * ptm
            plddt_avg = float(plddt.mean()) if plddt.size > 0 else 0.0
        except Exception as e:
            print(f"[错误] 无法读取JSON文件 {json_file}: {e}")
            continue

        if pae_matrix.shape[0] != total_res:
            print(f"[警告] PAE矩阵大小 {pae_matrix.shape} != 总残基数 {total_res} 对于文件 {pdb_file}")
        
        # Transform PAE to get the LIS map (LIA map)
        transformed_pae = transform_pae_matrix(pae_matrix, pae_cutoff)
        lia_map = (transformed_pae > 0).astype(int)
        
        # Calculate contact map (from the structure) and then cLIA map
        contact_map = calculate_contact_map_from_structure(structure, distance_cutoff)
        clia_map = ((transformed_pae > 0) & (contact_map == 1)).astype(int)
        
        # Compute mean values for each chain pair submatrix
        mean_lis = calculate_mean_lis(transformed_pae, chain_lengths)
        mean_clis = calculate_mean_lis(np.where(clia_map > 0, transformed_pae, 0), chain_lengths)
        iLIS_matrix = np.sqrt(mean_lis * mean_clis)
        
        n_sub = len(chain_lengths)
        cum_lengths = np.cumsum(chain_lengths)
        starts = np.concatenate(([0], cum_lengths[:-1]))
        
        for i in range(n_sub):
            for j in range(n_sub):
                si, ei = starts[i], cum_lengths[i]
                sj, ej = starts[j], cum_lengths[j]
                sub_lia = lia_map[si:ei, sj:ej]
                sub_clia = clia_map[si:ei, sj:ej]
                
                # Count-based metrics
                LIA_val = np.count_nonzero(sub_lia)
                LIR_val = len(np.unique(np.where(sub_lia > 0)[0])) + len(np.unique(np.where(sub_lia > 0)[1]))
                cLIA_val = np.count_nonzero(sub_clia)
                cLIR_val = len(np.unique(np.where(sub_clia > 0)[0])) + len(np.unique(np.where(sub_clia > 0)[1]))
                
                # Local residue indices (1-based)
                LIR_indices_A = np.unique(np.where(sub_lia > 0)[0].flatten() + 1).tolist()
                LIR_indices_B = np.unique(np.where(sub_lia > 0)[1].flatten() + 1).tolist()
                cLIR_indices_A = np.unique(np.where(sub_clia > 0)[0].flatten() + 1).tolist()
                cLIR_indices_B = np.unique(np.where(sub_clia > 0)[1].flatten() + 1).tolist()

                LIS_val = mean_lis[i, j]
                cLIS_val = mean_clis[i, j]
                iLIS_val = np.sqrt(LIS_val * cLIS_val)
                
                row_dict = {
                    'folder_name': folder_name,
                    'file_name': file_name,
                    'chain_1': i + 1,
                    'chain_2': j + 1,
                    'rank': rank,
                    'model': model_val,
                    'iLIS': iLIS_val,
                    'LIS': LIS_val,
                    'cLIS': cLIS_val,
                    'iptm': iptm,
                    'confidence': confidence,
                    'LIA': LIA_val,
                    'LIR': LIR_val,
                    'cLIA': cLIA_val,
                    'cLIR': cLIR_val,
                    'LIR_indices_A': LIR_indices_A,
                    'LIR_indices_B': LIR_indices_B,
                    'cLIR_indices_A': cLIR_indices_A,
                    'cLIR_indices_B': cLIR_indices_B,
                    'ptm': ptm,
                    'plddt': plddt_avg
                }
                all_rows.append(row_dict)
    
    if not all_rows:
        print("没有找到有效的文件对进行分析")
        return pd.DataFrame()
        
    df_merged = pd.DataFrame(all_rows)

    # Optional grouping to average symmetric pairs (e.g. (1,2) and (2,1))
    def union_of_lists(series):
        combined = set()
        for lst in series.dropna():
            combined.update(lst)
        return sorted(combined)
    
    # Add a helper column 'chain_pair' (sorted tuple)
    df_merged['chain_pair'] = df_merged.apply(lambda row: tuple(sorted((row['chain_1'], row['chain_2']))), axis=1)
    numeric_cols = ['iLIS', 'LIS', 'cLIS', 'LIA', 'cLIA', 'LIR', 'cLIR', 'iptm', 'confidence', 'ptm', 'plddt']
    list_cols = ['LIR_indices_A', 'LIR_indices_B', 'cLIR_indices_A', 'cLIR_indices_B']
    agg_dict = {col: 'mean' for col in numeric_cols}
    agg_dict.update({col: union_of_lists for col in list_cols})
    for c in ['file_name', 'rank', 'model', 'folder_name']:
        agg_dict[c] = 'first'
    group_cols = ['chain_pair', 'file_name', 'rank', 'model', 'folder_name']
    
    # Perform the grouping operation
    df_grouped = df_merged.groupby(group_cols, as_index=False).agg(agg_dict)

    # Round values as needed
    df_grouped['LIS'] = df_grouped['LIS'].round(3)
    df_grouped['cLIS'] = df_grouped['cLIS'].round(3)
    df_grouped['iLIS'] = df_grouped['iLIS'].round(3)
    df_grouped['plddt'] = df_grouped['plddt'].round(2)
    df_grouped['chain_1'] = df_grouped['chain_pair'].apply(lambda x: x[0])
    df_grouped['chain_2'] = df_grouped['chain_pair'].apply(lambda x: x[1])
    
    final_cols = [
        'file_name', 'chain_1', 'chain_2', 'rank', 'model',
        'iLIS', 'LIS', 'cLIS', 'LIA', 'cLIA', 'LIR', 'cLIR',
        'LIR_indices_A', 'LIR_indices_B', 'cLIR_indices_A', 'cLIR_indices_B',
        'iptm', 'confidence', 'ptm', 'plddt', 'folder_name'
    ]
    df_final = df_grouped[final_cols].sort_values(['folder_name', 'chain_1', 'chain_2'])
    
    return df_final

###############################################################################
# Find files function
###############################################################################
def find_colabfold_files(directory):
    """
    Find matching JSON and PDB files from ColabFold results in subdirectories.
    Updated to find rank_001 files regardless of model number.
    """
    file_pairs = []
    
    # 获取所有子文件夹
    subdirs = [d for d in os.listdir(directory) if os.path.isdir(os.path.join(directory, d)) and d.endswith('_batch_result')]
    
    print(f"找到 {len(subdirs)} 个子文件夹:")
    for subdir in subdirs:
        print(f"  - {subdir}")
    
    for subdir in subdirs:
        subdir_path = os.path.join(directory, subdir)
        print(f"\n处理子文件夹: {subdir}")
        
        json_files = {}
        pdb_files = {}
        
        # 在子文件夹中查找rank_001文件
        for file in os.listdir(subdir_path):
            if file.endswith('.json') and '_scores_rank_001_alphafold2_multimer_v3_' in file:
                # Extract base name (everything before _scores_rank_001)
                base_name = file.split('_scores_rank_001_alphafold2_multimer_v3_')[0]
                json_files[base_name] = file
                
            elif file.endswith('.pdb') and '_relaxed_rank_001_alphafold2_multimer_v3_' in file:
                # Extract base name (everything before _relaxed_rank_001)
                base_name = file.split('_relaxed_rank_001_alphafold2_multimer_v3_')[0]
                pdb_files[base_name] = file
        
        # Find matching pairs
        for base_name in json_files:
            if base_name in pdb_files:
                json_path = os.path.join(subdir_path, json_files[base_name])
                pdb_path = os.path.join(subdir_path, pdb_files[base_name])
                file_pairs.append((json_path, pdb_path))
                print(f"  ✓ 找到文件对:")
                print(f"    JSON: {json_files[base_name]}")
                print(f"    PDB: {pdb_files[base_name]}")
            else:
                print(f"  ⚠ 警告: 未找到匹配的PDB文件 {json_files[base_name]} 在 {subdir}")
        
        # Check for orphaned PDB files
        for base_name in pdb_files:
            if base_name not in json_files:
                print(f"  ⚠ 警告: 未找到匹配的JSON文件 {pdb_files[base_name]} 在 {subdir}")
    
    return file_pairs

def calculate_contact_map(pdb_file, distance_threshold=8.0):
    """
    直接从PDB文件计算接触图
    """
    from Bio.PDB import PDBParser
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    
    # 使用现有的函数计算接触图
    return calculate_contact_map_from_structure(structure, distance_threshold)

def extract_chain_info(pdb_file):
    """
    从PDB文件中提取链信息，返回每个链的起始和结束残基索引
    注意：按残基计数，不是原子计数
    """
    chain_info = {}
    current_chain = None
    current_residue = None
    chain_start = 0
    residue_count = 0
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                chain_id = line[21]
                residue_num = int(line[22:26].strip())
                
                # 如果是新链的开始
                if current_chain != chain_id:
                    # 如果不是第一个链，保存前一个链的信息
                    if current_chain is not None:
                        chain_info[current_chain] = (chain_start, residue_count - 1)
                    
                    # 开始新链
                    current_chain = chain_id
                    chain_start = residue_count
                    current_residue = residue_num
                    residue_count += 1
                    
                # 如果是同一链中的新残基
                elif residue_num != current_residue:
                    current_residue = residue_num
                    residue_count += 1
        
        # 保存最后一个链的信息
        if current_chain is not None:
            chain_info[current_chain] = (chain_start, residue_count - 1)
    
    return chain_info

def analyze_chain_pair(chain1, chain2, chain_info, transformed_pae, contact_map, distance_cutoff):
    """
    分析两个链之间的相互作用指标，严格按照原始代码逻辑
    """
    start1, end1 = chain_info[chain1]
    start2, end2 = chain_info[chain2] 
    
    # 提取链对之间的子矩阵
    pae_block = transformed_pae[start1:end1+1, start2:end2+1]
    contact_block = contact_map[start1:end1+1, start2:end2+1]
    
    # 创建LIA map (transformed_pae > 0的位置)
    lia_block = (pae_block > 0).astype(int)
    
    # 创建cLIA map (transformed_pae > 0 且 contact_map == 1的位置)
    clia_block = ((pae_block > 0) & (contact_block == 1)).astype(int)
    
    # LIS: 子矩阵中所有 > 0 的transformed_pae值的平均值
    positive_pae = pae_block[pae_block > 0]
    lis = positive_pae.mean() if len(positive_pae) > 0 else 0.0
    
    # cLIS: 接触区域中transformed_pae值的平均值
    contact_pae_values = pae_block[clia_block > 0]
    clis = contact_pae_values.mean() if len(contact_pae_values) > 0 else 0.0
    
    # iLIS: sqrt(LIS * cLIS)
    ilis = np.sqrt(lis * clis) if lis > 0 and clis > 0 else 0.0
    
    # LIA: transformed_pae > 0的位置数量
    lia = np.count_nonzero(lia_block)
    
    # cLIA: 接触区域中transformed_pae > 0的位置数量
    clia = np.count_nonzero(clia_block)
    
    # LIR: LIA区域中涉及的唯一残基数量 (行索引 + 列索引)
    lia_indices = np.where(lia_block > 0)
    lir = len(np.unique(lia_indices[0])) + len(np.unique(lia_indices[1])) if len(lia_indices[0]) > 0 else 0
    
    # cLIR: cLIA区域中涉及的唯一残基数量
    clia_indices = np.where(clia_block > 0)
    clir = len(np.unique(clia_indices[0])) + len(np.unique(clia_indices[1])) if len(clia_indices[0]) > 0 else 0
    
    # 获取残基索引 (1-based, 相对于链内)
    lia_residues_chain1 = np.unique(lia_indices[0] + 1).tolist() if len(lia_indices[0]) > 0 else []
    lia_residues_chain2 = np.unique(lia_indices[1] + 1).tolist() if len(lia_indices[0]) > 0 else []
    clia_residues_chain1 = np.unique(clia_indices[0] + 1).tolist() if len(clia_indices[0]) > 0 else []
    clia_residues_chain2 = np.unique(clia_indices[1] + 1).tolist() if len(clia_indices[0]) > 0 else []
    
    result = {
        'chain1': chain1,
        'chain2': chain2,
        'iLIS': ilis,
        'LIS': lis,
        'cLIS': clis,
        'LIA': lia,
        'cLIA': clia,
        'LIR': lir,
        'cLIR': clir,
        'total_contacts': np.sum(contact_block),
        'total_residues': (end1 - start1 + 1) * (end2 - start2 + 1),
        'chain1_length': end1 - start1 + 1,
        'chain2_length': end2 - start2 + 1,
        'LIR_indices_A': lia_residues_chain1,
        'LIR_indices_B': lia_residues_chain2,
        'cLIR_indices_A': clia_residues_chain1,
        'cLIR_indices_B': clia_residues_chain2
    }
    
    return result

def analyze_single_pair(json_file, pdb_file, distance_cutoff=8.0):
    """
    分析单个JSON/PDB文件对，严格按照原始代码逻辑
    """
    try:
        folder_name = os.path.basename(os.path.dirname(pdb_file))
        base_name = os.path.basename(pdb_file)

        # Extract file_name as the part before "_unrelaxed_" or "_relaxed_"
        if "_unrelaxed_" in base_name:
            file_name = base_name.split("_unrelaxed_")[0]
        elif "_relaxed_" in base_name:
            file_name = base_name.split("_relaxed_")[0]
        else:
            file_name = base_name

        # Extract rank
        rank = None
        if "_rank_" in base_name:
            try:
                rank = int(base_name.split("_rank_")[1].split("_")[0])
            except ValueError:
                rank = None

        # Extract model number
        model_val = None
        if "model_" in base_name:
            try:
                model_val = base_name.split("model_")[1].split("_")[0]
            except IndexError:
                model_val = None
        
        # 加载JSON数据
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        # 提取PAE数据 - ColabFold格式使用'pae'字段
        if 'pae' in data:
            pae_data = data['pae']
        elif 'predicted_aligned_error' in data:
            pae_data = data['predicted_aligned_error']
        else:
            raise ValueError("No PAE data found in JSON file")
            
        pae_array = np.array(pae_data)
        
        # 使用标准的LIS分析PAE转换公式：LIS = 1 - (PAE/12.0) for PAE < 12.0, else 0
        pae_cutoff = 12.0
        transformed_pae = np.zeros_like(pae_array)
        mask = (pae_array < pae_cutoff)
        transformed_pae[mask] = 1.0 - (pae_array[mask] / pae_cutoff)
        
        # 提取置信度分数
        iptm = data.get('iptm', 0)
        ptm = data.get('ptm', 0) 
        confidence = 0.8 * iptm + 0.2 * ptm
        
        # 计算平均plddt - ColabFold格式使用'plddt'字段
        if 'plddt' in data:
            plddt_array = np.array(data['plddt'])
            mean_plddt = np.mean(plddt_array)
        elif 'max_pae' in data and isinstance(data['max_pae'], list):
            plddt_array = np.array(data['max_pae'])
            mean_plddt = np.mean(plddt_array)
        else:
            mean_plddt = 0
        
        # 获取链信息
        chain_info = extract_chain_info(pdb_file)
        
        # 计算接触图
        contact_map = calculate_contact_map(pdb_file, distance_cutoff)
        
        results = []
        chains = list(chain_info.keys())
        
        # 分析所有链对
        for i, chain1 in enumerate(chains):
            for j, chain2 in enumerate(chains):
                if i < j:  # 避免重复分析
                    result = analyze_chain_pair(
                        chain1, chain2, 
                        chain_info, 
                        transformed_pae, 
                        contact_map, 
                        distance_cutoff
                    )
                    
                    # 添加文件信息和置信度分数，按照原始代码格式
                    result.update({
                        'folder_name': folder_name,
                        'file_name': file_name,
                        'chain_1': i + 1,  # 1-based chain numbering
                        'chain_2': j + 1,
                        'rank': rank,
                        'model': model_val,
                        'iptm': iptm,
                        'confidence': confidence,
                        'ptm': ptm,
                        'plddt': mean_plddt
                    })
                    
                    # 移除不需要的字段
                    result.pop('chain1', None)
                    result.pop('chain2', None)
                    
                    results.append(result)
        
        return results
        
    except Exception as e:
        print(f"错误处理文件 {json_file}: {str(e)}")
        return []

def save_results(results, output_file):
    """
    保存结果到CSV文件，按照原始代码格式
    """
    if not results:
        print("没有结果可保存")
        return
    
    # 转换为DataFrame
    df = pd.DataFrame(results)
    
    # Round values as needed (按照原始代码)
    if 'LIS' in df.columns:
        df['LIS'] = df['LIS'].round(3)
    if 'cLIS' in df.columns:
        df['cLIS'] = df['cLIS'].round(3)
    if 'iLIS' in df.columns:
        df['iLIS'] = df['iLIS'].round(3)
    if 'plddt' in df.columns:
        df['plddt'] = df['plddt'].round(2)
    
    # 按照原始代码的列顺序
    final_cols = [
        'file_name', 'chain_1', 'chain_2', 'rank', 'model',
        'iLIS', 'LIS', 'cLIS', 'LIA', 'cLIA', 'LIR', 'cLIR',
        'LIR_indices_A', 'LIR_indices_B', 'cLIR_indices_A', 'cLIR_indices_B',
        'iptm', 'confidence', 'ptm', 'plddt', 'folder_name'
    ]
    
    # 只保留存在的列
    available_cols = [col for col in final_cols if col in df.columns]
    df_final = df[available_cols].sort_values(['chain_1', 'chain_2'])
    
    # 保存到CSV
    df_final.to_csv(output_file, index=False)
    print(f"结果已保存到: {output_file}")
    print(f"包含 {len(df_final)} 行数据")

###############################################################################
# Main function
###############################################################################
def main():
    parser = argparse.ArgumentParser(description='批量分析ColabFold结果的LIS指标')
    parser.add_argument('--directory', '-d', default='.', 
                       help='要搜索的目录 (默认: 当前目录)')
    parser.add_argument('--output', '-o', default='colabfold_lis_analysis_summary.csv',
                       help='输出CSV文件名 (默认: colabfold_lis_analysis_summary.csv)')
    parser.add_argument('--cutoff', '-c', type=float, default=8.0,
                       help='接触距离阈值，单位Å (默认: 8.0)')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("ColabFold LIS 批量分析工具")
    print("=" * 60)
    print(f"搜索目录: {args.directory}")
    print(f"接触距离阈值: {args.cutoff} Å")
    print(f"输出文件: {args.output}")
    print("=" * 60)
    
    # 查找所有匹配的文件
    file_pairs = find_colabfold_files(args.directory)
    
    if not file_pairs:
        print("\n❌ 错误: 未找到任何匹配的文件对")
        print("请确保目录中包含以下格式的文件:")
        print("- *_scores_rank_001_alphafold2_multimer_v3_model_X_seed_000.json")
        print("- *_relaxed_rank_001_alphafold2_multimer_v3_model_X_seed_000.pdb")
        return
    
    print(f"\n✅ 找到 {len(file_pairs)} 个文件对")
    print("=" * 60)
    
    # 分析所有文件
    all_results = []
    successful_count = 0
    failed_count = 0
    
    for i, (json_file, pdb_file) in enumerate(file_pairs, 1):
        print(f"\n📁 处理文件 {i}/{len(file_pairs)}:")
        print(f"   JSON: {os.path.basename(json_file)}")
        print(f"   PDB:  {os.path.basename(pdb_file)}")
        
        try:
            results = analyze_single_pair(json_file, pdb_file, args.cutoff)
            if results:
                all_results.extend(results)
                successful_count += 1
                print(f"   ✅ 成功分析，得到 {len(results)} 个链对结果")
            else:
                failed_count += 1
                print(f"   ❌ 分析失败，未得到结果")
        except Exception as e:
            failed_count += 1
            print(f"   ❌ 处理文件时出错: {e}")
    
    # 保存结果
    print("\n" + "=" * 60)
    if all_results:
        save_results(all_results, args.output)
        print(f"✅ 批量分析完成!")
        print(f"📊 统计信息:")
        print(f"   - 成功处理: {successful_count} 个文件")
        print(f"   - 处理失败: {failed_count} 个文件")
        print(f"   - 总结果数: {len(all_results)} 行")
        print(f"   - 输出文件: {args.output}")
    else:
        print("❌ 没有成功分析任何文件")
    print("=" * 60)

if __name__ == "__main__":
    main() 