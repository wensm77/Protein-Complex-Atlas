import os
import re
import numpy as np
import pandas as pd
import csv
import pickle as pk
import scipy.sparse as sps
import argparse
from Bio import SeqIO


def parse_args1():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="Extract co-localized protein pairs.\nPlase cite: ", 
                        epilog="Co_localization Version: 1.0\n@ Cedric Liang (cedricliang@163.com)")
    parser.add_argument("-f","--fasta_path", type=str, required=True, help="The path of query fasta file. Type = string")
    parser.add_argument("-w","--work_path", type=str, required=True, help="The path of the co-localization output. Type = string")
    parser.add_argument("-o","--output_suffix", type=str, default="", help="The suffix of output file. Type = string [default = '']")
    parser.add_argument("-c","--cluster_size", type=int, default=2000, help="The minimum size of searched results. Type = integer [default = 2000]")
    parser.add_argument("-m","--min_points", type=float, default=0.6, help="The minimum score of searched results. Type = float [default = 0.6]")
    parser.add_argument("-s","--start_no", type=int, default=1, help="The start index of valid pairs. Type = integer [default = 1].")
    return parser.parse_args()


def min_distance_between_ranges(r1_start:int, r1_end:int, r2_start:int, r2_end:int)->int:
    # calculate the minimum distance between two ranges
    r1_start, r1_end = sorted([r1_start, r1_end])
    r2_start, r2_end = sorted([r2_start, r2_end])

    if r1_end < r2_start or r2_end < r1_start:
        dist1 = abs(r2_start - r1_end)
        dist2 = abs(r1_start - r2_end)
        return min(dist1, dist2)
    else:
        return 0

def min_distance_between_ranges1(r1_start:int, r1_end:int, r2_start:int, r2_end:int)->int:
    # calculate the minimum distance between two ranges considering the circularity. e.g. [45678, 123]
    max_num = 99999999
    if r1_start > r1_end and r2_start > r2_end:
        return 0
    if r1_start > r1_end:
        dist1 = min_distance_between_ranges(r1_start, max_num, r2_start, r2_end)
        dist2 = min_distance_between_ranges(1, r1_end, r2_start, r2_end)
        return min(dist1, dist2)
    if r2_start > r2_end:
        dist1 = min_distance_between_ranges(r1_start, r1_end, r2_start, max_num)
        dist2 = min_distance_between_ranges(r1_start, r1_end, 1, r2_end)
        return min(dist1, dist2)
    return min_distance_between_ranges(r1_start, r1_end, r2_start, r2_end)


def get_seq_pairs_from_query(argus):
    
    fasta_path = argus.fasta_path
    work_path = argus.work_path
    output_suffix = argus.output_suffix
    cluster_size = argus.cluster_size
    min_points = argus.min_points
    pairs_no = argus.start_no

    # Read the fasta file and create a dictionary to save it.
    sequence_dict = {}  # key = seq_id, value = sequence_string
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequence_dict[record.id] = str(record.seq)
    
    dict_path = os.path.join(work_path, "dict_folder")
    path4 = os.path.join(dict_path, 'seq_id_dict.pickle')
    path5 = os.path.join(dict_path, 'query_dict.pickle')
    with open(path4,'rb') as file:
        seq_id_dict = pk.load(file)
    with open(path5,'rb') as file:
        query_dict = pk.load(file)
       
    # Read the sparse matrix
    current_path = os.path.join(work_path, "output_folder_3")
    filename_mul = "ad_points_mul.npz"
    full_path_mul = os.path.join(current_path, filename_mul)

    ad_mat_mul_csr = sps.load_npz(full_path_mul)
    num_cols = ad_mat_mul_csr.shape[1]
    
    # Transform the CSR matrix to LIL matrix, in order to accelerate indexing
    ad_mat_mul = ad_mat_mul_csr.tolil()
    
    # Filter pairs and write to file
    pairs_name_to_be_written = []
    pairs_cnt = 0
    
    distance_thr = 3000
    row_indices, col_indices = ad_mat_mul.nonzero()
    for row, col in zip(row_indices.tolist(), col_indices.tolist()):
        if row < col:
            seqs_list1 = query_dict[row].split("_")
            seqs_list2 = query_dict[col].split("_")
            contig_1 = seqs_list1[0] if len(seqs_list1[0]) > 3 else '_'.join(seqs_list1[0:2])  # some contig name is like "NX_xxxxx"
            contig_2 = seqs_list2[0] if len(seqs_list2[0]) > 3 else '_'.join(seqs_list2[0:2])
            # Conditions：
            # 1. the mmseqs_search filtered result sequence numbers of two queries can not lower than cluster_size;
            # 2. the score matrix[i,j] >= min_points or the score matrix[j,i] >= min_points;
            # 3. two queries must share the same contig;
            # 4. the minimum distance of two queries <= distance_thr
            if len(seq_id_dict[row]) >= cluster_size and len(seq_id_dict[col]) >= cluster_size and \
                (ad_mat_mul[row, col] >= min_points or ad_mat_mul[col, row] >= min_points) and \
                contig_1 == contig_2 and min_distance_between_ranges1(int(seqs_list1[-2]), int(seqs_list1[-1]), int(seqs_list2[-2]), int(seqs_list2[-1])) <= distance_thr:
                list_tmp1 = [pairs_no, sequence_dict[query_dict[row]], sequence_dict[query_dict[col]], query_dict[row], query_dict[col], min_distance_between_ranges1(int(seqs_list1[-2]), int(seqs_list1[-1]), int(seqs_list2[-2]), int(seqs_list2[-1])), ad_mat_mul[row, col], ad_mat_mul[col, row], len(seq_id_dict[row]), len(seq_id_dict[col])]
                pairs_name_to_be_written.append(list_tmp1)
                pairs_cnt += 1
                pairs_no += 1 
    
    # write pairs into files
    out_suffix = f"output_seq_pairs_{output_suffix}.csv"
    save_file_path = os.path.join(work_path, out_suffix)
    header = ["Pair Index", "Sequence 1", "Sequence 2", "SequenceName 1", "SequenceName 2", "Min Distance", "Ad Points(row, col)", "Ad Points(col, row)", "ClusterSize1", "ClusterSize2"]
    if len(pairs_name_to_be_written) > 0:
        with open(save_file_path, 'w', newline='') as wcsvfile:
            writer = csv.writer(wcsvfile, delimiter='\t')
            writer.writerow(header)
            writer.writerows(pairs_name_to_be_written)                           
        pairs_name_to_be_written = []
        print(f"{pairs_cnt} type-pairs are filtered!")               
    else:
        print("No pairs are found!!!")

    return pairs_no


if __name__ == '__main__':

    args1 = parse_args1()    
    result = get_seq_pairs_from_query(args1)

    