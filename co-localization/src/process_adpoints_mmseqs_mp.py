import re
import os
import csv
import math
import numpy as np
import argparse
import pickle as pk
import scipy.sparse as sps
from multiprocessing import Process, Manager


def make_contig_dict_simple(seq_dict):
    # pre-process data: create a dictionary to save it
    clust_type_dict = {}  # key = contig, value = list of seq_id
    for key in seq_dict.keys():
        split_list = key.split('_')
        contig = split_list[0]
        if contig:
            if contig not in clust_type_dict:
                clust_type_dict[contig] = []
            clust_type_dict[contig].append(key)
        else:
            print(f"Can not take Acc Number from '{key}'!!!")
    return clust_type_dict


# 2.
def generate_filterred_clust_rdd(res_file_path, length_thr):
    # res_file_path: m8_file
    # length_thr: the min threshold of the length of the cluster(0 currently, we don't use it)
    result_dict = {} 
    clust_type_dict = {} 
    query_id_dict = {} 
    seq_identity_dict = {} 
    start_index = -1 

    with open(res_file_path, "r") as f_res:
        for line in f_res:
            columns = line.strip().split('\t')
            query_seq = columns[0]
            if query_seq not in clust_type_dict:
                start_index += 1
                clust_type_dict[query_seq] = start_index
                result_dict[start_index] = []
                query_id_dict[start_index] = query_seq
                seq_identity_dict[start_index] = {}

            result_dict[start_index].append(columns[1])
            seq_identity_dict[start_index][columns[1]] = float(columns[2])           

    # filter the cluster whose length < length_thr
    flt_dict = {}
    seq_type_dict = {}
    rm_rdd_dict = {}
    seq_id_dict = {}
    query_dict = {}
    tmp_num = 0
    for key, value in result_dict.items():
        if len(value) >= length_thr:
            flt_dict[tmp_num] = value
            seq_id_dict[tmp_num] = seq_identity_dict[key]
            query_dict[tmp_num] = query_id_dict[key]
            for val in value:
                if val not in seq_type_dict:
                    seq_type_dict[val] = []
                seq_type_dict[val].append(tmp_num)
            tmp_num += 1
    return flt_dict, seq_type_dict, seq_id_dict, query_dict


def main_func(argus):
    csv_path = argus.clust_path
    length_thr = argus.len_min
    process_num = argus.process_num
    if (process_num < 1):
        process_num = 1
    rrs = argus.remove_redundancy_th
    save_path = argus.output_folder

    flt_clust_dict, seq_type_dict, seq_id_dict, query_dict = generate_filterred_clust_rdd(csv_path, length_thr)
    contig_dict = make_contig_dict_simple(seq_type_dict)

    clust_type = len(flt_clust_dict)  # the number of searched sequences after filtering
    interval = clust_type / process_num
    quotient = math.ceil(interval / 5)
    interval = quotient * 5

    def process_func(pos_s, pos_e, p_idx):
        # single process function
        # input: pos_s: starting index, pos_e: ending index(not included), p_idx: process index
        pos_e = min(pos_e, clust_type)  # protect index out of range
        if (pos_s >= pos_e):
            print(f"!!!{p_idx}!!!: pos_s is larger than pos_e!!! No need to record or process!")
            return
        # part of score matrix initialization with dimention pos_e-pos_s × N
        co_factor_single = sps.lil_matrix((pos_e - pos_s, clust_type))
        co_factor_multiple = sps.lil_matrix((pos_e - pos_s, clust_type))
        print(f"==={p_idx}:=== matrices have been initialized!!!It totally contains {clust_type} types!")

        seq_id_rdd_list = []
        for i in range(clust_type):
            tmp = sum(1 for val in seq_id_dict[i].values() if val <= rrs)
            seq_id_rdd_list.append(tmp)
        assert len(seq_id_rdd_list) == clust_type  # make sure the length of the list matches the number of columns

        # only the first process need to save the dict
        if p_idx == 0:
            dict_path = os.path.join(save_path, "dict_folder")
            if not os.path.exists(dict_path):
                os.makedirs(dict_path)
            path1 = os.path.join(dict_path, 'flt_clust_dict.pickle')
            path2 = os.path.join(dict_path, 'seq_type_dict.pickle')
            path3 = os.path.join(dict_path, 'contig_dict.pickle')
            path4 = os.path.join(dict_path, 'seq_id_dict.pickle')
            path5 = os.path.join(dict_path, 'query_dict.pickle')        
            with open(path1,'wb') as file:
                pk.dump(flt_clust_dict, file)
            with open(path2,'wb') as file:
                pk.dump(seq_type_dict, file)
            with open(path3,'wb') as file:
                pk.dump(contig_dict, file)
            with open(path4,'wb') as file:
                pk.dump(seq_id_dict, file)
            with open(path5,'wb') as file:
                pk.dump(query_dict, file)      
        
        for i in range(pos_s, pos_e):
            for seq1 in flt_clust_dict[i]:
                contig1 = seq1.split('_')[0]
                seq1_flag_list = [False] * clust_type
                for raw_seq in contig_dict[contig1]:
                    tmp_list = seq_type_dict[raw_seq]
                    for type_id in tmp_list:
                        if type_id > i and seq_id_dict[i][seq1] <= rrs and seq_id_dict[type_id][raw_seq] <= rrs and is_protein_colocalization_simple(seq1, raw_seq):
                            co_factor_multiple[i-pos_s, type_id] += 1
                            if not seq1_flag_list[type_id]:
                                co_factor_single[i-pos_s, type_id] += 1
                                seq1_flag_list[type_id] = True                           
       
        output_file1 = f"ad_factor_single{p_idx}.npz"
        output_file2 = f"ad_factor_multiple{p_idx}.npz"
        output_path2 = os.path.join(save_path, "output_folder_2")
        if not os.path.exists(output_path2):
            os.makedirs(output_path2)
        output_file1 = os.path.join(output_path2, output_file1)
        output_file2 = os.path.join(output_path2, output_file2)

        co_factor_single_csr = co_factor_single.tocsr()
        co_factor_multiple_csr = co_factor_multiple.tocsr()
        sps.save_npz(output_file1, co_factor_single_csr)
        sps.save_npz(output_file2, co_factor_multiple_csr)
        print(f"==={p_idx}:=== Well Done! Processing Finished!")
    

    # start multiple processes
    with Manager() as manager:
        processes = []
        for i in range(process_num):
            p_start = i * interval
            p_end = (i + 1) * interval
            p = Process(target=process_func, args=(p_start, p_end, i,))
            p.start()
            processes.append(p)
    for p in processes:
        p.join()

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


def is_protein_colocalization_simple(seq1:str, seq2:str)->bool:
    # Decide whether two protein sequences are colocalized
    # input: sequence1 name, sequence2 name
    # output: true/false
    if len(seq1) * len(seq2) == 0 or seq1 == seq2:
        return False 
    seq1_seg = seq1.split('_')
    seq2_seg = seq2.split('_')
    start_pos = 1 if len(seq1_seg[0]) > 3 else 2
    dist = min_distance_between_ranges(int(seq1_seg[start_pos]), int(seq1_seg[start_pos+1]), int(seq2_seg[start_pos]), int(seq2_seg[start_pos+1]))
    dist_th = 3000
    return dist <= dist_th


if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, 
                        description="Co-localization pair count matrix calculation.\nPlase cite: ", 
                        epilog="Co_localization Version: 1.0\n@ Cedric Liang (cedricliang@163.com)")
    parser.add_argument("-c", "--clust_path", type=str, required=True,
                        help="mmseqs_searched .m8 file path. Type = string")
    parser.add_argument("-l", "--len_min", type=int, default=10,
                        help="Clust type minimum length. The type contains more than this num will be processed. Type = integer [default: 10]")
    parser.add_argument("-p", "--process_num", type=int, default=15, help="Process number used for parallel computing. Type = integer [default: 15]")
    parser.add_argument("-r", "--remove_redundancy_th", help="Remove redundancy maximum threshold. If the similarity is larger than this threshold, the redundancy will be removed. Type = float [default: 0.96]", 
                        type = float, 
                        default=0.96)
    parser.add_argument("-o", "--output_folder", help="Output directory. Type = string [default = current work directory]", type=str, default="./")
    args = parser.parse_args()

    main_func(args)

