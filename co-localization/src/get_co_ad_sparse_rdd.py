import os
import pickle as pk
import numpy as np
import math
import scipy.sparse as sps
import argparse
from datetime import datetime 

def main(argu):
    curr_folder = argu.work_path
    # read dictionaries
    dict_path = os.path.join(curr_folder, "dict_folder")
    path1 = os.path.join(dict_path, 'flt_clust_dict.pickle')
    path2 = os.path.join(dict_path, 'seq_type_dict.pickle')
    path3 = os.path.join(dict_path, 'contig_dict.pickle')
    path4 = os.path.join(dict_path, 'seq_id_dict.pickle')
    with open(path1,'rb') as file:
        flt_clust_dict = pk.load(file)
    with open(path2,'rb') as file:
        seq_type_dict = pk.load(file)
    with open(path3,'rb') as file:
        contig_dict = pk.load(file)
    with open(path4,'rb') as file:
        seq_id_dict = pk.load(file)
    
    # read co-pair count matrix, which is sparse and upper triangular
    non_zero_elements = []
    matrix_path = os.path.join(curr_folder, "output_folder_2")
    index = 0  # co-pair count sparse matrix index initialized to 0
    filename = f"ad_factor_multiple{index}.npz"
    # load the first matrix and get its shape
    full_path = os.path.join(matrix_path, filename)
    first_matrix = sps.load_npz(full_path)
    num_cols = first_matrix.shape[1]
    # initialize two sparse matrices to store the merged results
    ad_mat_sgl = sps.lil_matrix((num_cols, num_cols))  
    ad_mat_mul = sps.lil_matrix((num_cols, num_cols))
    start_row = 0  # the starting row index of the current matrix

    # loop through the count matrices and union them into merged matrix above-mentioned
    while True:
        filename_sgl = f"ad_factor_single{index}.npz"
        filename_mul = f"ad_factor_multiple{index}.npz"
        full_path_sgl = os.path.join(matrix_path, filename_sgl)
        full_path_mul = os.path.join(matrix_path, filename_mul)
        if os.path.exists(full_path_sgl):           
            matrix1 = sps.load_npz(full_path_sgl)
            matrix2 = sps.load_npz(full_path_mul)
            end_row = start_row + matrix1.shape[0]
            ad_mat_sgl[start_row:end_row, :] = matrix1
            ad_mat_mul[start_row:end_row, :] = matrix2
            start_row = end_row
            index += 1
        else:
            break 
    
    seq_id_rdd_list = []
    for i in range(num_cols):
        tmp = sum(1 for val in seq_id_dict[i].values() if val <= 0.9)
        seq_id_rdd_list.append(tmp)
    assert len(seq_id_rdd_list) == num_cols  # make sure the length of the list matches the number of columns
        
    # calculate and save co-localization score matrix, and save single npz every interval rows
    current_path = os.path.join(curr_folder, "output_folder_3")
    if not os.path.exists(current_path):
        os.makedirs(current_path)
    matrix_sgl = sps.lil_matrix((num_cols, num_cols), dtype=float)
    matrix_mul = sps.lil_matrix((num_cols, num_cols), dtype=float)
    row_indices, col_indices = ad_mat_mul.nonzero()  # get nonzero indices
    for row, col in zip(row_indices.tolist(), col_indices.tolist()):
        tmp1 = seq_id_rdd_list[row]
        tmp2 = seq_id_rdd_list[col]
        matrix_sgl[row, col] = ad_mat_sgl[row, col] / tmp1 if tmp1 > 0 else 0.0
        matrix_mul[row, col] = ad_mat_mul[row, col] / tmp1 if tmp1 > 0 else 0.0
        matrix_sgl[col, row] = ad_mat_sgl[row, col] / tmp2 if tmp2 > 0 else 0.0
        matrix_mul[col, row] = ad_mat_mul[row, col] / tmp2 if tmp2 > 0 else 0.0
    sgl_path = os.path.join(current_path, "ad_points_sgl.npz") 
    mul_path = os.path.join(current_path, "ad_points_mul.npz")
    matrix_sgl_csr = matrix_sgl.tocsr()
    matrix_mul_csr = matrix_mul.tocsr()
    sps.save_npz(sgl_path, matrix_sgl_csr)
    sps.save_npz(mul_path, matrix_mul_csr)
    formatted_time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    print(f"[{formatted_time}]: Co-localization score matrix of workdir {curr_folder} has been saved!")

if __name__ == '__main__':
    
    arg = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="Generate co-localization score matrix.\nPlase cite: ", 
                        epilog="Co_localization Version: 1.0\n@ Cedric Liang (cedricliang@163.com)")
    arg.add_argument("-w", "--work_path", help="The path of the directory to process. Type = string", type=str)
    args = arg.parse_args()
    main(args)


