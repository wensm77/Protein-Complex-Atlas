import os
import argparse
    
def filter_indentity_args(args):
    file_path = args.file_path
    min_seq_id = args.min_seq_id
    max_seq_id = args.max_seq_id
    outname_base = args.outname
    save_path = args.save_path

    if not os.path.exists(save_path):
        os.makedirs(save_path)
    
    file_extention = os.path.splitext(file_path)[1]
    out_name = f"{outname_base}_{min_seq_id}{file_extention}"
    output_file_path = os.path.join(save_path, out_name)
    cnt = 0
    cnt_all = 0
    search_dict = {}

    with open(file_path, "r") as input_file:
        with open(output_file_path, "w") as output_file:
            for line in input_file:
                cnt_all += 1
                columns = line.strip().split('\t')
                if float(columns[2]) >= min_seq_id and float(columns[2]) <= max_seq_id:
                    output_file.write(line)
                    cnt += 1
                    
    print(f"Filter with minimum sequence identity {min_seq_id} and maximum sequence identity {max_seq_id} is completed!!!")
    print(f"The result {file_extention} file contains {cnt_all} sequences! Sequences with identity between {min_seq_id} and {max_seq_id} are totally {cnt}!")

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="Filter sequences from searched results with sequence identity between -i and -a.\nPlase cite: ", 
                        epilog="Co_localization Version: 1.0\n@ Cedric Liang (cedricliang@163.com)")
    parser.add_argument("-p","--file_path", type=str, required=True, help="The path to the input file. e.g., *.m8. Type = string")
    parser.add_argument("-i","--min_seq_id", type=float, default=0.3, help="Minimum sequence identity threshold. Type = float [default = 0.3]")
    parser.add_argument("-a","--max_seq_id", type=float, default=1.0, help="Maximum sequence identity threshold. Type = float [default = 1.0]")
    parser.add_argument("-o","--outname", type=str, default="output", help="Prefix for the output file. Type = string [default = 'output']")
    parser.add_argument("-s", "--save_path", type=str, default="./", help="Directory to save the output file. Type = string [default = pwd]")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    filter_indentity_args(args)
