import os
import csv
import matplotlib.pyplot as plt

class Union_Bac():
    def __init__(self, deal_list:[], total_size:int, output_folder:str):
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        self.output_folder = output_folder
        self.total_size = total_size
        self.deal_list = deal_list
        self.bacteria_names = ["" for _ in range(len(self.deal_list))]
        tmp = [faa.split(".")[0] for faa in os.listdir("/home/madacheng/DATA1/NCBI/bac_36")]
        for bac_name in tmp:
            idx = int(bac_name.split("_")[0])
            self.bacteria_names[self.deal_list.index(idx)] = bac_name
        self.pair_dict_list = [{} for _ in range(len(self.deal_list))]

    #  Pair Index	Sequence 1	Sequence 2	SequenceName 1	SequenceName 2	Min Distance	Ad Points(row, col)	Ad Points(col, row)	ClusterSize1	ClusterSize2
    def add_pairs2dict(self, csv_file:str, dict_idx:int): 
        if dict_idx < 0:
            print("Such bacteria does not exist!")
            return     
        with open(csv_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader)
            for row in reader:
                tpl = tuple(sorted([row[3], row[4]]))
                if tpl not in self.pair_dict_list[dict_idx]:
                    self.pair_dict_list[dict_idx][tpl] = row[1:6]

    def add_pairs2dict_batch(self, bac_name:str, pair_index:int = 1) -> int:
        base_folder = "/home/madacheng/liangjq/colocalization_All"
        dict_idx = self.bacteria_names.index(bac_name)
        job_folder = f"Job_{bac_name}"
        csv_name = f"output_seq_pairs_{bac_name}.csv"
        for i in range(1, self.total_size + 1):
            colocal_folder = f"colocalization_{i}"
            csv_file = os.path.join(base_folder, colocal_folder, job_folder, csv_name)
            self.add_pairs2dict(csv_file, dict_idx)

        output_file = os.path.join(self.output_folder, f"{bac_name}_pair_dict.tsv")
        with open(output_file, 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(["Pair Index", "Sequence 1", "Sequence 2", "SequenceName 1", "SequenceName 2", "Min Distance"])
            for val in self.pair_dict_list[dict_idx].values():
                tmp_list = [pair_index] + val
                writer.writerow(tmp_list)
                pair_index += 1
        
        print(f"{bac_name} done! Totally pairs: {len(self.pair_dict_list[dict_idx])}")
        return pair_index

    def add_pairs2dict_batch_all(self, pair_index:int = 1) -> int:
        start_index = pair_index
        for bac_name in self.bacteria_names:
            pair_index = self.add_pairs2dict_batch(bac_name, pair_index)
        print(f"All done! Totally pairs: {pair_index - start_index}")
        return pair_index


if __name__ == "__main__":
    to_be_dealt = [1,3,4,6,13,5,7,8,9,12,15,17,33,2,11,36]
    another = range(80,100)
    to_be_dealt = to_be_dealt + list(another)

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, 
                        description="Union 13 split tsv.\nPlase cite: ", 
                        epilog="Co_localization Version: 1.0\n@ Cedric Liang (cedricliang@163.com)")
    parser.add_argument("-o", "--output_dir", type=str, required=True,
                        help="The directory of union tsv output. Type = string")
    args = parser.parse_args()

    output_path = args.output_dir
    union_bac = Union_Bac(to_be_dealt, 13, output_path)

    union_bac.add_pairs2dict_batch_all()
