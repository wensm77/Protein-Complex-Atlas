import os
import pandas as pd
import string
import random

def generate_unique_name(existing_names):
    while True:
        name = ''.join(random.choices(string.ascii_lowercase + string.digits, k=6))
        if name not in existing_names:
            return name

folder_path = '.../MPNN_data/foldseek/pdb'
files = os.listdir(folder_path)

mapping = {}

# 生成新文件名
for file in files:
    if file.endswith('.pdb'):
        new_name = generate_unique_name(mapping.values())
        mapping[file] = new_name + '.pdb'
        os.rename(os.path.join(folder_path, file), os.path.join(folder_path, new_name + '.pdb'))

# 保存映射到CSV文件
mapping_df = pd.DataFrame(list(mapping.items()), columns=['Original Name', 'New Name'])
mapping_df.to_csv('.../MPNN_data/foldseek/mapping.csv', index=False)
