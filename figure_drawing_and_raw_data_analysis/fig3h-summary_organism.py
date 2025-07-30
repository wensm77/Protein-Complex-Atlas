import pandas as pd
from collections import Counter

df = pd.read_csv('0724cluster_with_seq_org.csv')
# 每个类（pdb_id_1）下的organism_2种类集合
organism_sets = df.groupby('pdb_id_1')['organism_2'].apply(lambda x: set(x.dropna()))
# 统计每种分布情况出现的类数
summary = Counter([','.join(sorted(list(s))) for s in organism_sets])
# 输出统计结果
pd.DataFrame(summary.items(), columns=['organism_2_types', 'class_count']).to_csv('0724cluster_class_organism_summary.csv', index=False)
print('统计完成，结果已写入0724cluster_class_organism_summary.csv') 