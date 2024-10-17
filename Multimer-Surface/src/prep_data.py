import torch
from torch.utils.data import DataLoader, Dataset
from transformers import EsmTokenizer, EsmForMaskedLM
from tqdm import tqdm

# 初始化设备
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model_path = "/mnt/wensm/SaProt-main/SaProt_650M_AF2"
tokenizer = EsmTokenizer.from_pretrained(model_path)
saprot = EsmForMaskedLM.from_pretrained(model_path)

saprot.eval()

# 使用 DataParallel 来支持多张卡并行
saprot = torch.nn.DataParallel(saprot, device_ids=[0, 1, 2, 3])  # 指定卡的 ID
saprot.to(device)

# 创建原始数据集类，用于加载序列和标签
class RawProteinDataset(Dataset):
    def __init__(self, csv_file, tokenizer, max_length=1024):
        import pandas as pd
        self.data = pd.read_csv(csv_file)
        self.tokenizer = tokenizer
        self.max_length = max_length

        self.sequences = self.data['sequence'].tolist()
        self.labels = []
        for code_str in self.data['code'].astype(str).tolist():
            # 处理 labels，过滤非数字字符
            code_str = ''.join(filter(str.isdigit, code_str.strip()))
            label = [int(char) for char in code_str]
            self.labels.append(label)

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, idx):
        sequence = self.sequences[idx]
        labels = self.labels[idx]

        # Tokenize the sequence
        inputs = self.tokenizer(
            sequence,
            padding='max_length',
            max_length=self.max_length,
            return_tensors="pt",
        )

        inputs = {key: val.squeeze(0) for key, val in inputs.items()}
        if len(labels) > self.max_length:
            labels = labels[:self.max_length]
        else:
            labels = labels + [0] * (self.max_length - len(labels))
        labels = torch.tensor(labels, dtype=torch.float)
        return inputs, labels

max_length = 1024
batch_size = 128

train_dataset = RawProteinDataset('/mnt/wensm/SaProt-main/saport_data/test.csv', tokenizer, max_length=max_length)
train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=False)

train_data_save_path = 'test_data.pt'
all_data = []  # 用于保存每个样本的特征、mask 和标签

with torch.no_grad():
    for batch in tqdm(train_loader, desc='Extracting features'):
        inputs, labels = batch
        inputs = {k: v.to(device) for k, v in inputs.items()}
        
        outputs = saprot(**inputs)
        features = outputs.logits.cpu()  # (batch_size, sequence_length, hidden_size)
        masks = inputs['attention_mask'].cpu()      # (batch_size, sequence_length)

        # 将每个样本的数据保存为一个字典
        for i in range(features.size(0)):
            sample_data = {
                'feature': features[i],  # (sequence_length, hidden_size)
                'mask': masks[i],        # (sequence_length)
                'label': labels[i]       # (sequence_length)
            }
            all_data.append(sample_data)

# 保存所有数据到一个文件
torch.save(all_data, train_data_save_path)
