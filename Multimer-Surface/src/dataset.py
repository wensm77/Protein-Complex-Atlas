import torch
from torch.utils.data import Dataset

class PrecomputedProteinDataset(Dataset):
    def __init__(self, data_file):
        self.data = torch.load(data_file)  # 加载所有数据，data 是一个列表，每个元素是一个字典

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        sample_data = self.data[idx]
        feature = sample_data['feature']  # (sequence_length, hidden_size)
        mask = sample_data['mask']        # (sequence_length)
        label = sample_data['label']      # (sequence_length)
        return feature, mask, label
