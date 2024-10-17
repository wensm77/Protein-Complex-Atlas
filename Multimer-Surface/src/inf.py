from transformers import EsmTokenizer, EsmForMaskedLM
from model import ResidueClassifierWithTransformer
import torch

def remove_module_prefix(state_dict):
    new_state_dict = {}
    for key, value in state_dict.items():
        new_key = key.replace('module.', '') if key.startswith('module.') else key
        new_state_dict[new_key] = value
    return new_state_dict


model = ResidueClassifierWithTransformer()
best_model_path = '/mnt/SaProt-main/mytrain/checkpoint_Transformer/best_model.pth'
state_dict = torch.load(best_model_path, map_location=torch.device('cpu'))
state_dict = remove_module_prefix(state_dict)
model.load_state_dict(state_dict)
model.cuda()

threshold = 0.369
seq = 'SdLdRdQdEdDfFdPwPkRdIwVpEdHfPfSaDaLdIeVeSaKfGfEaPkAdTkLtNaCtKaAiEdGgRpPvTfPwTdIkEwWkYaKfGvGnEhRtVdEdTaDcKvDnDpPvRpSfHqRwMhLqLdPpSrGrSmLiFiFrLrRgIaVhHdGdRpKpSdRhPpDpEfGgVkYmVwCmVwAtRgNdYpLhGgEiAdViShHdDiAyStLyEhVhAd'

model_path = "/mnt/SaProt-main/SaProt_650M_AF2"
tokenizer = EsmTokenizer.from_pretrained(model_path)
saprot = EsmForMaskedLM.from_pretrained(model_path)
device = "cuda"
saprot.to(device)
tokens = tokenizer.tokenize(seq)
inputs = tokenizer(seq, return_tensors="pt")

inputs = {k: v.to(device) for k, v in inputs.items()}
outputs = saprot(**inputs)
features = outputs.logits
features = features.cuda() 

logits = model(features)
probs = torch.sigmoid(logits)
predictions = (probs >= threshold).long().squeeze(0)
str_list = [str(int(i)) for i in predictions]
result_string = ''.join(str_list)[1:-1]
print(threshold)