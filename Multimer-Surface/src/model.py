import torch.nn as nn

class ResidueClassifier(nn.Module):
    def __init__(self, input_dim=446, hidden_dim=2048):
        super(ResidueClassifier, self).__init__()
        self.classifier = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 1)
        )
        self._initialize_weights()
    def forward(self, x):
        logits = self.classifier(x).squeeze(-1)  # (B, L=1024)
        return logits

    def _initialize_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.kaiming_uniform_(m.weight, nonlinearity='relu')
                if m.bias is not None:
                    nn.init.zeros_(m.bias)
    


class ResidueClassifierWithTransformer(nn.Module):
    def __init__(self, origin_dim=446, input_dim=1024, hidden_dim=2048, num_heads=8, num_layers=4):
        super(ResidueClassifierWithTransformer, self).__init__()
        
        self.adjust = nn.Linear(origin_dim, input_dim)
        encoder_layer = nn.TransformerEncoderLayer(d_model=input_dim, nhead=num_heads, dim_feedforward=hidden_dim)
        self.transformer_encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)
        
        self.classifier = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 1)  
        )
        
        self._initialize_weights()

    def forward(self, x):
        x = self.adjust(x) # [B, L, 1024]
        x = self.transformer_encoder(x) 
        logits = self.classifier(x).squeeze(-1)  
        return logits 

    def _initialize_weights(self):
        for m in self.classifier.modules():
            if isinstance(m, nn.Linear):
                nn.init.kaiming_uniform_(m.weight, nonlinearity='relu')
                if m.bias is not None:
                    nn.init.zeros_(m.bias)

