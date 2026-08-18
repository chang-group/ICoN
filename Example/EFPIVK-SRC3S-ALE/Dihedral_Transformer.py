import torch
import torch.nn as nn

class ICONTransformer(nn.Module):
    def __init__(self, n_features, n_head=8, n_layers=6, d_model=512, d_latent=3):
        super(ICONTransformer, self).__init__()
        self.d_model = d_model
        self.n_features = n_features

        # Input embedding
        self.input_embedding = nn.Linear(2, d_model)  # 4 features per token: bond, angle, sin, cos

        # Positional encoding
        self.pos_encoder = PositionalEncoding(d_model)

        # Transformer Encoder
        encoder_layers = nn.TransformerEncoderLayer(d_model=d_model, nhead=n_head, activation='gelu')
        self.transformer_encoder = nn.TransformerEncoder(encoder_layers, num_layers=n_layers)

        # Projection to latent space
        self.to_latent = nn.Linear(d_model * n_features // 2, d_latent)

        # Projection from latent space
        self.from_latent = nn.Linear(d_latent, d_model * n_features // 2)

        # Transformer Decoder
        decoder_layers = nn.TransformerEncoderLayer(d_model=d_model, nhead=n_head, activation='gelu')
        self.transformer_decoder = nn.TransformerEncoder(decoder_layers, num_layers=n_layers)

        # Output projection
        self.output_projection = nn.Linear(d_model, 2)

    def encode(self, src):
        # src shape: [batch_size, n_features, 4]
        src = src.transpose(0, 1)  # [n_features, batch_size, 4]
        src = self.input_embedding(src)  # [n_features, batch_size, d_model]
        src = self.pos_encoder(src)
        memory = self.transformer_encoder(src)
        memory = memory.transpose(0, 1)  # [batch_size, n_features, d_model]
        memory = memory.reshape(memory.size(0), -1)  # [batch_size, n_features * d_model]
        latent = self.to_latent(memory)  # [batch_size, d_latent]
        return latent

    def decode(self, latent):
        # latent shape: [batch_size, d_latent]
        hidden = self.from_latent(latent)  # [batch_size, n_features * d_model]
        hidden = hidden.reshape(hidden.size(0), self.n_features // 2, self.d_model)
        hidden = hidden.transpose(0, 1)  # [n_features, batch_size, d_model]
        hidden = self.pos_encoder(hidden)
        output = self.transformer_decoder(hidden)
        output = output.transpose(0, 1)  # [batch_size, n_features, d_model]
        output = self.output_projection(output)  # [batch_size, n_features, 4]
        return output

    def forward(self, src):
        latent = self.encode(src)
        output = self.decode(latent)
        return output, latent

class PositionalEncoding(nn.Module):
    def __init__(self, d_model, max_len=5000):
        super(PositionalEncoding, self).__init__()
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * (-torch.log(torch.tensor(10000.0)) / d_model))
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0).transpose(0, 1)
        self.register_buffer('pe', pe)

    def forward(self, x):
        x = x + self.pe[:x.size(0), :]
        return x

def init_weights(m):
    if type(m) == nn.Linear:
        torch.nn.init.xavier_uniform_(m.weight)
        if m.bias is not None:
            torch.nn.init.zeros_(m.bias)

if __name__ == '__main__':
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(device)

    batch_size = 1
    n_torsions = 178
    n_feats = 2 * n_torsions

    inp = torch.rand(batch_size, n_torsions, 2).to(device)

    model = ICONTransformer(n_torsions).to(device)
    model.apply(init_weights)

    out, latent = model(inp)

    print('Number of features:', n_feats)
    print('Input shape:', inp.shape)
    print('Latent shape:', latent.shape)
    print('Output shape:', out.shape)
    print('Number of Model parameters:', sum(p.numel() for p in model.parameters() if p.requires_grad))
