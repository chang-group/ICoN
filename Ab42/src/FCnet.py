import os
import sys
import torch
from torch import nn
from torch.nn.modules.module import Module
import numpy as np


def init_weights(m):
    """
    
    """

    if type(m) == nn.Conv3d:
        torch.nn.init.xavier_uniform_(m.weight)
    if type(m) == nn.Linear:
        torch.nn.init.xavier_uniform_(m.weight)


class fc_encode(Module):
    def __init__(self, n_feats):
        super(fc_encode, self).__init__()
        bias = False
        self.fc = nn.Sequential(
            nn.Linear(n_feats, n_feats//4, bias=bias),
            nn.LayerNorm(n_feats//4),
            nn.LeakyReLU(),
            nn.Linear(n_feats//4, n_feats//8, bias=bias),
            nn.LeakyReLU(),
            nn.Dropout(p=0.1),
            nn.Linear(n_feats//8, n_feats//16, bias=bias),
            nn.LeakyReLU(),
            nn.Linear(n_feats//16, n_feats//32, bias=bias),
            nn.LeakyReLU(),
            nn.Linear(n_feats//32, n_feats//64, bias=bias),
            nn.LeakyReLU(),
            nn.Linear(n_feats//64, n_feats//64, bias=bias),
            nn.LeakyReLU(),
            nn.Linear(n_feats//64, 3, bias=bias),
        )
        

        self.fc.apply(init_weights)
    
    def forward(self, input):
        return self.fc(input)
    
    
class fc_decode(Module):
    def __init__(self, n_feats):
        super(fc_decode, self).__init__()
        bias = False
        self.fc = nn.Sequential(
            nn.Linear(3, n_feats//64, bias=bias),
            nn.LeakyReLU(),
            nn.Linear(n_feats//64, n_feats//64, bias=bias),
            nn.LeakyReLU(),
            nn.Linear(n_feats//64, n_feats//32, bias=bias),
            nn.LeakyReLU(),
            nn.Linear(n_feats//32, n_feats//16, bias=bias),
            nn.LeakyReLU(),
            #nn.Dropout(p=0.1),
            nn.Linear(n_feats//16, n_feats//8, bias=bias),
            nn.LeakyReLU(),
            nn.Linear(n_feats//8, n_feats//4, bias=bias),
            nn.LeakyReLU(),
            nn.Linear(n_feats//4, n_feats, bias=bias),
        )

        self.fc.apply(init_weights)
    
    def forward(self, input):
        return self.fc(input)
    
    
if __name__=='__main__':

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(device)
    
    batch_size = 2
    n_torsions = 624
    n_feats =  4 * 3 * n_torsions
    
    inp = torch.rand(batch_size, n_feats).to(device)    
    
    modelE = fc_encode(n_feats).to(device).eval()
    modelD = fc_decode(n_feats).to(device).eval()
    
    n_params = sum(p.numel() for p in modelE.parameters() if p.requires_grad) + \
               sum(p.numel() for p in modelD.parameters() if p.requires_grad)
    
    latent = modelE(inp)  
    out = modelD(latent) 
    
    print('number of features -- ', n_feats)
    print('input dim', inp.shape)
    print('Latent dim', latent.shape)
    print('Output dim', out.shape)
    print('Number of Model parameters -- ', n_params)
