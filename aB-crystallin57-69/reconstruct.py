import os
import sys
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import random
import tqdm
import gc
import MDAnalysis as mda
from MDAnalysis import Universe
from MDAnalysis.analysis.bat import BAT

from src import BATT
from src import load_xyz
from src import fc_encode, fc_decode


def eval_model(dat, modelE, modelD, n_torsions, n_frames, device='cuda'):
       
    n1, n2, va, vb =  dat[0], dat[1], dat[2], dat[3]    
    batch_size = n1.shape[0]
    n_tors = n_torsions
    
    n1 = torch.reshape(n1, shape = (batch_size, n_tors*3, 1)).squeeze()
    n2 = torch.reshape(n2, shape = (batch_size, n_tors*3, 1)).squeeze()
    va = torch.reshape(va, shape = (batch_size, n_tors*3, 1)).squeeze()
    vb = torch.reshape(vb, shape = (batch_size, n_tors*3, 1)).squeeze()
       
    feats = torch.cat([n1, n2, va, vb], dim=-1)    
      
    modelE.eval()
    modelD.eval()

    z = modelE(feats)
    out = modelD(z)    

    feats_rmsd = (feats - out).pow(2).mean(axis=0).sqrt()
   
    n1 = out[:,:3*n_tors]
    n2 = out[:,3*n_tors:6*n_tors]
    va = out[:,6*n_tors:9*n_tors]
    vb = out[:,9*n_tors:]
    
    n1 = torch.reshape(n1.unsqueeze(-1), shape = (batch_size, n_tors, 3))
    n2 = torch.reshape(n2.unsqueeze(-1), shape = (batch_size, n_tors, 3))
    va = torch.reshape(va.unsqueeze(-1), shape = (batch_size, n_tors, 3))
    vb = torch.reshape(vb.unsqueeze(-1), shape = (batch_size, n_tors, 3))
    
    return z, n1, n2, va, vb, feats_rmsd


if __name__=='__main__':
    
    par_id = 5000
    device = 'cuda'
    dev_id = 0    
    
    eps = 1e-7  #toll
    out_path = 'output/'
    
    
    params_file_name = 'net_params'
    top_path = '../TRAJ/AB13.prmtop'
    path = '../TRAJ/AB13.dcd'

    xyz, (u, R), n_atoms = load_xyz(path, top_path)

    
    B = BATT(R, device=device)
    n_torsions = n_atoms - 3
    n_frames = xyz.shape[0]
    dat = B.Coords2MainVecs(xyz) # 4 vectors + bonds
    bonds = dat[4]
        
        
    n_feats = 4 * 3 * n_torsions    
    torch.cuda.set_device(dev_id)
    modelE = fc_encode(n_feats).to(device)
    modelD = fc_decode(n_feats).to(device)


    #--load pre-trained model
    checkpoint = torch.load(out_path + str(par_id) + params_file_name)
    modelE.load_state_dict(checkpoint['modelEncode'])
    modelD.load_state_dict(checkpoint['modelDecode'])

    z, n1, n2, va, vb, feats_rmsd = eval_model(dat, modelE, modelD, n_torsions, n_frames, device='cuda')

    ######################################################################
    z = z.detach().cpu().numpy()
    np.save('visual/' + 'latent' + '.npy', z)
    print('Pred rmsd--', feats_rmsd.mean().item())
    
    #set all root based bat to initial
    #root_based = B.GetRoot(xyz)
    #root_based[:,:] = root_based[0,:] 
    root_3_xyz = B.GetRootXYZ(xyz)
    coordsP = B.BatV2Coords(n1, n2, va, vb, bonds, root_3_xyz).detach().cpu().numpy()
        
    with mda.Writer('output/pdb/' + 'AB13' + '.dcd', n_atoms=u.atoms.n_atoms) as w:
        i = 0
        for ts in u.trajectory:
            ts.positions = coordsP[i,:,:]
            i = i + 1
            w.write(u.atoms)
            if i == n_frames:
                break
        
