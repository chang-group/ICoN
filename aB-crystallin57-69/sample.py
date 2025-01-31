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


def eval_model(z, modelD, n_tors, n_frames, device='cuda'):
       
    
    batch_size = n_frames
    modelD.eval()
    out = modelD(z)    

    n1 = out[:,:3*n_tors]
    n2 = out[:,3*n_tors:6*n_tors]
    va = out[:,6*n_tors:9*n_tors]
    vb = out[:,9*n_tors:]
    
    n1 = torch.reshape(n1.unsqueeze(-1), shape = (batch_size, n_tors, 3))
    n2 = torch.reshape(n2.unsqueeze(-1), shape = (batch_size, n_tors, 3))
    va = torch.reshape(va.unsqueeze(-1), shape = (batch_size, n_tors, 3))
    vb = torch.reshape(vb.unsqueeze(-1), shape = (batch_size, n_tors, 3))
    
    return n1, n2, va, vb



def PathSample(z1, z2, n_steps, seed, amp):
    """
    
    """
    
    np.random.seed(seed=seed)
    Center = (z1 + z2) /2.0
    Center = Center - amp*(np.random.rand(3)-0.5)
    
    R1 = np.sqrt(((z1 - Center)**2).sum())
    R2 = np.sqrt(((z2 - Center)**2).sum())
    dR = (R2 - R1)/n_steps
    
    
    theta1 = np.arccos((z1[2] - Center[2]) / R1)
    phi1 = np.arctan2(z1[1]-Center[1], z1[0]-Center[0])
    
    
    theta2 = np.arccos((z2[2] -Center[2]) / R2)
    phi2 = np.arctan2(z2[1]-Center[1],z2[0]-Center[0])


    dphi = (phi2 - phi1)/n_steps    
    dtheta = (theta2 - theta1)/n_steps    


    pos = []
    phi = phi1
    theta = theta1
    R = R1
    for i in range(n_steps):
    
        phi = phi + dphi         
        R = R + dR
        theta = theta + dtheta
        
        path = Center + np.array([R * np.sin(theta) * np.cos(phi),
                               R * np.sin(theta) * np.sin(phi),
                               R * np.cos(theta)])
        pos.append(path)

    pos = np.array(pos)
    noise = 0.2*(np.random.rand(n_steps,3)-0.5)
    pos[1:-1,:] = pos[1:-1,:] + noise[1:-1,:]  
    return pos


def LinInter(z1, z2, n):
    """
    Linear interpolation per small interval
    """
    dz = (z2 - z1)/n
    z_out = np.zeros(shape=(n, 3))
    
    for i in range(n):
        z_out[i,:] = z1 + i*dz
    
    return z_out



if __name__=='__main__':
  
    
    par_id = 5000
    device = 'cuda'
    dev_id = 0    
     
    eps = 1e-7  #toll
    out_path = 'output/'
     
     
    params_file_name = 'net_params'
    top_path = '../TRAJ/AB13.prmtop'
    path = '../TRAJ/0-frame.pdb'
    out_path = 'output/'
    pdb_name = "AB13"

    #id_list = [50, 56]

    xyz, (u, R), n_atoms = load_xyz(path, top_path, isBAT=False)
    B = BATT(R, device=device)
    n_torsions = n_atoms - 3
    n_feats = 4 * 3 * n_torsions    
    dat = B.Coords2MainVecs(xyz) # 4 vectors + bonds
    bonds = dat[4]

    #model
    torch.cuda.set_device(dev_id)
    modelD = fc_decode(n_feats).to(device)

    #--load pre-trained model
    checkpoint = torch.load(out_path + str(par_id) + params_file_name)
    modelD.load_state_dict(checkpoint['modelDecode'])
   
    
    #first_min, sec_min = 16627, 12917
    z = np.load('visual/' +  'latent.npy')
    #z = torch.as_tensor(z, device=device, dtype=torch.float32)
    
    n_samples = 10 #per interval
    seed = 10
    amp = .5

    tmp = np.zeros((1,3))
    n_tot = z.shape[0]
    for i in range(n_tot-1):
        
        z1 = z[i] #1st point
        z2 = z[i+1] # 2nd point
            
        z_int = PathSample(z1, z2, n_samples, seed, amp)
        tmp = np.concatenate([tmp, z_int])
        
    z_new = torch.as_tensor(tmp[1:,:], device=device, dtype=torch.float32)
    n_frames = n_samples * (n_tot-1) 
 
    #evaluate with decoder
    n1, n2, va, vb = eval_model(z_new, modelD, n_torsions, n_frames, device=device)
    bondsP = torch.zeros(n_frames, n_torsions, device=device, dtype=torch.float32)
    bondsP[:,:] = bonds
    
    #get root based atom positions from ref
    root_based_ref = B.GetRootXYZ(xyz)
    p1 = torch.zeros(n_frames, 3, device=device, dtype=torch.float32)
    p2 = torch.zeros(n_frames, 3, device=device, dtype=torch.float32)
    p3 = torch.zeros(n_frames, 3, device=device, dtype=torch.float32)
    p1[:,:] = root_based_ref[0];  p2[:,:] = root_based_ref[1];  p3[:,:] = root_based_ref[2]
    root_3_xyz = (p1, p2, p3)

    coordsP = B.BatV2Coords(n1, n2, va, vb, bondsP, root_3_xyz).detach().cpu().numpy()
    
        
    u.trajectory.n_frames=n_frames
    with mda.Writer('output/pdb/idps/' + 'NonLinIntAll.dcd',
                    n_atoms = u.atoms.n_atoms) as w:

        for i in range(n_frames):
            for ts in u.trajectory:
                ts.positions = coordsP[i,:,:]
                w.write(u.atoms)
       
    
