import os
import sys
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import random
import tqdm
import MDAnalysis as mda
from MDAnalysis import Universe
from MDAnalysis.analysis.bat import BAT

from src import BATT, load_xyz
from src import fc_encode, fc_decode

def run_model(dat, modelE, modelD, indx, train = True, device='cuda'):
    
    n1, n2, va, vb =  dat[0][indx,:], dat[1][indx,:], dat[2][indx,:], dat[3][indx,:]
  
    batch_size = n1.shape[0]
    n_tors = n1.shape[1]
    
    
    n1_f = torch.reshape(n1, shape = (batch_size, n_tors*3, 1)).squeeze()
    n2_f = torch.reshape(n2, shape = (batch_size, n_tors*3, 1)).squeeze()
    va_f = torch.reshape(va, shape = (batch_size, n_tors*3, 1)).squeeze()
    vb_f = torch.reshape(vb, shape = (batch_size, n_tors*3, 1)).squeeze()
       
    feats = torch.cat([n1_f, n2_f, va_f, vb_f], dim=-1)
    del n1, n2
    
    if train:
        modelE.train();      modelD.train()    
        modelE.zero_grad();  modelD.zero_grad()
    else:
        modelE.eval()
        modelD.eval()

    z = modelE(feats)
    out = modelD(z)    
  
    criterion =  nn.SmoothL1Loss(reduction='mean', beta=1.0) 
    loss = criterion(out, feats) 
    
    if train:
        loss.backward()
        optimizer.step()
        #scheduler.step()
        
    return loss.item()



if __name__=='__main__':
        
    lrt = 0.0002 # learning rate
    #lrt = 0.0001
    max_epoch = 10000
    start = 0
    device='cuda'
    dev_id = 0    
    #ReTrain = True
    ReTrain = False
 
    batch_size = 200
    spacing ="\t"
    out_path = 'output/'
    
    name_list = ['idps/1z0q_1', 'idps/1z0q_2','idps/1z0q_3', 'idps/1z0q_4',
                 'idps/1z0q_5', 'idps/1z0q_6', 'idps/2nao_mono_1',
                 'idps/2nao_mono_2', 'idps/2nao_mono_3', 'idps/2nao_mono_4']

    Traj_ID = 9 #int(sys.argv[1]) # which traj to use from above list 
    params_file_name = 'net_params' + name_list[Traj_ID][5:]
    top_path = '../TRAJ/top.prmtop'
    
    print('Training --', name_list[Traj_ID])
    path = '../TRAJ/' + name_list[Traj_ID] + '.dcd'
    xyz, (u,R), n_atoms = load_xyz(path, top_path)

    T_id = torch.arange(0, xyz.shape[0], 2)
    V_id = torch.arange(1, xyz.shape[0], 2)

    xyzT = xyz[T_id,:,:]
    xyzV = xyz[V_id,:,:]
    
    n_train = xyzT.shape[0] 
    n_valid = xyzV.shape[0]

    n_torsions = n_atoms - 3
    n_feats = 4 * 3 * n_torsions # 4 vector feature   
    
    # get features
    B = BATT(R, device=device)
    datT = B.Coords2MainVecs(xyzT) # 4 vectors
    datV = B.Coords2MainVecs(xyzV)

    del xyz, xyzT, xyzV
    
    torch.cuda.set_device(dev_id)
    modelE = fc_encode(n_feats).to(device)
    modelD = fc_decode(n_feats).to(device)
      
    if ReTrain:
        checkpoint = torch.load(out_path + str(start) + params_file_name)
        modelE.load_state_dict(checkpoint['modelEncode'])
        modelD.load_state_dict(checkpoint['modelDecode'])
       
    optimizer = optim.Adam([{'params': modelE.parameters()},
                            {'params': modelD.parameters()}], lr = lrt)

    iStart = start + 1
    iEnd = max_epoch + 1
          
    
    for epoch in tqdm.tqdm(range(iStart, iEnd)):

        IT = np.random.randint(0, n_train, batch_size)
        lossT = run_model(datT, modelE, modelD, IT, train = True, device=device)
        
        IV = np.random.randint(0, n_valid, batch_size)
        lossV = run_model(datV, modelE, modelD, IV, train = False, device=device)
         
        
        if epoch % 1000 == 0:
            torch.save({'modelEncode': modelE.state_dict(), 'modelDecode': modelD.state_dict(),},
                        out_path + str(epoch) + params_file_name)

 
        #save loss   
        if epoch % 10 == 0:
            with open(out_path + 'log_train_val.txt', 'a') as fout:
                fout.write('%d\t%f\t%f\n'%(epoch, lossT, lossV))

            #save gradients to as text file
            ave_gradsE = []
            layersE = []
            for n, p in modelE.named_parameters():
                
                if(p.requires_grad) and ("bias" not in n) and (p.grad != None):
                    if epoch == 10:
                        layersE.append(n)#[15:25])
                    ave_gradsE.append(round(p.grad.abs().sum().item(),5))
                    
            if epoch == 10:
                line = spacing.join(char for char in layersE)
            else:
                line = spacing.join(str(e) for e in ave_gradsE)
            with open(out_path + 'gradsE.txt', 'a') as fg:
                fg.write(line + "\n") 

            ave_gradsD = []
            layersD = []
            for n, p in modelD.named_parameters():
                if(p.requires_grad) and ("bias" not in n) and (p.grad != None):

                    if epoch == 10:
                        layersD.append(n)#[15:25])
                    ave_gradsD.append(round(p.grad.abs().sum().item(),5))
                    
            if epoch == 10:
                line = spacing.join(char for char in layersD)
            else:
                line = spacing.join(str(e) for e in ave_gradsD)
            with open(out_path + 'gradsD.txt', 'a') as fg:
                fg.write(line + "\n")
 
