#!/usr/bin/env python

import pytraj as pt
import numpy as np

name = '1z0q_1'
TOP = '../../../top.prmtop'
TRAJ = '../step2/' + name + 'SyntheticE.dcd'
traj = pt.iterload(TRAJ, TOP, stride=1)


def SortFrames(traj, cut=1.5, freq=50, mask = ['!@H*']):
    """
    
    """
   
    ind1 = []
    nf = traj.shape[0]
    bol = np.ones(nf)
    
    for i in range(0, nf, freq):
        print('Processing conformation -- ', i)
        rmsd = pt.rmsd(traj, ref=i, mask = mask) #uses i_th frame as ref
       
        bol[rmsd < cut] = 0
        ind1.append(i)
    
    ind1 = np.array(ind1)
    ind2 = np.where(bol > 0)[0]
    ind = np.concatenate([ind1, ind2])
    ind.sort()
    
    return traj[ind], ind

traj_new , ind = SortFrames(traj, cut=1.0, freq=50, mask=['!@H*'])

#output trajs
pt.write_traj('SyntheticSR_Removed.dcd', traj_new, overwrite=True)

