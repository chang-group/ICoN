#!/usr/bin/env python
# coding: utf-8

import pytraj as pt
import numpy as np
import matplotlib.pyplot as plt


#mask = ['@N,@CA,@C,@O']
mask = ['!@H*']

Traj_ID = 9
name_list = ['1z0q_1',      '1z0q_2',
             '1z0q_3',      '1z0q_4',
             '1z0q_5',      '1z0q_6',
             '2nao_mono_1', '2nao_mono_2',
             '2nao_mono_3', '2nao_mono_4']

#path to MD
path = '../../../TRAJ/idps/'
TOP = './top.prmtop'
TRAJ = path + name_list[Traj_ID]  + '.dcd'
traj = pt.iterload(TRAJ, TOP, stride=1)

#path to DL reconstruction
path_r = './idps/'
TRAJ_r = path_r + name_list[Traj_ID] + '.dcd'
traj_r = pt.iterload(TRAJ_r, TOP, stride=1)



n_frames = traj.shape[0]
pair_rmsd = []
for i in range(0,n_frames,1):
    p_rmsd = pt.rmsd(traj[i:i+1], ref=traj_r[i:i+1], mask = mask)
    pair_rmsd.append(p_rmsd)
    
pair_rmsd = np.array(pair_rmsd)
pair_rmsd = pair_rmsd[:,0]

#split train and val 
pair_rmsdV = pair_rmsd[np.arange(1,10000,2)]
pair_rmsdT = pair_rmsd[np.arange(0,10000,2)]


print('Validatio pair rmsd mean -- ', pair_rmsdV.mean())



tcl = 'yellow'
t_size=16
bg_cl = '#1b212c'#'#1A1A1A'#'black'
fig = plt.figure(figsize=(10,6), facecolor=bg_cl)
plt.rcParams['axes.facecolor'] = bg_cl
plt.hist(pair_rmsdT, bins=50, density=True, color='cyan')
plt.hist(pair_rmsdV, bins=50, density=True, color='orange')

plt.title('pair rmsd distribution' , color=tcl, size=t_size)

plt.xticks(color=tcl,size=15)
plt.yticks(color=tcl,size=15)
plt.xlabel('pair rmsd', color=tcl, size=t_size)
plt.ylabel('P', color=tcl, size=t_size)

plt.show()

