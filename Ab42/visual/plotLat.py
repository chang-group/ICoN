import numpy as np
import matplotlib.pyplot as plt
import sys


name_list = ['1z0q_1',      '1z0q_2',
             '1z0q_3',      '1z0q_4',
             '1z0q_5',      '1z0q_6',
             '2nao_mono_1', '2nao_mono_2',
             '2nao_mono_3', '2nao_mono_4']

Traj_ID = int(sys.argv[1])
path = './idps/'

dat = np.load(path + name_list[Traj_ID] + 'latent.npy')

fig = plt.figure(figsize = [8, 6])
ax = plt.axes(projection='3d')
p = ax.scatter3D(dat[:,0], dat[:,1], dat[:,2], c=range(int(len(dat))), alpha=0.9)
fig.colorbar(p, ax=ax)
ax.set_xlabel('Z1'); ax.set_ylabel('Z2'); ax.set_zlabel('Z3')
fig.tight_layout()
plt.savefig('figs/3D_latent.png')
plt.show()
