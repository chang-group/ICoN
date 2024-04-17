import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.spatial.distance import cdist

name_list = ['1z0q_1',      '1z0q_2',
             '1z0q_3',      '1z0q_4',
             '1z0q_5',      '1z0q_6',
             '2nao_mono_1', '2nao_mono_2',
             '2nao_mono_3', '2nao_mono_4']

Traj_ID1 = 9#int(sys.argv[1])
Traj_ID2 = 7#int(sys.argv[2])

path = './idps/'

dat1 = np.load(path + name_list[Traj_ID1] + 'latent.npy')
dat2 = np.load(path + name_list[Traj_ID2] + 'latent.npy')

print(dat1.shape, dat2.shape)
dist = cdist(dat1, dat2)
print(dist.shape)

plt.show()
