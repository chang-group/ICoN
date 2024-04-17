import numpy as np
import matplotlib.pyplot as plt
import sys

dat = np.load('out/latent.npy')
n = dat.shape[0]
fig = plt.figure(figsize = [8, 6])

ax = plt.axes(projection='3d')
p = ax.scatter3D(dat[:,0], dat[:,1], dat[:,2], c=range(int(len(dat))), alpha=0.9)
fig.colorbar(p, ax=ax)
ax.set_xlabel('Z1'); ax.set_ylabel('Z2'); ax.set_zlabel('Z3')
fig.tight_layout()
plt.show()
