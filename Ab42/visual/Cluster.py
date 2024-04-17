import numpy as np
import matplotlib.pyplot as plt
import sys
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering



name_list = ['1z0q_1', '1z0q_2','1z0q_3', '1z0q_4',
             '1z0q_5', '1z0q_6', '2nao_mono_1',
             '2nao_mono_2', '2nao_mono_3', '2nao_mono_4']

Traj_ID = int(sys.argv[1])
path = './idps/'

dat = np.load(path + name_list[Traj_ID] + 'latent.npy')
cmap='tab20b'

n_clus = 20
cb_labs = list(range(n_clus))

#kmeans = KMeans(n_clusters=20, random_state=10, n_init=7).fit(dat)
#labels = kmeans.labels_
agg = AgglomerativeClustering(n_clusters=n_clus).fit(dat)
labels = agg.labels_



fig = plt.figure(figsize = [10, 8])
ax = plt.axes(projection='3d')

for i in cb_labs:
    line = dat[labels==i].mean(axis=0) 
    ax.text(line[0], line[1], line[2], s = str(i))


p = ax.scatter3D(dat[:,0], dat[:,1], dat[:,2], c=labels, alpha=0.8, cmap=cmap)

#ax.set_xlabel('Z1'); ax.set_ylabel('Z2'); ax.set_zlabel('Z3')
#cbar = fig.colorbar(p, ax=ax)
#cbar.ax.set_yticklabels("")

ax.set_axis_off()

fig.tight_layout()
plt.savefig('figs/Cluster_lat' + name_list[Traj_ID] +  '.png', dpi=800)
plt.show()
