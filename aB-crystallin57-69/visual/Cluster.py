import numpy as np
import matplotlib.pyplot as plt
import sys
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering

path = './out/'
dat = np.load(path + 'latent.npy')

cmap='tab20b'
n_clus = 20
cb_labs = list(range(n_clus))

#kmeans = KMeans(n_clusters=20, random_state=10, n_init=7).fit(dat)
#labels = kmeans.labels_
agg = AgglomerativeClustering(n_clusters=n_clus).fit(dat)
labels = agg.labels_



fig = plt.figure(figsize = [8, 6])
ax = plt.axes(projection='3d')
for i in cb_labs:
    line = dat[labels==i].mean(axis=0) 
    ax.text(line[0], line[1], line[2], s = str(i))

p = ax.scatter3D(dat[:,0], dat[:,1], dat[:,2], c=labels, alpha=0.8, cmap=cmap)
ax.set_xlabel('Z1'); ax.set_ylabel('Z2'); ax.set_zlabel('Z3')

cbar = fig.colorbar(p, ax=ax)
cbar.ax.set_yticklabels("")

fig.tight_layout()
plt.savefig('figs/Cluster_lat' +  '.png', dpi=800)
plt.show()
