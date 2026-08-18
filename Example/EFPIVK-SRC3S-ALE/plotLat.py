import numpy as np
import matplotlib.pyplot as plt
import matplotlib

print(matplotlib.__version__)  # Make sure the Matplotlib version is at least 3.5
import sys
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

from sympy.physics.quantum.circuitplot import matplotlib

matplotlib.use('TkAgg')
from sklearn.manifold import TSNE

dat = np.load('../CyclicPeptideOutput/visual/Driver_Cyclic_Peptide_Predictionlatent.npy')
dat_PCA = np.load('../CyclicPeptideOutput/visual/PC1_Driver_Cyclic_Peptide_Predictionlatent.npy')
dat_min = np.load('../CyclicPeptideOutput/visual/Driver_Cyclic_Peptide_Prediction_Minimizedlatent.npy')
dat_inter = np.load('../CyclicPeptideOutput/visual/interpolated_latent_points_slerp.npy')
#dat_inter = np.load('../CyclicPeptideOutput/visual/spherical_expansion_latent.npy')
centers = np.load('../CyclicPeptideOutput/visual/cluster_centers.npy')

sub_id = np.arange(0, dat.shape[0])
dat = dat[sub_id, :]

fig = plt.figure(figsize=[8, 6])
ax = plt.axes(projection='3d')
indices = np.arange(len(dat))
p = ax.scatter3D(dat[:, 0], dat[:, 1], dat[:, 2], c=indices, alpha=0.1, picker=5)
pInt = ax.scatter3D(dat_inter[:, 0], dat_inter[:, 1], dat_inter[:, 2], c='red', alpha=1, picker=5)
#pMin = ax.scatter3D(dat_min[:, 0], dat_min[:, 1], dat_min[:, 2], c='orange', alpha=0.9, picker=5)
# p2 = ax.scatter3D(dat_PCA[:, 0], dat_PCA[:, 1], dat_PCA[:, 2], c='red', alpha=0.9)

# Plot centers
pCenters = ax.scatter3D(centers[:, 0], centers[:, 1], centers[:, 2], c='green', s=100, alpha=1)

for i, center in enumerate(centers):
    ax.text(center[0], center[1], center[2], str(i), fontsize=12, fontweight='bold')

fig.colorbar(p, ax=ax, label='Index')
ax.set_xlabel('Z1');
ax.set_ylabel('Z2');
ax.set_zlabel('Z3')


def onpick(event):
    ind = event.ind
    if len(ind) > 0:
        closest_point_index = ind[0]
        print(f"Clicked near point with index: {closest_point_index}")
        print(f"Data for this point: {dat[closest_point_index]}")


fig.canvas.mpl_connect('pick_event', onpick)
fig.tight_layout()
plt.savefig('3D_latent_with_centers.png')
plt.show()
