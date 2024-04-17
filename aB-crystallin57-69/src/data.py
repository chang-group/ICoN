import numpy as np
import MDAnalysis as mda
from MDAnalysis import Universe
from MDAnalysis.analysis.bat import BAT
import torch


def load_xyz(path, top_path, device="cuda", isBAT=False):
    """

    """

    u = mda.Universe(top_path, path)
    sel = u.select_atoms("protein")
    n_atoms = sel.n_atoms
    R = BAT(sel)

    xyz = []
    for i in R._trajectory:
        xyz.append(i.positions)
    xyz = np.array(xyz)
    xyz = torch.as_tensor(xyz, device=device, dtype=torch.float32)

    if isBAT:
        R.run()
        bat = R.results.bat
        return xyz, (u, R), n_atoms, bat
    else:
        return xyz, (u, R), n_atoms



if __name__=='__main__':

    data_path = '../../../data/14sbff_alpha/sorted/train/14sbff_alpha_train.dcd'
    prmtop_path = '../../../data/14sbff_alpha/14sbff_alpha.prmtop'

    xyz, _, _ = load_xyz(data_path, prmtop_path)