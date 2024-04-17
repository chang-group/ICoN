import MDAnalysis as mda
import numpy as np
from MDAnalysis import Universe
from MDAnalysis.analysis.bat import BAT
import matplotlib.pyplot as plt
import torch
import time


class BATT:
   
    def __init__(self, R, device='cuda'):
        """
            
        """
        self.R = R
        self.tor_indxs = torch.tensor(R._torsion_XYZ_inds,
                                      dtype=torch.long,
                                      device=device)
        self.device = device        
        self.n_tors = R.atoms.n_atoms - 3


    
    def GetRootXYZ(self, xyz):
        """
        
        """
        
        n_frames = xyz.shape[0]
          
        p0 = xyz[:,self.R._root_XYZ_inds[0],:] 
        p1 = xyz[:,self.R._root_XYZ_inds[1],:] 
        p2 = xyz[:,self.R._root_XYZ_inds[2],:] 
        
        return p0, p1, p2
    

    def Coords2MainVecs(self, xyz):
        """
        is differentiable -- Yes
        
        """       
        
        #get all positions for torsion calc
        p1 = xyz[:,self.tor_indxs[:,0]]
        p2 = xyz[:,self.tor_indxs[:,1]]
        p3 = xyz[:,self.tor_indxs[:,2]]
        p4 = xyz[:,self.tor_indxs[:,3]]
        
        # get all bond vecs for all bat
        va = p2 - p1
        vb = p3 - p2 # middle bond vector
        vc = p4 - p3

        #get bond length and store
        bonds = va.pow(2).sum(dim=-1).sqrt()
        
        #keep all as unit vectors
        va = va / (va*va).sum(dim=-1, keepdims=True).sqrt()
        vb = vb / (vb*vb).sum(dim=-1, keepdims=True).sqrt()
        vc = vc / (vc*vc).sum(dim=-1, keepdims=True).sqrt()


        n1 = torch.cross(-va, vb)  # n1 is normal vector to -va, vb
        n2 = torch.cross(-vb, vc)  # n2 is normal vector to -vb, vc

        n1 = n1 / (n1*n1).sum(dim=-1, keepdims=True).sqrt()
        n2 = n2 / (n2*n2).sum(dim=-1, keepdims=True).sqrt()


        return n1, n2, va, vb, bonds


    def AngleTorsion(self, n1, n2, va, vb):
        """
        is differentiable -- Yes
        
        """
        
        
        #get bond angles
        x1 = (-va*vb).sum(dim=-1)
        y1 = torch.cross(va, vb).pow(2).sum(axis=-1).sqrt()
        angles = torch.atan2(y1,x1)
        
        #get torsions
        xp = torch.cross(n1, n2)
        x2 = (n1 * n2).sum(dim=-1)
        y2  = (xp * vb).sum(dim=-1) #/ (vb*vb).sum(dim=-1).sqrt()
        torsions = torch.atan2(y2,x2)
    
        return angles, torsions


    def BatV2Coords(self, n1, n2, va, vb, bonds, root_3_xyz):
        """
    
        """          
        
        n_frames = n1.shape[0]
        p0, p1, p2 = root_3_xyz

        XYZ = torch.zeros(n_frames, self.n_tors + 3, 3, device=self.device)
        XYZ[:,self.R._root_XYZ_inds[0],:] = p0
        XYZ[:,self.R._root_XYZ_inds[1],:] = p1
        XYZ[:,self.R._root_XYZ_inds[2],:] = p2
                 
        angles, torsions = self.AngleTorsion(n1, n2, va, vb)
        angles = angles.clamp(min=1.52, max=2.47) #??????????????????
        
        cs_ang, sn_ang  = torch.cos(angles), torch.sin(angles)
        cs_tor, sn_tor  = torch.cos(torsions), torch.sin(torsions)
        
        i = 0
        for (a0,a1,a2,a3) in self.tor_indxs:
            
            pos1 = XYZ[:,a1,:]
            pos2 = XYZ[:,a2,:]
            pos3 = XYZ[:,a3,:]
        
            #update in position affects all features, so recalculate
            v21 = (pos1 - pos2) 
            v21 /= v21.pow(2).sum(dim=-1, keepdims=True).sqrt()
            v32 = (pos2 - pos3)
            v32 /= v32.pow(2).sum(dim=-1, keepdims=True).sqrt()
        
            vp = torch.cross(v32, v21)
            cs = (v21 * v32).sum(dim=-1, keepdims=True)
            
            sn = (1.0 - cs * cs).sqrt()
            vp = vp / sn
            vu = torch.cross(vp, v21)

            XYZ[:,a0,:] = pos1 + \
                bonds[:,i:i+1]*(vu *sn_ang[:,i:i+1] * cs_tor[:,i:i+1] + \
                                vp *sn_ang[:,i:i+1] * sn_tor[:,i:i+1] - \
                                v21*cs_ang[:,i:i+1])
          
            i = i + 1
                  
        return XYZ


if __name__=='__main__':

     
    data_path = '../../data/14sbff_alpha/sorted/train/14sbff_alpha_train.dcd'
    prmtop_path = './14sbff_alpha.prmtop'
    
    device = 'cuda'
    u = mda.Universe(prmtop_path,data_path)
    selected = u.select_atoms("protein")
    R = BAT(selected)
       
    xyz = []
    for i in u._trajectory:
        xyz.append(i.positions)
    xyz = np.array(xyz)
    xyz = torch.as_tensor(xyz, device=device, dtype=torch.float32)
    xyz.requires_grad=True


    B = BATT(R, device=device)    
    root_3_xyz = B.GetRootXYZ(xyz)
    n1, n2, va, vb, bonds = B.Coords2MainVecs(xyz)
    xyz_new = B.BatV2Coords(n1, n2, va, vb, bonds, root_3_xyz)

    rmsd = (xyz - xyz_new).pow(2).mean().sqrt()
    print('BAT to Coords reconstruction Err --' , rmsd.item())
   
    xyz_new = xyz_new.detach().cpu().numpy()
    with mda.Writer('out.dcd', n_atoms=u.atoms.n_atoms) as w:
        i = 0
        for ts in u.trajectory:
            ts.positions = xyz_new[i,:,:]
            i = i + 1
            w.write(u.atoms)