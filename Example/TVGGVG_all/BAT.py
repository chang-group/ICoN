import numpy as np
import torch
import math

torch.pi = math.pi
import pytraj as pt


def Coords2MainVecs(xyz, tor_indxs):
    """

    """

    # get all positions for torsion calc
    p1 = xyz[:, tor_indxs[:, 0]]
    p2 = xyz[:, tor_indxs[:, 1]]
    p3 = xyz[:, tor_indxs[:, 2]]
    p4 = xyz[:, tor_indxs[:, 3]]

    # get all bond vecs for all bat
    va = p2 - p1
    vb = p3 - p2  # middle bond vector
    vc = p4 - p3

    n1 = torch.cross(-va, vb)  # n1 is normal vector to -va, vb
    n2 = torch.cross(-vb, vc)  # n2 is normal vector to -vb, vc

    return n1, n2, va, vb


def BondAngleTorsion(n1, n2, va, vb, primary_torsion_indices=None):

    # get bond length
    bonds = va.pow(2).sum(dim=-1).sqrt()

    # get bond angles
    x1 = (-va * vb).sum(dim=-1)
    y1 = torch.cross(va, vb).pow(2).sum(axis=-1).sqrt()
    angles = torch.atan2(y1, x1)

    # get torsions
    xp = torch.cross(n1, n2)
    x2 = (n1 * n2).sum(dim=-1)
    y2 = (xp * vb).sum(dim=-1) / (vb * vb).sum(dim=-1).sqrt()
    torsions = torch.atan2(y2, x2)

    # torsions = _ApplyTorShift(torsions)

    if primary_torsion_indices is not None:
        shift = torsions[:, primary_torsion_indices]
        unique_primary_torsion_indices = list(set(primary_torsion_indices))
        shift[:, unique_primary_torsion_indices] = 0.0
        torsions -= shift
        torsions = ((torsions + torch.pi) % (2 * torch.pi)) - torch.pi

    return bonds, angles, torsions


def GetRoot(xyz, root_XYZ_inds, device='cuda'):
    """

    """

    n_frames = xyz.shape[0]

    p0 = xyz[:, root_XYZ_inds[0], :]
    p1 = xyz[:, root_XYZ_inds[1], :]
    p2 = xyz[:, root_XYZ_inds[2], :]
    p3 = xyz[:, root_XYZ_inds[3], :]
    p4 = xyz[:, root_XYZ_inds[4], :]

    v01 = p1 - p0
    v21 = p1 - p2
    v31 = p1 - p3
    v41 = p1 - p4

    # Internal coordinates
    r01 = v01.pow(2).sum(dim=-1).sqrt()  # Distance between first two root atoms
    r12 = v21.pow(2).sum(dim=-1).sqrt()  # Distance between second two root atoms
    r13 = v31.pow(2).sum(dim=-1).sqrt()  # Distance between third two root atoms
    r14 = v41.pow(2).sum(dim=-1).sqrt()  #Distance between fourth two root atoms

    # Angle between root atoms
    a012 = torch.acos((v01 * v21).sum(dim=-1) / (r01 * r12))
    a013 = torch.acos((v01 * v31).sum(dim=-1) / (r01 * r13))
    a014 = torch.acos((v01 * v41).sum(dim=-1) / (r01 * r14))

    # External coordinates
    e = v01 / r01.unsqueeze(dim=-1)
    phi = torch.atan2(e[:, 1], e[:, 0])  # Polar angle
    theta = torch.acos(e[:, 2])  # Azimuthal angle

    # Rotation to the z axis
    cp = torch.cos(phi);
    sp = torch.sin(phi)
    ct = torch.cos(theta);
    st = torch.sin(theta)

    Rz = torch.zeros(n_frames, 3, 3, device=device)
    Rz[:, 0, 0] = cp * ct;
    Rz[:, 0, 1] = ct * sp;
    Rz[:, 0, 2] = -st
    Rz[:, 1, 0] = -sp;
    Rz[:, 1, 1] = cp;
    Rz[:, 1, 2] = 0.0
    Rz[:, 2, 0] = cp * st;
    Rz[:, 2, 1] = sp * st;
    Rz[:, 2, 2] = ct

    pos2 = torch.zeros(n_frames, 3, device=device)
    for i in range(n_frames):
        pos2[i, :] = torch.matmul(Rz[i, :, :], -v21[i, :])

    # Angle about the rotation axis
    omega = torch.atan2(pos2[:, 1], pos2[:, 0])

    #         root_based = torch.cat([p0, phi.unsqueeze(dim=-1),
    #                                     theta.unsqueeze(dim=-1),
    #                                     omega.unsqueeze(dim=-1),
    #                                     r01.unsqueeze(dim=-1),
    #                                     r12.unsqueeze(dim=-1),
    #                                     a012.unsqueeze(dim=-1)], dim=-1)
    #                                     #a013.unsqueeze(dim=-1)

    root_based = torch.cat([p0, r01.unsqueeze(dim=-1),
                            r12.unsqueeze(dim=-1),
                            a012.unsqueeze(dim=-1),
                            p3, p4], dim=-1)
    # root_based = torch.cat([p0, r01.unsqueeze(dim=-1),
    #                         r12.unsqueeze(dim=-1),
    #                         a012.unsqueeze(dim=-1),
    #                         p3], dim=-1)

    return root_based


def Bat2Coords(bat, root_XYZ_inds, tor_indxs, primary_torsion_indices, root_atom_xyz, device='cuda'):
    print("bat2coords called")
    root_buffer = (len(root_XYZ_inds) - 1) * 3
    n_frames = bat.shape[0]
    n_tors = int((bat.shape[1] - 12) / 3)

    # Extract bat components
    bonds = bat[:, root_buffer:n_tors + root_buffer]
    angles = bat[:, n_tors + root_buffer:2 * n_tors + root_buffer]
    torsions = bat[:, 2 * n_tors + root_buffer:]

    sn_ang = torch.sin(angles)
    cs_ang = torch.cos(angles)
    sn_tor = torch.sin(torsions)
    cs_tor = torch.cos(torsions)



    XYZ = torch.zeros(n_frames, n_tors + len(root_XYZ_inds)-2, 3, device=device)
    XYZ[:, root_XYZ_inds[0], :] = root_atom_xyz[0]
    XYZ[:, root_XYZ_inds[1], :] = root_atom_xyz[1]
    XYZ[:, root_XYZ_inds[2], :] = root_atom_xyz[2]
    XYZ[:, root_XYZ_inds[3], :] = root_atom_xyz[3]
    XYZ[:, root_XYZ_inds[4], :] = root_atom_xyz[4]

    eps = 1e-8
    for i, (a0, a1, a2, a3) in enumerate(tor_indxs):
        pos1 = XYZ[:, a1, :]
        pos2 = XYZ[:, a2, :]
        pos3 = XYZ[:, a3, :]

        v21 = pos1 - pos2
        v32 = pos2 - pos3

        v21_norm = torch.norm(v21, dim=-1, keepdim=True)
        v32_norm = torch.norm(v32, dim=-1, keepdim=True)
        v21 = v21 / (v21_norm + eps)
        v32 = v32 / (v32_norm + eps)

        
        #vp = torch.cross(v32, v21)
        vp = torch.linalg.cross(v32, v21)
        
        cs = torch.clamp((v21 * v32).sum(dim=-1, keepdim=True), -1.0 + eps, 1.0 - eps)
        sn = torch.sqrt(1.0 - cs * cs + eps)
        vp = vp / (sn + eps)
        
        #vu = torch.cross(vp, v21)
        vu = torch.linalg.cross(vp, v21)

        new_pos = pos1 + bonds[:, i:i + 1] * (
                vu * sn_ang[:, i:i + 1] * cs_tor[:, i:i + 1] +
                vp * sn_ang[:, i:i + 1] * sn_tor[:, i:i + 1] -
                v21 * cs_ang[:, i:i + 1]
        )

        XYZ[:, a0, :] = new_pos

    return XYZ
