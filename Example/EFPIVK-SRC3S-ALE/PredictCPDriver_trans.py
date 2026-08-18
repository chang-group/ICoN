import torch
import numpy as np
from tqdm import tqdm
from GetTorsionList import *
from BAT import *
import pytraj as pt
#from FCNetDrvier import fc_encode, fc_decode
from Transform_test import TransformerModel


def normalize_sin_cos(sin_t, cos_t):
    norm = torch.sqrt(torch.pow(sin_t, 2) + torch.pow(cos_t, 2))
    return torch.div(sin_t, norm), torch.div(cos_t, norm)
def reconstruct_full_bat(primary_bat, angle_differences, all_bonds, all_angles, primary_indices, non_primary_bonds,
                         non_primary_angles, non_primary_indices, torsion_XYZ_inds):
    batch_size, n_primary = primary_bat.shape[0], primary_bat.shape[1] // 4
    n_total = len(all_bonds)

    # Initialize tensors for each component
    bonds = torch.zeros((batch_size, n_total), device=primary_bat.device)
    angles = torch.zeros((batch_size, n_total), device=primary_bat.device)
    sin_torsions = torch.zeros((batch_size, n_total), device=primary_bat.device)
    cos_torsions = torch.zeros((batch_size, n_total), device=primary_bat.device)

    bonds[:, primary_indices] = primary_bat[:, :n_primary]
    angles[:, primary_indices] = primary_bat[:, n_primary:2*n_primary]

    bonds[:, non_primary_indices] = non_primary_bonds.unsqueeze(0).expand(batch_size, -1)
    angles[:, non_primary_indices] = non_primary_angles.unsqueeze(0).expand(batch_size, -1)

    # Set primary torsions
    sin_torsions[:, primary_indices] = primary_bat[:, 2*n_primary:3*n_primary]
    cos_torsions[:, primary_indices] = primary_bat[:, 3*n_primary:]

    # Reconstruct other torsions
    for key, diffs in angle_differences.items():
        primary_idx = primary_indices[list(angle_differences.keys()).index(key)]
        primary_torsion = torch.atan2(sin_torsions[:, primary_idx], cos_torsions[:, primary_idx])
        for idx, diff in diffs.items():
            reconstructed_torsion = (primary_torsion + diff + np.pi) % (2 * np.pi) - np.pi
            sin_torsions[:, idx] = torch.sin(reconstructed_torsion)
            cos_torsions[:, idx] = torch.cos(reconstructed_torsion)

    full_bat = torch.cat([bonds, angles, sin_torsions, cos_torsions], dim=1)

    return full_bat

def eval_model(dat, modelE, n_tors, angle_differences, all_bonds, all_angles, primary_indices,
               non_primary_bonds, non_primary_angles, non_primary_indices, device='cuda'):
    bond, angle, torsion = dat
    bond = dat[0]
    angle = dat[1]
    torsion = dat[2]
    batch_size = torsion.shape[0]
    n_bonds = bond.shape[1]
    n_angles = angle.shape[1]
    n_tors = torsion.shape[1]
    feats = torch.cat([bond, angle, torch.sin(torsion), torch.cos(torsion)], dim=-1)

    feats = feats.to(device)
    feats = feats.unsqueeze(1)

    modelE.eval()
    #modelD.eval()

    #z = modelE(feats)  # Latent Representation
    #out = modelD(z)  # Only has BAT for primary atoms
    out = modelE(feats)
    
    bonds = out[:, :n_bonds]
    angles = out[:, n_bonds:n_bonds + n_angles]
    sin_t = out[:, n_bonds + n_angles: n_bonds + n_angles + n_tors]
    cos_t = out[:, n_bonds + n_angles + n_tors:]
    sin_t_norm, cos_t_norm = normalize_sin_cos(sin_t, cos_t)
    out = torch.cat([bonds, angles, sin_t_norm, cos_t_norm], dim=-1)

    full_bat = reconstruct_full_bat(out, angle_differences, all_bonds, all_angles, primary_indices,
                                    non_primary_bonds, non_primary_angles, non_primary_indices, torsion_XYZ_inds)
    # Reconstruct the BAT for all atoms from the BAT for primary atoms

    return out, full_bat

if __name__ == '__main__':
    par_id = 20000
    device = 'cuda'
    dev_id = 0
    batch_size = 64  # Make this as large as your PC can handle

    out_path = "./out/"
    out_path_params = out_path + "Driver_XYZ_params/"
    params_file_name = 'Driver_XYZ_Params'
    top_path = "../traj/cyclicpeptide/top.prmtop"
    path = "../traj/cyclicpeptide/Cyclic_Peptide_MinimizedSorted.dcd"
    pdb_name = "Driver_Cyclic_Peptide_Prediction"


    traj = pt.load(path, top_path)
    n_frames = traj.n_frames

    print(f"Total frames to process: {n_frames}")

    torsion_XYZ_inds, n_atoms, root = getTorsionList(path, top_path)
    torsion_XYZ_inds = np.array(torsion_XYZ_inds)

    # Load pre-trained model and saved torsion information
    checkpoint = torch.load(out_path_params + str(par_id) + params_file_name)
    primary_torsions = checkpoint['primary_torsions']
    angle_differences = checkpoint['angle_differences']
    primary_indices = checkpoint['primary_indices']
    all_bonds = checkpoint['all_bonds'].to(device)
    all_angles = checkpoint['all_angles'].to(device)
    non_primary_bonds = torch.tensor(checkpoint['non_primary_bonds'], device=device)  # From only first frame
    non_primary_angles = torch.tensor(checkpoint['non_primary_angles'], device=device)  # From only first frame
    non_primary_indices = checkpoint['non_primary_indices']

    n_torsions = len(primary_torsions)
    n_feats = 4 * n_torsions  # bond, angle, sin(torsion), cos(torsion)
    

    torch.cuda.set_device(dev_id)
    #modelE = fc_encode(n_feats).to(device)
    #modelD = fc_decode(n_feats).to(device)
   
    #feature_dim = bond.shape[1] + angle.shape[1] + 2 * torsion.shape[1]  # Adjust as needed
    n_heads = 4  # Number of attention heads
    num_layers = 4  # Number of transformer layers
    hidden_dim = 128  # Size of hidden layer
    dropout = 0.1  # Dropout rate

    # Initialize the transformer model
    modelE = TransformerModel(n_feats, n_heads, num_layers, hidden_dim, dropout).to(device)

    modelE.load_state_dict(checkpoint['modelEncode'])
    #modelD.load_state_dict(checkpoint['modelDecode'])

    root_XYZ_inds = [root[0].index, root[1].index, root[2].index, root[3].index, root[4].index]
    prior_atoms = [sorted([a1, a2]) for (a0, a1, a2, a3) in torsion_XYZ_inds]
    primary_torsion_indices = [prior_atoms.index(prior_atoms[n]) for n in range(len(prior_atoms))]

    all_z = []
    all_coords = []

    for i in tqdm(range(0, n_frames, batch_size), desc="Processing batches"):
        batch_frames = slice(i, min(i + batch_size, n_frames))
        xyz = torch.as_tensor(traj.xyz[batch_frames], device=device, dtype=torch.float32)

        # Calculate bonds and angles for the current batch
        n1, n2, va, vb = Coords2MainVecs(xyz, torsion_XYZ_inds)
        current_bat = BondAngleTorsion(n1, n2, va, vb)
        current_bonds = current_bat[0][1]
        current_angles = current_bat[1][1]

        # Now get features only for primary torsions
        n1, n2, va, vb = Coords2MainVecs(xyz, primary_torsions)
        dat = BondAngleTorsion(n1, n2, va, vb)

        z, full_bat = eval_model(dat, modelE, n_torsions, angle_differences, current_bonds, current_angles,
                                 primary_indices, non_primary_bonds, non_primary_angles, non_primary_indices,
                                 device=device)

        all_z.append(z.detach().cpu().numpy())

        n_total = len(all_bonds)
        bondP = full_bat[:, :n_total]
        angleP = full_bat[:, n_total:2*n_total]
        sin_torsionP = full_bat[:, 2*n_total:3*n_total]
        cos_torsionP = full_bat[:, 3*n_total:]
        torsionP = torch.atan2(sin_torsionP, cos_torsionP)

        root_based = GetRoot(xyz, root_XYZ_inds)
        root_atoms_XYZ = [xyz[:, i, :] for i in root_XYZ_inds]

        bat = torch.cat([root_based, bondP, angleP, torsionP], dim=-1)

        xyzbat = Bat2Coords(bat, root_XYZ_inds, torsion_XYZ_inds, primary_torsion_indices, root_atoms_XYZ)

        coords = xyzbat.detach().cpu().numpy()
        all_coords.append(coords)

    all_z = np.concatenate(all_z, axis=0)
    np.save(out_path + 'visual/' + pdb_name + 'latent' + '.npy', all_z)

    all_coords = np.concatenate(all_coords, axis=0)

    print('Min max X Coordinate:', np.min(all_coords[:, :, 0]), np.max(all_coords[:, :, 0]))
    print('Min max Y Coordinate:', np.min(all_coords[:, :, 1]), np.max(all_coords[:, :, 1]))
    print('Min max Z Coordinate:', np.min(all_coords[:, :, 2]), np.max(all_coords[:, :, 2]))

    traj1 = pt.Trajectory(top=top_path)
    traj1.xyz = all_coords
    pt.write_traj(out_path + pdb_name + '.dcd', traj1, overwrite=True)

    print("Prediction complete. Output saved to:", out_path + pdb_name + '.dcd')
