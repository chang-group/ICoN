import numpy as np
import torch
from tqdm import tqdm
import pytraj as pt
from PredictCPDriver import eval_model, reconstruct_full_bat, fc_decode
from BAT import *
from GetTorsionList import *


def spherical_expansion(modelD, latent_coords, num_points=11544, radius=8, device='cuda'):
    # Find the center of the latent space
    center = np.mean(latent_coords, axis=0)

    # Generate random points in a sphere
    phi = np.random.uniform(0, 2 * np.pi, num_points)
    costheta = np.random.uniform(-1, 1, num_points)
    u = np.random.uniform(0, 1, num_points)

    theta = np.arccos(costheta)
    r = radius * np.cbrt(u)

    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    # Combine the coordinates and shift by the center
    sphere_points = np.column_stack((x, y, z)) + center

    # Convert to torch tensor
    sphere_points_tensor = torch.tensor(sphere_points, dtype=torch.float32, device=device)

    # Decode the points
    modelD.eval()
    with torch.no_grad():
        decoded_points = modelD(sphere_points_tensor)

    return sphere_points, decoded_points


# Main execution
if __name__ == '__main__':
    
    par_id = 26000
    device = 'cuda'
    dev_id = 0
    batch_size = 64  # Make this as large as your PC can handle
    
    out_path = "./out/"
    out_path_params = out_path + "Driver_XYZ_params/"
    params_file_name = 'Driver_XYZ_Params'
    top_path = "../traj/cyclicpeptide/top.prmtop"
    path = "../traj/cyclicpeptide/Cyclic_Peptide_MinimizedSorted.dcd"
    pdb_name = "Driver_Cyclic_Peptide_Prediction"


    # Load the latent coordinates
    latent_coords = np.load(out_path + 'visual/Driver_Cyclic_Peptide_Predictionlatent.npy')

    # Load the model and other necessary data
    checkpoint = torch.load(out_path + "Driver_XYZ_params/" + str(par_id) + params_file_name)
    modelD = fc_decode(4 * len(checkpoint['primary_torsions'])).to(device)
    modelD.load_state_dict(checkpoint['modelDecode'])

    # Other necessary data
    angle_differences = checkpoint['angle_differences']
    all_bonds = checkpoint['all_bonds'].to(device)
    all_angles = checkpoint['all_angles'].to(device)
    primary_indices = checkpoint['primary_indices']
    non_primary_bonds = torch.tensor(checkpoint['non_primary_bonds'], device=device)
    non_primary_angles = torch.tensor(checkpoint['non_primary_angles'], device=device)
    non_primary_indices = checkpoint['non_primary_indices']

    # Perform spherical expansion
    sphere_points, decoded_points = spherical_expansion(modelD, latent_coords)

    # Save the new latent points
    np.save(out_path + 'visual/spherical_expansion_latent.npy', sphere_points)

    torsion_XYZ_inds, n_atoms, root = getTorsionList(path, top_path)
    # Reconstruct full BAT and convert to coordinates
    full_bat = reconstruct_full_bat(decoded_points,
                                    angle_differences, all_bonds, all_angles,
                                    primary_indices, non_primary_bonds,
                                    non_primary_angles, non_primary_indices,
                                    torsion_XYZ_inds)

    root_XYZ_inds = [root[0].index, root[1].index, root[2].index, root[3].index, root[4].index]
    prior_atoms = [sorted([a1, a2]) for (a0, a1, a2, a3) in torsion_XYZ_inds]
    primary_torsion_indices = [prior_atoms.index(prior_atoms[n]) for n in range(len(prior_atoms))]
    traj = pt.load(path, top_path)
    xyz = torch.as_tensor(traj.xyz, device=device, dtype=torch.float32)
    root_based = GetRoot(xyz, root_XYZ_inds)
    root_atoms_XYZ = [xyz[:, i, :] for i in root_XYZ_inds]
    n_total = len(all_bonds)
    #n_total = 10000
    
    bondP = full_bat[:, :n_total]
    angleP = full_bat[:, n_total:2 * n_total]
    sin_torsionP = full_bat[:, 2 * n_total:3 * n_total]
    cos_torsionP = full_bat[:, 3 * n_total:]
    torsionP = torch.atan2(sin_torsionP, cos_torsionP)

    print(root_atoms_XYZ[0][0])
    print(root_XYZ_inds)

    #rb = root_based[0:10000,:]
    #bat = torch.cat([rb, bondP, angleP, torsionP], dim=-1)
    bat = torch.cat([root_based, bondP, angleP, torsionP], dim=-1)
    coords = Bat2Coords(bat, root_XYZ_inds, torsion_XYZ_inds, primary_torsion_indices, root_atoms_XYZ)


    # Save as trajectory
    traj = pt.Trajectory(top=top_path)
    traj.xyz = coords.cpu().numpy()
    pt.write_traj(out_path + 'spherical_expansion.dcd', traj, overwrite=True)

    print("Spherical expansion complete. Output saved to:", out_path + 'spherical_expansion.dcd')
    print("New latent points saved to:", out_path + 'visual/spherical_expansion_latent.npy')
