import os
import math

os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'max_split_size_mb:32'
import gc

import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import tqdm

from BAT import *
from Atom import *
import pytraj as pt
from GetTorsionList import *
from calculate_Energy import *
from Dihedral_Transformer import *
#from Driver import reconstruct_full_bat
#from pcgrad import PCGrad


class WarmupCosineScheduler:
    def __init__(self, optimizer, warmup_iterations, total_iterations, min_lr=1e-6):
        self.optimizer = optimizer
        self.warmup_iterations = warmup_iterations
        self.total_iterations = total_iterations
        self.min_lr = min_lr
        self.base_lr = optimizer.param_groups[0]['lr']
        self.current_iteration = 0

    def step(self):
        self.current_iteration += 1

        if self.current_iteration < self.warmup_iterations:
            # Linear warmup
            lr = self.base_lr * (self.current_iteration / self.warmup_iterations)
        else:
            # Cosine decay
            progress = (self.current_iteration - self.warmup_iterations) / (
                        self.total_iterations - self.warmup_iterations)
            lr = self.min_lr + 0.5 * (self.base_lr - self.min_lr) * (1 + np.cos(progress * np.pi))

        for param_group in self.optimizer.param_groups:
            param_group['lr'] = lr

        return lr


def setup_XYZ_dat(all_dat, IT, primary_indices, all_bonds, all_angles):
    # save non-primary bat for each frame during training
    non_primary_indices = [i for i in range(len(torsion_XYZ_inds)) if i not in primary_indices]
    non_primary_bonds = torch.stack([all_dat[0][IT, i] for i in non_primary_indices], dim=-1)
    non_primary_angles = torch.stack([all_dat[1][IT, i] for i in non_primary_indices], dim=-1)
    bonds = all_dat[0]
    angles = all_dat[1]
    torsion_angles = all_dat[2]
    primary_torsions = []
    angle_differences_all = {}
    primary_indices = []
    all_bonds = bonds[IT, :]
    all_angles = angles[IT, :]

    for i, (a, b, c, d) in enumerate(torsion_XYZ_inds):

        key = (b, c, d)
        if key not in angle_differences_all:
            primary_torsions.append([a, b, c, d])
            angle_differences_all[key] = {}
            primary_indices.append(i)
        else:
            primary_idx = primary_indices[list(angle_differences_all.keys()).index(key)]
            diff_all = (torsion_angles[IT, i] - torsion_angles[IT, primary_idx] + np.pi) % (2 * np.pi) - np.pi
            angle_differences_all[key][i] = diff_all

    return primary_torsions, angle_differences_all, primary_indices, all_bonds, all_angles, non_primary_bonds, non_primary_angles, non_primary_indices


def reconstruct_XYZ(xyz, out, xyz_data, name):
    # Very similar to reconstruct_full_bat, but bonds,angles,and saved indices are frame specific

    primary_torsions = xyz_data[0]  # All saved values are frame specific rather than the first frame
    angle_differences_all = xyz_data[1]
    primary_indices = xyz_data[2]
    all_bonds = xyz_data[3]
    all_angles = xyz_data[4]
    non_primary_bonds = xyz_data[5]
    non_primary_angles = xyz_data[6]
    non_primary_indices = xyz_data[7]

    batch_size, n_primary = out.shape[0], out.shape[1] // 4
    n_total = len(all_bonds[1])

    # Initialize tensors for each component
    bonds = torch.zeros((batch_size, n_total), device=out.device)
    angles = torch.zeros((batch_size, n_total), device=out.device)
    sin_torsions = torch.zeros((batch_size, n_total), device=out.device)
    cos_torsions = torch.zeros((batch_size, n_total), device=out.device)

    bonds[:, primary_indices] = out[:, :n_primary]
    angles[:, primary_indices] = out[:, n_primary:2 * n_primary]

    bonds[:, non_primary_indices] = non_primary_bonds
    angles[:, non_primary_indices] = non_primary_angles

    # Set primary torsions
    sin_torsions[:, primary_indices] = out[:, 2 * n_primary:3 * n_primary]
    cos_torsions[:, primary_indices] = out[:, 3 * n_primary:]

    # Reconstruct other torsions
    for key, diffs in angle_differences_all.items():
        primary_idx = primary_indices[list(angle_differences_all.keys()).index(key)]
        primary_torsion = torch.atan2(sin_torsions[:, primary_idx].clone(), cos_torsions[:, primary_idx].clone())
        for idx, diff in diffs.items():
            reconstructed_torsion = (primary_torsion + diff + np.pi) % (2 * np.pi) - np.pi
            sin_torsions[:, idx] = torch.sin(reconstructed_torsion)
            cos_torsions[:, idx] = torch.cos(reconstructed_torsion)

    full_bat = torch.cat([bonds, angles, sin_torsions, cos_torsions], dim=1)

    root_XYZ_inds = [root[0].index, root[1].index, root[2].index, root[3].index, root[4].index]
    n_total = len(all_bonds[1])
    bondP = full_bat[:, :n_total]
    angleP = full_bat[:, n_total:2 * n_total]
    sin_torsionP = full_bat[:, 2 * n_total:3 * n_total]
    cos_torsionP = full_bat[:, 3 * n_total:]
    torsionP = torch.atan2(sin_torsionP, cos_torsionP)

    root_based = GetRoot(xyz, root_XYZ_inds)
    root_atoms_XYZ = [xyz[:, i, :] for i in root_XYZ_inds]

    bat = torch.cat([root_based, bondP, angleP, torsionP], dim=-1)

    xyz_bat = Bat2Coords(bat, root_XYZ_inds, torsion_XYZ_inds, primary_indices, root_atoms_XYZ)
    xyz_bat = xyz_bat.clone()

    coords = xyz_bat.detach().cpu().numpy()
    traj1 = pt.Trajectory(top=top_path)
    traj1.xyz = coords
    pdb_name = name
    pt.write_traj(out_path + pdb_name + '.dcd', traj1, overwrite=True)
    print("Prediction complete. Output saved to:", out_path + pdb_name + '.dcd')
    path_prediction = out_path + pdb_name + '.dcd'
    return xyz_bat, path_prediction

def get_inter_pts(z, batch_size):
    z_idx = torch.arange(batch_size, device=z.device)
    pairs = torch.combinations(z_idx, r=2)

    ras_endpts = pairs[torch.randperm(len(pairs))[:batch_size]]

    start_points = z[ras_endpts[:, 0]]
    end_points = z[ras_endpts[:, 1]]

    return start_points, end_points

def slerp_batch(p0, p1, t):
    p0_norm = torch.nn.functional.normalize(p0, dim=-1)
    p1_norm = torch.nn.functional.normalize(p1, dim=-1)

    dot = torch.sum(p0_norm * p1_norm, dim=-1, keepdim=True)
    dot = torch.clamp(dot, -1.0 + 1e-7, 1.0 - 1e-7)
    omega = torch.acos(dot)

    small_angle_mask = omega < 1e-3
    t = t.unsqueeze(1)

    lerp = p0 * (1.0 - t) + p1 * t

    sin_omega = torch.sin(omega)
    sin_omega = torch.where(sin_omega < 1e-7, torch.ones_like(sin_omega) * 1e-7, sin_omega)

    # Calculate SLERP
    scale0 = torch.sin((1.0 - t) * omega)
    scale1 = torch.sin(t * omega)
    slerp = (scale0 * p0 + scale1 * p1) / sin_omega

    result = torch.where(small_angle_mask, lerp, slerp)

    return result


def sample_inter_pts(z, batch_size):
    t_range = 1000

    start_points, end_points = get_inter_pts(z, batch_size)
    ras_t = torch.rand(batch_size, device=z.device)
    t = ras_t * (t_range - 1) / t_range

    interpolated_points = slerp_batch(start_points, end_points, t)

    return interpolated_points


def reconstrcut_xyz_inter(xyz, out, xyz_data, f0_bonds, f0_angles, f0_tors, name):

    primary_torsions = xyz_data[0]  # All saved values are frame specific rather than the first frame
    angle_differences_all = xyz_data[1]
    primary_indices = xyz_data[2]
    all_bonds = xyz_data[3]
    all_angles = xyz_data[4]
    non_primary_bonds = xyz_data[5]
    non_primary_angles = xyz_data[6]
    non_primary_indices = xyz_data[7]

    batch_size, n_primary = out.shape[0], out.shape[1] // 4
    n_total = len(all_bonds[1])

    # Initialize tensors for each component
    bonds = torch.zeros((batch_size, n_total), device=out.device)
    angles = torch.zeros((batch_size, n_total), device=out.device)
    sin_torsions = torch.zeros((batch_size, n_total), device=out.device)
    cos_torsions = torch.zeros((batch_size, n_total), device=out.device)

    bonds[:, primary_indices] = out[:, :n_primary]
    angles[:, primary_indices] = out[:, n_primary:2 * n_primary]

    bonds[:, non_primary_indices] = non_primary_bonds
    angles[:, non_primary_indices] = non_primary_angles

    # Set primary torsions
    sin_torsions[:, primary_indices] = out[:, 2 * n_primary:3 * n_primary]
    cos_torsions[:, primary_indices] = out[:, 3 * n_primary:]

    # Reconstruct other torsions
    for key, diffs in angle_differences_all.items():
        primary_idx = primary_indices[list(angle_differences_all.keys()).index(key)]
        primary_torsion = torch.atan2(sin_torsions[:, primary_idx].clone(), cos_torsions[:, primary_idx].clone())
        for idx, diff in diffs.items():
            reconstructed_torsion = (primary_torsion + diff + np.pi) % (2 * np.pi) - np.pi
            sin_torsions[:, idx] = torch.sin(reconstructed_torsion)
            cos_torsions[:, idx] = torch.cos(reconstructed_torsion)

    full_bat = torch.cat([bonds, angles, sin_torsions, cos_torsions], dim=1)

    root_XYZ_inds = [root[0].index, root[1].index, root[2].index, root[3].index, root[4].index]
    n_total = len(all_bonds[1])
    bondP = full_bat[:, :n_total]
    angleP = full_bat[:, n_total:2 * n_total]
    sin_torsionP = full_bat[:, 2 * n_total:3 * n_total]
    cos_torsionP = full_bat[:, 3 * n_total:]
    torsionP = torch.atan2(sin_torsionP, cos_torsionP)

    root_based = GetRoot(xyz, root_XYZ_inds)
    root_atoms_XYZ = [xyz[:, i, :] for i in root_XYZ_inds]

    bat = torch.cat([root_based, bondP, angleP, torsionP], dim=-1)

    xyz_bat = Bat2Coords(bat, root_XYZ_inds, torsion_XYZ_inds, primary_indices, root_atoms_XYZ)
    xyz_bat = xyz_bat.clone()

    coords = xyz_bat.detach().cpu().numpy()
    traj1 = pt.Trajectory(top=top_path)
    traj1.xyz = coords
    pdb_name = name
    pt.write_traj(out_path + pdb_name + '.dcd', traj1, overwrite=True)
    print("Prediction complete. Output saved to:", out_path + pdb_name + '.dcd')
    path_prediction = out_path + pdb_name + '.dcd'
    return xyz_bat


def process_torsions(torsion_XYZ_inds, torsion_angles, bonds, angles):
    primary_torsions = []
    angle_differences = {}
    primary_indices = []
    all_bonds = bonds[0, :]  # Save all bonds from the first frame
    all_angles = angles[0, :]  # Save all angles from the first frame

    for i, (a, b, c, d) in enumerate(torsion_XYZ_inds):
        key = (b, c, d)
        if key not in angle_differences:
            primary_torsions.append([a, b, c, d])
            angle_differences[key] = {}
            primary_indices.append(i)
        else:
            primary_idx = primary_indices[list(angle_differences.keys()).index(key)]
            diff = (torsion_angles[0, i] - torsion_angles[0, primary_idx] + np.pi) % (2 * np.pi) - np.pi
            angle_differences[key][i] = diff

    return np.array(primary_torsions), angle_differences, primary_indices, all_bonds, all_angles
def calculate_loss_weights(iteration, max_iterations, rmsd):
    progress = iteration / max_iterations

    # used to calculate scaling penalty for energy terms based on RMSD
    inv_rmsd = 1.0 / (10.0 ** rmsd)
    exp_rmsd = torch.exp(rmsd)
    bonded_penalty = exp_rmsd * inv_rmsd
    nb_penalty = inv_rmsd

    bonded_energy_weight = torch.tensor(0.0, device=device)
    non_bonded_energy_weight = torch.tensor(0.0, device=device)

    # Gradually introduce bonded energy (20-70%), scale down if RMSD becomes too high
    if progress > 0.3:
        weight = torch.tensor(min(1.0, (progress - 0.3) / 0.5), device=device)
        bonded_energy_weight = weight * bonded_penalty
        non_bonded_energy_weight = weight * nb_penalty

    return bonded_energy_weight, non_bonded_energy_weight

def run_model(dat, xyz_data, xyz, model, indx, iterations, iteration, fixed_tors_ind, flexibility, curr_batch_size, f0_bonds,
              f0_angles, f0_tors, train=True, device='cuda'):

    xyz = xyz[indx, :, :]
    bonds = dat[0][indx, :]
    angles = dat[1][indx, :]
    torsion = dat[2][indx, :]
    batch_size = torsion.shape[0]
    n_tors = torsion.shape[1]

    input_dihedrals = torch.cat([torch.sin(torsion), torch.cos(torsion)], dim=-1)
    flexible_mask = torch.ones(n_tors, dtype=torch.bool, device=device)
    flexible_mask[fixed_tors_ind] = False
    fixed_torsions = torsion[:, fixed_tors_ind]
    flexible_torsions = torsion[:, flexible_mask]

    transformer_input = torch.stack([torch.sin(flexible_torsions), torch.cos(flexible_torsions)], dim=-1)

    if train:
        model.train()
    else:
        model.eval()
    criterion = nn.MSELoss(reduction='mean')

    transformer_output, z = model(transformer_input)

    if (epoch <= 1):
        loss = criterion(transformer_output, transformer_input)
        print("Dihedral Loss:", loss)
        energy, rmsd = torch.tensor(0.0, device='cuda'), torch.tensor(0.0, device='cuda')

    else:
        flexible_sin_torsions = transformer_output[:, :, 0]
        flexible_cos_torsions = transformer_output[:, :, 1]

        sin_torsions = torch.zeros((curr_batch_size, n_tors), device=device)
        cos_torsions = torch.zeros((curr_batch_size, n_tors), device=device)

        sin_torsions[:, fixed_tors_ind] = torch.sin(fixed_torsions)
        cos_torsions[:, fixed_tors_ind] = torch.cos(fixed_torsions)
        sin_torsions[:, flexible_mask] = flexible_sin_torsions
        cos_torsions[:, flexible_mask] = flexible_cos_torsions
        out = torch.cat([bonds, angles, sin_torsions, cos_torsions], dim=-1)

        out_reshaped = out.reshape(curr_batch_size, -1)
        xyz_pred, path_prediction = reconstruct_XYZ(xyz, out_reshaped, xyz_data, 'Training_Driver_CP_Prediction')

        loss_xyz = criterion(xyz_pred, xyz)
        rmsd = loss_xyz.sqrt().detach()

        z_sample = sample_inter_pts(z, curr_batch_size)
        out_inter = model.decode(z_sample)
        flexible_sin_torsions_inter = out_inter[:, :, 0]
        flexible_cos_torsions_inter = out_inter[:, :, 1]
        sin_torsions_inter = torch.zeros((curr_batch_size, n_tors), device=device)
        cos_torsions_inter = torch.zeros((curr_batch_size, n_tors), device=device)

        sin_torsions_inter[:, fixed_tors_ind] = torch.sin(fixed_torsions)
        cos_torsions_inter[:, fixed_tors_ind] = torch.cos(fixed_torsions)
        sin_torsions_inter[:, flexible_mask] = flexible_sin_torsions_inter
        cos_torsions_inter[:, flexible_mask] = flexible_cos_torsions_inter

        out_inter = torch.cat([bonds, angles, sin_torsions_inter, cos_torsions_inter], dim=-1)
        out_inter_reshaped = out_inter.reshape(curr_batch_size, -1)
        xyz_inter = reconstrcut_xyz_inter(xyz, out_inter_reshaped, xyz_data, f0_bonds, f0_angles, f0_tors, "Sample")

        energy, bonded_energy, non_bonded_energy = calculate_energy(xyz_pred, top_path)
        energy_inter, bonded_energy_inter, non_bonded_energy_inter = calculate_energy(xyz_inter, top_path)
        bonded_weight, nb_weight = calculate_loss_weights(iteration, iterations, loss_xyz.detach())
        bonded_energy = bonded_energy * bonded_weight
        non_bonded_energy = non_bonded_energy * nb_weight
        bonded_energy_inter = bonded_energy_inter * bonded_weight
        non_bonded_energy_inter = non_bonded_energy_inter * nb_weight


        loss = loss_xyz + bonded_energy + non_bonded_energy + bonded_energy_inter + non_bonded_energy_inter
        #loss = loss_xyz + non_bonded_energy + non_bonded_energy_inter
        print("\nBonded Energy Weight:", bonded_weight.item(),
            "\nNon-bonded Energy Weight:", nb_weight.item(), )

        print("\n")
        print("Potential Energy:", energy.item())
        print("Interpolation Potential Energy:", energy_inter.item())
        print("XYZ Loss (RMSD):", rmsd.item())
        print("\n")

    if train:
        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        optimizer.step()
        current_lr = scheduler.step()
        # print("Learning Rate:", current_lr)

    return energy.item(), rmsd.item()


def get_bonds(atom, bond_inds, bond_k):
    num_bonds = 0
    bond_data = []
    for (atom1, atom2, bond_type_idx) in bond_inds:
        if atom1 == atom or atom2 == atom:
            k_value = bond_k[bond_type_idx - 1]
            if atom1 == atom:
                bond_data.append((atom1, atom2, k_value))
            else:
                bond_data.append((atom2, atom1, k_value))
            num_bonds += 1

    return num_bonds, bond_data


def is_correct_pi(bonds, atom1, atom2):
    sorted_bonds = sorted(bonds, key=lambda x: x[2])
    strongest_strength = sorted_bonds[-1][2]
    strongest_bonds = [bond for bond in sorted_bonds if bond[2] == strongest_strength]
    for bond in strongest_bonds:
        if ((bond[0] == atom1 and bond[1] == atom2) or
                (bond[0] == atom2 and bond[1] == atom1)):
            return True
    return False


def find_rigid_torsions(torsions):
    """
    Finds all dihedral angles that can't rotate due to a double or triple bond. Bonds rigid due to resonance structures
    (peptide c-n) are not counted. Rewrite if there is a better way to do this
    """
    torsion_tensor = torsions
    ttt = torch.tan(torsion_tensor)
    tt = torch.atan(ttt)
    tt_std = torch.std(tt, dim=0)
    parser = PRMTOPParser(top_path)
    prmtop_params = parser.get_prmtop_params()
    names = prmtop_params['atom_names']
    bond_indices = torch.tensor(prmtop_params['bond_indices'], dtype=torch.long, device='cuda')
    bond_k = prmtop_params['bond_k'].clone().detach().to(device='cuda', dtype=torch.float32)
    names_np = np.array(names)
    fixed_dihedral_inds = []

    for i, (dihedral, std) in enumerate(zip(primary_torsions, tt_std)):
        atom2_name = names_np[dihedral[1]]
        atom3_name = names_np[dihedral[2]]

        if 'O' in atom2_name or 'O' in atom3_name:
            continue
        elif 'N' in atom2_name or 'N' in atom3_name:
            N = dihedral[1] if atom2_name == 'N' else dihedral[2]
            num_bonds, bond_parms = get_bonds(N, bond_indices, bond_k)
            if num_bonds == 2:
                bonded_atom = None
                bonded_atom_name = None
                if 'N' in atom2_name:
                    bonded_atom_name = atom3_name
                    bonded_atom = dihedral[2]
                elif 'N' in atom3_name:
                    bonded_atom_name = atom2_name
                    bonded_atom = dihedral[1]
                if 'N' in bonded_atom_name:
                    num_bondsN, _ = get_bonds(bonded_atom, bond_indices, bond_k)
                    if num_bondsN == 1:
                        fixed_dihedral_inds.append(bonded_atom_name)
                elif is_correct_pi(bond_parms, dihedral[1], dihedral[2]):
                    fixed_dihedral_inds.append(i)
        elif 'C' in atom2_name or 'C' in atom3_name:
            C = None
            bonded_atom_name = None
            if 'C' in atom2_name:
                C = dihedral[1]
                bonded_atom_name = atom3_name
            elif 'C' in atom3_name:
                C = dihedral[2]
                bonded_atom_name = atom2_name
            num_bonds, bond_parms = get_bonds(C, bond_indices, bond_k)
            if num_bonds == 2:
                if 'N' in bonded_atom_name:
                    if is_correct_pi(bond_parms, dihedral[1], dihedral[2]):
                        fixed_dihedral_inds.append(i)
                elif tt_std[i] < 0.1:
                    fixed_dihedral_inds.append(i)
                else:
                    if is_correct_pi(bond_parms, dihedral[1], dihedral[2]):
                        fixed_dihedral_inds.append(i)
            elif num_bonds == 3:
                if is_correct_pi(bond_parms, dihedral[1], dihedral[2]):
                    fixed_dihedral_inds.append(i)
        elif 'P' in atom2_name or 'P' in atom3_name:
            P = dihedral[1] if atom2_name == 'N' else dihedral[2]
            num_bonds, bond_parms = get_bonds(P, bond_indices, bond_k)
            if num_bonds == 4:
                if is_correct_pi(bond_parms, dihedral[1], dihedral[2]):
                    fixed_dihedral_inds.append(i)
        elif 'S' in atom2_name or 'S' in atom3_name:
            S = None
            bonded_atom_name = None
            if 'S' in atom2_name:
                S = dihedral[1]
                bonded_atom_name = atom3_name
            elif 'S' in atom3_name:
                S = dihedral[2]
                bonded_atom_name = atom2_name
            num_bonds, bond_parms = get_bonds(S, bond_indices, bond_k)
            if num_bonds == 3 or num_bonds == 4:
                if 'N' in bonded_atom_name:
                    if is_correct_pi(bond_parms, dihedral[1], dihedral[2]):
                        fixed_dihedral_inds.append(i)

    mask = torch.ones(len(tt_std), dtype=torch.bool, device=tt_std.device)
    mask[fixed_dihedral_inds] = False
    training_dihedral_stds = tt_std[mask]

    # normalize from 0.1 - 1
    min_val = training_dihedral_stds.min()
    max_val = training_dihedral_stds.max()
    flexibility_score = 0.9 * (training_dihedral_stds - min_val) / (max_val - min_val) + 0.1

    print("\nFixed Dihedrals:")
    for i, dihedral in enumerate(primary_torsions):
        if i in fixed_dihedral_inds:
            print(f"Dihedral {i} [{names_np[dihedral[0]]}({dihedral[0]})-"
                  f"{names_np[dihedral[1]]}({dihedral[1]})-"
                  f"{names_np[dihedral[2]]}({dihedral[2]})-"
                  f"{names_np[dihedral[3]]}({dihedral[3]})]")

    return fixed_dihedral_inds, flexibility_score


def split_dataset_repeating(data, train_size=4, val_size=1):
    split_size = train_size + val_size  # Total size of each split
    num_splits = len(data) // split_size

    train_indices_list = []
    val_indices_list = []

    for i in range(num_splits):
        start_idx = i * split_size
        end_idx = start_idx + split_size

        train_indices = torch.arange(start_idx, start_idx + train_size)
        val_indices = torch.arange(start_idx + train_size, end_idx)

        # Append the indices to the respective lists
        train_indices_list.append(train_indices)
        val_indices_list.append(val_indices)

    train_indices = torch.cat(train_indices_list)
    val_indices = torch.cat(val_indices_list)

    return train_indices, val_indices


if __name__ == '__main__':
    torch.cuda.empty_cache()
    gc.collect()
    lrt = 0.0001
    max_epoch = 200
    start = 0
    device = 'cuda'
    dev_id = 0 

    batch_size = 64 
    spacing = "\t"

    
    out_path = "./out/"
    out_path_params = out_path + "Driver_XYZ_params/"
    params_file_name = 'Driver_XYZ_Params'
    # pdb_name = "C-pep"
    top_path = "../../traj/cyclicpeptide/ExoR2/ExoR2.prmtop"
    path = "../../traj/cyclicpeptide/ExoR2/Cyclic_Peptide_MinimizedSorted.dcd"
    pdb_name = "Driver_Cyclic_Peptide_Prediction"

    torch.cuda.set_device(dev_id)

    traj = pt.load(path, top_path)
    xyz = traj.xyz
    torch.cuda.empty_cache()
    xyz = torch.as_tensor(xyz, device='cuda', dtype=torch.float32)
    n_frames = xyz.shape[0]



    train_indices_list, val_indices_list = split_dataset_repeating(xyz)
    xyzT = xyz[train_indices_list]
    xyzV = xyz[val_indices_list]

    n_train = xyzT.shape[0]
    n_valid = xyzV.shape[0] 
    
    #n_train = int(0.8 * n_frames)
    #n_valid = n_frames - n_train

    #all_indices = torch.randperm(n_frames)
    #train_indices = all_indices[:n_train]
    #valid_indices = all_indices[n_train:]

    #xyzT = xyz[train_indices, :, :]
    #xyzV = xyz[valid_indices, :, :]

    print(f"Total frames: {n_frames}")
    print(f"Training frames: {n_train}")
    print(f"Validation frames: {n_valid}")

    torsion_XYZ_inds, n_atoms, root = getTorsionList(path, top_path, multifragment=False)
    torsion_XYZ_inds = np.array(torsion_XYZ_inds)

    # Get all torsion angles first
    n1T, n2T, vaT, vbT = Coords2MainVecs(xyzT, torsion_XYZ_inds)
    n1V, n2V, vaV, vbV = Coords2MainVecs(xyzV, torsion_XYZ_inds)

    all_datT = BondAngleTorsion(n1T, n2T, vaT, vbT)
    all_datV = BondAngleTorsion(n1V, n2V, vaV, vbV)

    print("all", all_datT[0][1])
    print(len(all_datT[0][1]))

    # Process torsions to get primary ones
    primary_torsions, angle_differences, primary_indices, all_bonds, all_angles = process_torsions(
        torsion_XYZ_inds, all_datT[2], all_datT[0], all_datT[1])

    print("Number of original torsions:", len(torsion_XYZ_inds))
    print("Number of primary torsions:", len(primary_torsions))
    np.set_printoptions(threshold=np.inf)
    print(torsion_XYZ_inds)
    print("\n\n\n")
    print(primary_torsions)

    # Now get features only for primary torsions

    n1T, n2T, vaT, vbT = Coords2MainVecs(xyzT, primary_torsions)
    n1V, n2V, vaV, vbV = Coords2MainVecs(xyzV, primary_torsions)

    datT = BondAngleTorsion(n1T, n2T, vaT, vbT)
    datV = BondAngleTorsion(n1V, n2V, vaV, vbV)

    print("primary", datT[0][1])
    print(len(datT[0][1]))

    bond = datT[0]
    angle = datT[1]
    torsion = datT[2]

    fixed_tors_ind, flexibility = find_rigid_torsions(torsion)

    n_tors = len(primary_torsions)
    f0_bonds = datT[0][0, :]
    f0_angles = datT[1][0, :]
    flexible_mask = torch.ones(n_tors, dtype=torch.bool, device=device)
    flexible_mask[fixed_tors_ind] = False
    f0_tors = datT[2][0, fixed_tors_ind]

    n_bonds = bond.shape[1]
    n_angles = angle.shape[1]
    n_tors = torsion.shape[1] - len(fixed_tors_ind)
    n_feats = 2 * n_tors  # sin and cos of primary torsions
    print("Number of features:", n_feats)

    #torch.cuda.set_device(dev_id)

    model = ICONTransformer(n_feats).to(device)
    model.apply(init_weights)

    optimizer = optim.Adam(model.parameters(), lr=lrt)

    iterations_per_epoch = math.ceil(n_train / batch_size)
    total_iterations = iterations_per_epoch * max_epoch
    print("Total Iterations", total_iterations)
    warmup_iterations = int(0.05 * total_iterations)
    scheduler = WarmupCosineScheduler(
        optimizer,
        warmup_iterations=warmup_iterations,
        total_iterations=total_iterations,
        min_lr=1e-6
    )

    # data for the first frame
    non_primary_indices = [i for i in range(len(torsion_XYZ_inds)) if i not in primary_indices]
    non_primary_bonds = [all_datT[0][0, i].item() for i in non_primary_indices]
    non_primary_angles = [all_datT[1][0, i].item() for i in non_primary_indices]

    iStart = start + 1
    iEnd = max_epoch + 1

    iteration = 0
    torch.autograd.set_detect_anomaly(True)
    with open(out_path + 'log_train_val_Cyclic_Peptide_Driver.txt', 'w') as fout:
        fout.write("Epoch\tTrain_RMSD\tVal_RMSD\tTrain_Energy\tVal_Energy\n")

        for epoch in tqdm.tqdm(range(iStart, iEnd)):
            print(f"\nEpoch {epoch}/{max_epoch}")

            # Training phase
            model.train()
            train_energy = 0.0
            train_rmsd = 0.0
            n_batches_train = (n_train + batch_size - 1) // batch_size

            # Create random permutation for the entire training set
            train_perm = torch.randperm(n_train)

            for batch in range(n_batches_train):
                print("Epoch", epoch)
                start_idx = batch * batch_size
                end_idx = min((batch + 1) * batch_size, n_train)
                IT = train_perm[start_idx:end_idx]

                # Setup data for this batch
                xyz_dataT = setup_XYZ_dat(all_datT, IT, primary_indices, all_bonds, all_angles)
                curr_batch_size = len(IT)
                f0_bonds_net = f0_bonds.unsqueeze(0).expand(curr_batch_size, -1)
                f0_angles_net = f0_angles.unsqueeze(0).expand(curr_batch_size, -1)
                f0_tors_net = f0_tors.unsqueeze(0).expand(curr_batch_size, -1)

                batch_energy, batch_rmsd = run_model(
                    datT, xyz_dataT, xyzT, model, IT, total_iterations, iteration,
                    fixed_tors_ind, flexibility, curr_batch_size, f0_bonds_net, f0_angles_net, f0_tors_net, train=True,
                    device=device
                )

                if curr_batch_size == batch_size:
                    train_energy += batch_energy
                    train_rmsd += batch_rmsd
                iteration += 1

            train_energy /= n_batches_train
            train_rmsd /= n_batches_train

            # Validation phase
            model.eval()
            val_energy = 0.0
            val_rmsd = 0.0
            n_batches_val = (n_valid + batch_size - 1) // batch_size

            # Create random permutation for the entire validation set
            val_perm = torch.randperm(n_valid)

            with torch.no_grad():
                for batch in range(n_batches_val):
                    start_idx = batch * batch_size
                    end_idx = min((batch + 1) * batch_size, n_valid)
                    IV = val_perm[start_idx:end_idx]

                    # Setup data for this batch
                    xyz_dataV = setup_XYZ_dat(all_datV, IV, primary_indices, all_bonds, all_angles)
                    curr_batch_size = len(IV)
                    f0_bonds_net = f0_bonds.unsqueeze(0).expand(curr_batch_size, -1)
                    f0_angles_net = f0_angles.unsqueeze(0).expand(curr_batch_size, -1)
                    f0_tors_net = f0_tors.unsqueeze(0).expand(curr_batch_size, -1)

                    batch_energy, batch_rmsd = run_model(
                        datV, xyz_dataV, xyzV, model, IV, total_iterations, iteration,
                        fixed_tors_ind, flexibility, curr_batch_size, f0_bonds_net, f0_angles_net, f0_tors_net,
                        train=False, device=device
                    )

                    if curr_batch_size == batch_size:
                        val_energy += batch_energy
                        val_rmsd += batch_rmsd

            val_energy /= n_batches_val
            val_rmsd /= n_batches_val

            print(f"\nEpoch {epoch} Summary:")
            print(f"Train - RMSD: {train_rmsd:.4f}, Energy: {train_energy:.4f}")
            print(f"Valid - RMSD: {val_rmsd:.4f}, Energy: {val_energy:.4f}")

            fout.write(f"{epoch}\t{train_rmsd:.6f}\t{val_rmsd:.6f}\t{train_energy:.6f}\t{val_energy:.6f}\n")

            # Save model checkpoint
            torch.save({
                'model': model.state_dict(),
                'primary_torsions': primary_torsions,
                'angle_differences': angle_differences,
                'primary_indices': primary_indices,
                'all_bonds': all_bonds,
                'all_angles': all_angles,
                'non_primary_bonds': non_primary_bonds,
                'non_primary_angles': non_primary_angles,
                'non_primary_indices': non_primary_indices,
                'fixed_torsions': fixed_tors_ind
            }, out_path_params + str(epoch) + params_file_name)
