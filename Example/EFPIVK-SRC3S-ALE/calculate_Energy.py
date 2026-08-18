import torch
import math

torch.pi = math.pi
import numpy as np
import pytraj as pt
from Atom import PRMTOPParser, Atom

"""
Some future things: 
    - Right now, I have a separate loop to calculate the non-bonded forces, the solvation bond radii, and the
      solvation energy, all of these can be combined into one loop
    - Move all of the parameter parsing and assignments from sander to another function that can be called once
      so the parameters don't have to be re-parsed every iteration
"""


def calculate_energy_test(path, top_path):
    # to test dcd files, function is not used in DL training
    traj = pt.load(path, top_path)
    xyz = traj.xyz
    xyz = torch.as_tensor(xyz, device='cuda', dtype=torch.float32)

    # Calculate mean energy of all frames
    #energies = calculate_energy(xyz, top_path)
    bonded_energy, lennard_jones_energy, electrostatic_energy, solvation_energy = calculate_energy(xyz, top_path)

#    return energies
    return bonded_energy, lennard_jones_energy, electrostatic_energy, solvation_energy

def calculate_energy(coords, top_path, saltcon=150.0):
    parser = PRMTOPParser(top_path)
    prmtop_params = parser.get_prmtop_params()

    bond_energy = calculate_bond_energy(coords, prmtop_params)
    angle_energy = calculate_angle_energy(coords, prmtop_params)
    torsion_energy, lj_14_energy, ee_14_energy = calculate_torsion_energy(coords, prmtop_params)
    lennard_jones_energy, electrostatic_energy = calculate_non_bonded_energy(coords, prmtop_params)
    solvation_energy = calculate_solvation_energy(coords, prmtop_params, saltcon)

#    print("Bond Energy: ", bond_energy)
#    print("\nAngle Energy: ", angle_energy)
#    print("\nTorsion Energy: ", torsion_energy)
#    print("\nNon-Bonded Lennard-Jones Energy: ", lennard_jones_energy)
#    print("\nNon-Bonded Electrostatic Energy: ", electrostatic_energy)
#    print("\n1-4 Lennard Jones Energy: ", lj_14_energy)
#    print("\n1-4 Electrostatic Energy: ", ee_14_energy)
    print("\nTotal Lennard Jones Energy: ", lj_14_energy + lennard_jones_energy)
#    print("\nTotal Electrostatic Energy: ", ee_14_energy + electrostatic_energy)
#    print("\nSolvation Energy: ", solvation_energy)

    total_energy = (bond_energy + angle_energy + torsion_energy + lennard_jones_energy + electrostatic_energy
                    + lj_14_energy + ee_14_energy + solvation_energy)
#    print("\nTotal Energy: ", total_energy)

#    return torch.mean(total_energy)
    #return torch.mean(bond_energy + angle_energy + torsion_energy), torch.mean(lennard_jones_energy+lj_14_energy), torch.mean(electrostatic_energy+ee_14_energy), torch.mean(solvation_energy) 
    return torch.mean(total_energy), torch.mean(bond_energy + angle_energy + torsion_energy), torch.mean(lennard_jones_energy+lj_14_energy+electrostatic_energy+ee_14_energy+solvation_energy) 




def calculate_bond_energy(coords, prmtop_params):
    # _______________________________________________________________________________________________
    # Bond indices are stored with the two atoms of the bond as the first two entries and a pointer
    # to the k and r0 values as the third entry. Energy is calculated according to the
    # formula: k(r-r0)^2 with r being the bond length. Indices are read from prmtop
    # _______________________________________________________________________________________________

    bond_indices = torch.tensor(prmtop_params['bond_indices'], dtype=torch.long, device='cuda')
    bond_k = prmtop_params['bond_k'].clone().detach().to(device='cuda', dtype=torch.float32)
    bond_r0 = prmtop_params['bond_r0'].clone().detach().to(device='cuda', dtype=torch.float32)

    atom1 = bond_indices[:, 0]
    atom2 = bond_indices[:, 1]
    param_idx = bond_indices[:, 2] - 1  # PRMTOP starts with index 1 not 0

    vec = coords[:, atom2] - coords[:, atom1]
    bond_length = vec.pow(2).sum(dim=-1).sqrt()
    energy = bond_k[param_idx] * (bond_length - bond_r0[param_idx]) ** 2
    sum_bond_energy = torch.sum(energy, dim=1)
    return sum_bond_energy


def calculate_angle_energy(coords, prmtop_params):
    # _______________________________________________________________________________________________
    # Angle indices are stored with the three atoms of the angle as the first entries and a pointer
    # to the k and theta0 values as the third entry. Energy is calculated according to the
    # formula: k(theta-theta0)^2 with theta being the bond angle. Indices are read from prmtop
    # _______________________________________________________________________________________________

    angle_indices = torch.tensor(prmtop_params['angle_indices'], dtype=torch.long, device='cuda')
    angle_k = prmtop_params['angle_k'].clone().detach().to(device='cuda', dtype=torch.float32)
    angle_theta0 = prmtop_params['angle_theta0'].clone().detach().to(device='cuda', dtype=torch.float32)

    atom1 = angle_indices[:, 0]
    atom2 = angle_indices[:, 1]
    atom3 = angle_indices[:, 2]
    param_idx = angle_indices[:, 3] - 1  # PRMTOP starts with index 1 not 0

    vec1 = coords[:, atom2] - coords[:, atom1]  # two vector to represent the angle
    vec2 = coords[:, atom3] - coords[:, atom2]

    x1 = (-vec1 * vec2).sum(dim=-1)
    y1 = torch.cross(vec1, vec2).pow(2).sum(axis=-1).sqrt()
    theta = torch.atan2(y1, x1)  # find the angle between two vectors to get bond angle

    angle_energy = angle_k[param_idx] * (theta - angle_theta0[param_idx]) ** 2
    sum_angle_energy = torch.sum(angle_energy, dim=1)
    return sum_angle_energy


def calculate_torsion_energy(coords, prmtop_params):
    # _______________________________________________________________________________________________
    # Torsion indices are stored with the four atoms of the torsion as the first entries and a pointer
    # to the n, v and gamma values as the fourth entry. Energy is calculated according to the
    # formula: v * (1 + cos(n * phi - gamma) with phi being the bond angle. The indices in prmtop
    # are read twice, once with the sign and once without the signs and stored in different
    # lists. The lists without the negative values are used to find the torsional energy using the
    # above formula, and the list with negative values is used to find the torsion angles that are
    # not considered when calculating the 1-4 electrostatic and 1-4 van der wall energies that are
    # also calculated in this function.
    # _______________________________________________________________________________________________

    torsion_indices = torch.tensor(prmtop_params['torsion_indices'], dtype=torch.long, device='cuda')
    torsion_indices_with_sign = torch.tensor(prmtop_params['torsion_indices_with_sign'], dtype=torch.long,
                                             device='cuda')
    torsion_v = prmtop_params['torsion_v'].clone().detach().to(device='cuda', dtype=torch.float32)
    torsion_n = prmtop_params['torsion_n'].clone().detach().to(device='cuda', dtype=torch.float32)
    torsion_gamma = prmtop_params['torsion_gamma'].clone().detach().to(device='cuda', dtype=torch.float32)
    atom_type_index = prmtop_params['atom_type_index'].clone().detach().to(device='cuda', dtype=torch.long)
    nonbonded_parm_index = prmtop_params['nonbonded_parm_index'].clone().detach().to(device='cuda', dtype=torch.long)
    lennard_jones_acoef = prmtop_params['lennard_jones_acoef'].clone().detach().to(device='cuda', dtype=torch.float32)
    lennard_jones_bcoef = prmtop_params['lennard_jones_bcoef'].clone().detach().to(device='cuda', dtype=torch.float32)
    charges = prmtop_params['charges'].clone().detach().to(device='cuda', dtype=torch.float32)
    scnb = prmtop_params['scnb_scale_factor'].clone().detach().to(device='cuda', dtype=torch.float32)
    scee = prmtop_params['scee_scale_factor'].clone().detach().to(device='cuda', dtype=torch.float32)
    scnb = torch.where(scnb == 0, torch.tensor(-1.0, device='cuda', dtype=torch.float32), scnb)
    scee = torch.where(scee == 0, torch.tensor(-1.0, device='cuda', dtype=torch.float32), scee)

    MIN_BOND_LENGTH = 1 
    MIN = -1e6
    MAX = 1e6
    ntypes = int(len(nonbonded_parm_index) ** 0.5)  # From prmtop

    param_idx = torsion_indices[:, 4] - 1  # PRMTOP starts with index 1 not 0

    atom1 = torsion_indices[:, 0]
    atom2 = torsion_indices[:, 1]
    atom3 = torsion_indices[:, 2]
    atom4 = torsion_indices[:, 3]

    p1 = coords[:, atom1]
    p2 = coords[:, atom2]
    p3 = coords[:, atom3]
    p4 = coords[:, atom4]

    v21 = p1 - p2
    v32 = p2 - p3
    v43 = p3 - p4

    n1 = torch.cross(-v21, v32)
    n2 = torch.cross(-v32, v43)


    xp = torch.cross(n1, n2)
    x2 = (n1 * n2).sum(dim=-1)
    y2 = (xp * v32).sum(dim=-1) / torch.norm(v32, dim=-1)
    assert not torch.isnan(y2).any(), "NaN found in y2"
    assert not torch.isnan(x2).any(), "NaN found in x2"

    phi = torch.atan2(y2, x2)
    phi = ((phi + torch.pi) % (2 * torch.pi)) - torch.pi
    print(phi.grad)

    assert not torch.isnan(phi).any(), "NaN found in phi"

    torsion_energy = torsion_v[param_idx] * (1 + torch.cos(torsion_n[param_idx] * phi - torsion_gamma[param_idx]))

    assert not torch.isnan(torsion_energy).any(), "NaN found in torsion_energy"

    # _______________________________________________________________________________________________
    # This section is used to calculate the 1-4 van_der_wall and 1-4 electrostatic energies. Torsion
    # angles that have a negative third of fourth atom index in prmtop are ignored. The method
    # used to calculate the 1-4 van der wall and 1-4 electrostatic energies are detailed in the
    # "calculate_non_bonded_energy" function. The 1-4 electrostatic and 1-4 van der wall
    # interactions are scaled down by values given in the scee and scnb lists respectively as defined
    # in the prmtop file
    # _______________________________________________________________________________________________

    valid_14_mask = (torsion_indices_with_sign[:, 2] > 0) & (torsion_indices_with_sign[:, 3] > 0)

    abcoef_indx = nonbonded_parm_index[ntypes * (atom_type_index[atom1] - 1) + atom_type_index[atom4] - 1] - 1
    a_ij = lennard_jones_acoef[abcoef_indx]
    b_ij = lennard_jones_bcoef[abcoef_indx]
    qiqj = charges[atom1] * charges[atom4]

    bond_vec = coords[:, atom1] - coords[:, atom4]
    bond_length_squared = torch.sum(bond_vec ** 2, dim=-1)
    bond_length_squared = torch.clamp(bond_length_squared, min=MIN_BOND_LENGTH ** 2)
    bond_length = torch.sqrt(bond_length_squared)
    bond_length = torch.clamp(bond_length, min=MIN_BOND_LENGTH)

    inv_rij_12 = 1.0 / torch.pow(bond_length, 12)
    inv_rij_6 = 1.0 / torch.pow(bond_length, 6)
    inv_rij_6 = torch.clamp(inv_rij_6, max=MAX)
    inv_rij_12 = torch.clamp(inv_rij_12, max=MAX)


    lj_14 = torch.where(valid_14_mask,
                        (a_ij * inv_rij_12 - b_ij * inv_rij_6) / scnb[param_idx],
                        torch.zeros_like(bond_length)
                        )

    ee_14 = torch.where(valid_14_mask,
                        (qiqj / bond_length) / scee[param_idx],
                        torch.zeros_like(bond_length)
                        )


    sum_lj_14 = torch.sum(lj_14, dim=1)
    sum_ee_14 = torch.sum(ee_14, dim=1)
    sum_torsion_energy = torch.sum(torsion_energy, dim=1)

    return sum_torsion_energy, sum_lj_14, sum_ee_14


def calculate_non_bonded_energy(coords, prmtop_params):
    # _______________________________________________________________________________________________
    # This section is used to calculate the Van_der_wall and Electrostatic energies. Note that unlike
    # Sander, this section does not recalculate the exclusion list based on a cutoff and instead
    # reads the exclusion list directly from prmtop. Due to the large number of possible interactions
    # the calculations for this section are done in batches of frames. The number of frames that are
    # calculated in parallel can be tuned based on the device capacity
    # _______________________________________________________________________________________________

    atom_type_index = prmtop_params['atom_type_index'].clone().detach().to(device='cuda', dtype=torch.long)
    number_excluded_atoms = prmtop_params['number_excluded_atoms'].clone().detach().to(device='cuda', dtype=torch.long)
    nonbonded_parm_index = prmtop_params['nonbonded_parm_index'].clone().detach().to(device='cuda', dtype=torch.long)
    excluded_atoms_list = prmtop_params['excluded_atoms_list'].clone().detach().to(device='cuda', dtype=torch.long)
    lennard_jones_acoef = prmtop_params['lennard_jones_acoef'].clone().detach().to(device='cuda', dtype=torch.float32)
    lennard_jones_bcoef = prmtop_params['lennard_jones_bcoef'].clone().detach().to(device='cuda', dtype=torch.float32)
    charges = prmtop_params['charges'].clone().detach().to(device='cuda', dtype=torch.float32)
    MIN = -1e6
    MAX = 1e6
    MIN_BOND_LENGTH = 1 

    n_atoms = coords.shape[1]
    n_frames = coords.shape[0]
    ntypes = int(len(nonbonded_parm_index) ** 0.5)

    excluded_mask = torch.ones((n_atoms, n_atoms), device='cuda', dtype=torch.bool)
    exclusion_start = 0

    # _______________________________________________________________________________________________
    # Number_Excluded_Atoms has the number of atoms excluded in Lennard Jones and Electrostatic
    # interactions for each atom because they form a bond, angle or torsion with the atom.
    # The exact atom indices of the excluded atoms are in the excluded_atoms_list
    # _______________________________________________________________________________________________

    for i in range(n_atoms):  # parse through prmtop file to set excluded pairs to False
        exclusion_end = exclusion_start + number_excluded_atoms[i]
        excluded_atoms = excluded_atoms_list[exclusion_start:exclusion_end]
        excluded_mask[i, excluded_atoms[excluded_atoms != 0] - 1] = False
        exclusion_start = exclusion_end

    # Matrix of all interactions with self interactions and repeat interactions removed
    i, j = torch.triu_indices(n_atoms, n_atoms, offset=1, device='cuda')

    # Ignore excluded interactions
    valid_pairs = excluded_mask[i, j]
    i = i[valid_pairs]
    j = j[valid_pairs]

    # _______________________________________________________________________________________________
    # non_bonded_parm index is a pointer into lennard_jones_acoef and lennard_jones_bcoef which
    # will yield the A and B value for each non-bonded pair interaction. Using these parameters,
    # the Lennard-Jones energy can be calculated using the formula: A/(rij)^12 - B/(rij)^6. The
    # Electrostatic interactions are calculated using the same interaction pair indices and the
    # energy of the interaction is given by the formula: qi*qj/rij
    # _______________________________________________________________________________________________

    abcoef_indx = nonbonded_parm_index[ntypes * (atom_type_index[i] - 1) + atom_type_index[j] - 1] - 1

    a_ij = lennard_jones_acoef[abcoef_indx]
    b_ij = lennard_jones_bcoef[abcoef_indx]

    sum_vdw_energy = torch.zeros(n_frames, device='cuda', dtype=torch.float32)
    sum_ee_energy = torch.zeros(n_frames, device='cuda', dtype=torch.float32)

    batch_size = 1024  # Frame batch size, adjust this value based on your device capacity
    for batch_start in range(0, n_frames, batch_size):
        batch_end = min(batch_start + batch_size, n_frames)
        batch_coords = coords[batch_start:batch_end]
        curr_batch = batch_end - batch_start

        bond_vec = batch_coords[:, j] - batch_coords[:, i]
        qiqj = (charges[i] * charges[j]).unsqueeze(0).expand(curr_batch, -1)
        bond_length_squared = torch.sum(bond_vec ** 2, dim=-1)
        bond_length_squared = torch.clamp(bond_length_squared, min=MIN_BOND_LENGTH ** 2)
        bond_length = torch.sqrt(bond_length_squared)
        bond_length = torch.clamp(bond_length, min=MIN_BOND_LENGTH)

        inv_rij_12 = 1.0 / torch.pow(bond_length, 12)
        inv_rij_6 = 1.0 / torch.pow(bond_length, 6)
        #inv_rij_12 = 1.0 / torch.pow(bond_length, 6)
        #inv_rij_6 = 1.0 / torch.pow(bond_length, 3)
        inv_rij_6 = torch.clamp(inv_rij_6, max=MAX)
        inv_rij_12 = torch.clamp(inv_rij_12, max=MAX)

        
        E_lj = a_ij * inv_rij_12 - b_ij * inv_rij_6
        sign_E_lj = torch.sign(E_lj.sum(dim=1))
        abs_E_lj = torch.abs(E_lj.sum(dim=1))
        sqrt_abs_E_lj = torch.pow(abs_E_lj,0.5)
        signed_sqrt_E_lj = sign_E_lj * sqrt_abs_E_lj

        #E_lj = a_ij * inv_rij_12 - b_ij * inv_rij_6
        
        E_ee = qiqj / bond_length
        #sum_vdw_energy[batch_start:batch_end] = signed_sqrt_E_lj 
        sum_vdw_energy[batch_start:batch_end] = E_lj.sum(dim=1)
        sum_ee_energy[batch_start:batch_end] = E_ee.sum(dim=1)

    return sum_vdw_energy, sum_ee_energy


def calculate_solvation_energy(coords, prmtop_params, saltcon):
    # ___________________________________________________________________________________
    # This section calculates the solvation energy using the GB model flag igb = 8.
    # All parameters are either from the prmtop file or sander. Most variables that are
    # from sander are taken from the file mdread1 and mdread2. GBNECKCUT is taken from
    # gbhead. The screen parameters from prmtop are ignored and recalculated.
    # Note that the cutoff used to ignore atom interactions used by sander in the
    # calculation of the VDW Integral is not employed. Again, due to the large number of
    # interactions, this section in performed in batches
    # ___________________________________________________________________________________

    charges = prmtop_params['charges'].clone().detach().to(device='cuda', dtype=torch.float32)
    rborn = prmtop_params['rborn'].clone().detach().to(device='cuda', dtype=torch.float32)
    names = prmtop_params['atom_names']
    n_atoms = coords.shape[1]

    n_frames = coords.shape[0]
    MIN = 1e-5
    MAX = 1e5
    offset = 0.195141
    rborn_offset = rborn.add(-offset)  # subtract the offset
    kappa = 0.73 * (0.10806 * saltcon / 1000) ** 0.5  # salt concentration conversion
    extdiel = 78.5  # dielectric constant of water
    gbneckcut = 6.8  # "neck" - diameter at which a solvent molecule is considered unable to go through
    gbneckscale = 0.826836
    gbalphaH = 0.788440
    gbbetaH = 0.798699
    gbgammaH = 0.437334
    gbalphaC = 0.733756
    gbbetaC = 0.506378
    gbgammaC = 0.205844
    gbalphaN = 0.503364
    gbbetaN = 0.316828
    gbgammaN = 0.192915
    gbalphaOS = 0.867814
    gbbetaOS = 0.876635
    gbgammaOS = 0.387882
    gbalphaP = 0.418365
    gbbetaP = 0.290054
    gbgammaP = 0.1064245  # P parameters are not optimized for protein
    Sh = 1.425952
    Sc = 1.058554
    Sn = 0.733599
    So = 1.061039
    Ss = -0.703469
    Sp = 0.5  # P parameters are not optimized for proteins
    atom_scale = torch.zeros(n_atoms, device='cuda', dtype=torch.float32)
    gb_alpha = torch.zeros(n_atoms, device='cuda', dtype=torch.float32)
    gb_beta = torch.zeros(n_atoms, device='cuda', dtype=torch.float32)
    gb_gamma = torch.zeros(n_atoms, device='cuda', dtype=torch.float32)

    # _____________________________________________________________________________________________________
    # Set the atom specific scale factor, alpha, beta, and gamma values for each
    # specific atom. The scale of each atom is used to calculate the new screen
    # parameters that replace the ones from the prmtop file according to this
    # equation: [screen = scale * (rborn - offset)]. The alpha, beta, and gamma values
    # will eventually be used to find the effective radii according to the equation
    # effective_radius-1 = (radius - offset)-1 - radius-1 * tanh(alpha * PSI - beta * PSI^2 * gamma PSI^3)
    # with PSI being equal to (radius - offset) * Ivdw, with Ivdw being the de-screening contribution
    # _____________________________________________________________________________________________________

    names_np = np.array(names)
    H_condition = np.char.find(names_np, 'H') >= 0
    C_condition = np.char.find(names_np, 'C') >= 0
    N_condition = np.char.find(names_np, 'N') >= 0
    O_condition = np.char.find(names_np, 'O') >= 0
    S_condition = np.char.find(names_np, 'S') >= 0
    P_condition = np.char.find(names_np, 'P') >= 0

    atom_scale[H_condition] = Sh
    atom_scale[C_condition] = Sc
    atom_scale[N_condition] = Sn
    atom_scale[O_condition] = So
    atom_scale[S_condition] = Ss
    atom_scale[P_condition] = Sp

    gb_alpha[H_condition] = gbalphaH
    gb_alpha[C_condition] = gbalphaC
    gb_alpha[N_condition] = gbalphaN
    gb_alpha[O_condition] = gbalphaOS
    gb_alpha[S_condition] = gbalphaOS
    gb_alpha[P_condition] = gbalphaP

    gb_beta[H_condition] = gbbetaH
    gb_beta[C_condition] = gbbetaC
    gb_beta[N_condition] = gbbetaN
    gb_beta[O_condition] = gbbetaOS
    gb_beta[S_condition] = gbbetaOS
    gb_beta[P_condition] = gbbetaP

    gb_gamma[H_condition] = gbgammaH
    gb_gamma[C_condition] = gbgammaC
    gb_gamma[N_condition] = gbgammaN
    gb_gamma[O_condition] = gbgammaOS
    gb_gamma[S_condition] = gbgammaOS
    gb_gamma[P_condition] = gbgammaP

    fs = atom_scale * rborn_offset  # recalculate screening parameters
    fsmax = torch.max(fs)

    # _____________________________________________________________________________________________________
    # Lookup table and index for each atom in gbneck correction, taken from gbneck.h file of sander
    # _____________________________________________________________________________________________________

    neck_index = torch.round((rborn - 1.0) * 20.0).clone().detach().to(device='cuda', dtype=torch.long)

    neckMaxPos = torch.tensor([
        [2.26685, 2.32548, 2.38397, 2.44235, 2.50057, 2.55867, 2.61663, 2.67444, 2.73212, 2.78965, 2.84705, 2.9043,
         2.96141, 3.0184, 3.07524, 3.13196, 3.18854, 3.24498, 3.30132, 3.35752, 3.4136],
        [2.31191, 2.37017, 2.4283, 2.48632, 2.5442, 2.60197, 2.65961, 2.71711, 2.77449, 2.83175, 2.88887, 2.94586,
         3.00273, 3.05948, 3.1161, 3.1726, 3.22897, 3.28522, 3.34136, 3.39738, 3.45072],
        [2.35759, 2.41549, 2.47329, 2.53097, 2.58854, 2.646, 2.70333, 2.76056, 2.81766, 2.87465, 2.93152, 2.98827,
         3.0449, 3.10142, 3.15782, 3.21411, 3.27028, 3.32634, 3.3823, 3.43813, 3.49387],
        [2.4038, 2.46138, 2.51885, 2.57623, 2.63351, 2.69067, 2.74773, 2.80469, 2.86152, 2.91826, 2.97489, 3.0314,
         3.08781, 3.1441, 3.20031, 3.25638, 3.31237, 3.36825, 3.42402, 3.4797, 3.53527],
        [2.45045, 2.50773, 2.56492, 2.62201, 2.679, 2.7359, 2.7927, 2.8494, 2.90599, 2.9625, 3.0189, 3.07518, 3.13138,
         3.18748, 3.24347, 3.29937, 3.35515, 3.41085, 3.46646, 3.52196, 3.57738],
        [2.4975, 2.5545, 2.61143, 2.66825, 2.72499, 2.78163, 2.83818, 2.89464, 2.95101, 3.00729, 3.06346, 3.11954,
         3.17554, 3.23143, 3.28723, 3.34294, 3.39856, 3.45409, 3.50952, 3.56488, 3.62014],
        [2.54489, 2.60164, 2.6583, 2.71488, 2.77134, 2.8278, 2.88412, 2.94034, 2.9965, 3.05256, 3.10853, 3.16442,
         3.22021, 3.27592, 3.33154, 3.38707, 3.44253, 3.49789, 3.55316, 3.60836, 3.66348],
        [2.59259, 2.6491, 2.70553, 2.76188, 2.81815, 2.87434, 2.93044, 2.98646, 3.04241, 3.09827, 3.15404, 3.20974,
         3.26536, 3.32089, 3.37633, 3.4317, 3.48699, 3.54219, 3.59731, 3.65237, 3.70734],
        [2.64054, 2.69684, 2.75305, 2.80918, 2.86523, 2.92122, 2.97712, 3.03295, 3.0887, 3.14437, 3.19996, 3.25548,
         3.31091, 3.36627, 3.42156, 3.47677, 3.5319, 3.58695, 3.64193, 3.69684, 3.75167],
        [2.68873, 2.74482, 2.80083, 2.85676, 2.91262, 2.96841, 3.02412, 3.07976, 3.13533, 3.19082, 3.24623, 3.30157,
         3.35685, 3.41205, 3.46718, 3.52223, 3.57721, 3.63213, 3.68696, 3.74174, 3.79644],
        [2.73713, 2.79302, 2.84884, 2.90459, 2.96027, 3.01587, 3.0714, 3.12686, 3.18225, 3.23757, 3.29282, 3.34801,
         3.40313, 3.45815, 3.51315, 3.56805, 3.6229, 3.67767, 3.73237, 3.78701, 3.84159],
        [2.78572, 2.84143, 2.89707, 2.95264, 3.00813, 3.06356, 3.11892, 3.17422, 3.22946, 3.28462, 3.33971, 3.39474,
         3.44971, 3.5046, 3.55944, 3.61421, 3.66891, 3.72356, 3.77814, 3.83264, 3.8871],
        [2.83446, 2.89, 2.94547, 3.00088, 3.05621, 3.11147, 3.16669, 3.22183, 3.27689, 3.33191, 3.38685, 3.44174,
         3.49656, 3.55132, 3.60602, 3.66066, 3.71523, 3.76975, 3.82421, 3.8786, 3.93293],
        [2.88335, 2.93873, 2.99404, 3.04929, 3.10447, 3.15959, 3.21464, 3.26963, 3.32456, 3.37943, 3.43424, 3.48898,
         3.54366, 3.5983, 3.65287, 3.70737, 3.76183, 3.81622, 3.87056, 3.92484, 3.97905],
        [2.93234, 2.9876, 3.04277, 3.09786, 3.15291, 3.20787, 3.26278, 3.31764, 3.37242, 3.42716, 3.48184, 3.53662,
         3.591, 3.64551, 3.69995, 3.75435, 3.80867, 3.86295, 3.91718, 3.97134, 4.02545],
        [2.98151, 3.0366, 3.09163, 3.14659, 3.20149, 3.25632, 3.3111, 3.36581, 3.42047, 3.47507, 3.52963, 3.58411,
         3.63855, 3.69293, 3.74725, 3.80153, 3.85575, 3.90991, 3.96403, 4.01809, 4.07211],
        [3.03074, 3.08571, 3.14061, 3.19543, 3.25021, 3.30491, 3.35956, 3.41415, 3.46869, 3.52317, 3.57759, 3.63196,
         3.68628, 3.74054, 3.79476, 3.84893, 3.90303, 3.95709, 4.01111, 4.06506, 4.11897],
        [3.08008, 3.13492, 3.1897, 3.2444, 3.29905, 3.35363, 3.40815, 3.46263, 3.51704, 3.57141, 3.62572, 3.67998,
         3.73418, 3.78834, 3.84244, 3.8965, 3.95051, 4.00447, 4.05837, 4.11224, 4.16605],
        [3.12949, 3.18422, 3.23888, 3.29347, 3.348, 3.40247, 3.45688, 3.51124, 3.56554, 3.6198, 3.674, 3.72815,
         3.78225, 3.83629, 3.8903, 3.94425, 3.99816, 4.05203, 4.10583, 4.15961, 4.21333],
        [3.17899, 3.23361, 3.28815, 3.34264, 3.39706, 3.45142, 3.50571, 3.55997, 3.61416, 3.66831, 3.72241, 3.77645,
         3.83046, 3.8844, 3.93831, 3.99216, 4.04598, 4.09974, 4.15347, 4.20715, 4.26078],
        [3.22855, 3.28307, 3.33751, 3.39188, 3.4462, 3.50046, 3.55466, 3.6088, 3.6629, 3.71694, 3.77095, 3.82489,
         3.8788, 3.93265, 3.98646, 4.04022, 4.09395, 4.14762, 4.20126, 4.25485, 4.3084]
    ], device='cuda', dtype=torch.float32)

    neckMaxVal = torch.tensor([
        [0.0381511, 0.0338587, 0.0301776, 0.027003, 0.0242506, 0.0218529, 0.0197547, 0.0179109, 0.0162844, 0.0148442,
         0.0135647, 0.0124243, 0.0114047, 0.0104906, 0.00966876, 0.008928, 0.0082587, 0.00765255, 0.00710237,
         0.00660196, 0.00614589],
        [0.0396198, 0.0351837, 0.0313767, 0.0280911, 0.0252409, 0.0227563, 0.0205808, 0.0186681, 0.0169799, 0.0154843,
         0.014155, 0.0129696, 0.0119094, 0.0109584, 0.0101031, 0.00933189, 0.0086348, 0.00800326, 0.00742986,
         0.00690814, 0.00643255],
        [0.041048, 0.0364738, 0.0325456, 0.0291532, 0.0262084, 0.0236399, 0.0213897, 0.0194102, 0.0176622, 0.0161129,
         0.0147351, 0.0135059, 0.0124061, 0.0114192, 0.0105312, 0.00973027, 0.00900602, 0.00834965, 0.0077535,
         0.00721091, 0.00671609],
        [0.0424365, 0.0377295, 0.0336846, 0.0301893, 0.0271533, 0.0245038, 0.0221813, 0.0201371, 0.018331, 0.0167295,
         0.0153047, 0.014033, 0.0128946, 0.0118727, 0.0109529, 0.0101229, 0.00937212, 0.00869147, 0.00807306,
         0.00751003, 0.00699641],
        [0.0437861, 0.0389516, 0.0347944, 0.0311998, 0.0280758, 0.0253479, 0.0229555, 0.0208487, 0.0189864, 0.0173343,
         0.0158637, 0.0145507, 0.0133748, 0.0123188, 0.0113679, 0.0105096, 0.0097329, 0.00902853, 0.00838835,
         0.00780533, 0.0072733],
        [0.0450979, 0.0401406, 0.0358753, 0.0321851, 0.0289761, 0.0261726, 0.0237125, 0.0215451, 0.0196282, 0.017927,
         0.0164121, 0.0150588, 0.0138465, 0.0127573, 0.0117761, 0.0108902, 0.0100882, 0.00936068, 0.00869923,
         0.00809665, 0.00754661],
        [0.0463729, 0.0412976, 0.0369281, 0.0331456, 0.0298547, 0.026978, 0.0244525, 0.0222264, 0.0202567, 0.0185078,
         0.0169498, 0.0155575, 0.0143096, 0.0131881, 0.0121775, 0.0112646, 0.010438, 0.00968781, 0.00900559, 0.00838388,
         0.00781622],
        [0.0476123, 0.0424233, 0.0379534, 0.034082, 0.0307118, 0.0277645, 0.0251757, 0.0228927, 0.0208718, 0.0190767,
         0.0174768, 0.0160466, 0.0147642, 0.0136112, 0.0125719, 0.0116328, 0.0107821, 0.0100099, 0.00930735, 0.00866695,
         0.00808206],
        [0.0488171, 0.0435186, 0.038952, 0.0349947, 0.0315481, 0.0285324, 0.0258824, 0.0235443, 0.0214738, 0.0196339,
         0.0179934, 0.0165262, 0.0152103, 0.0140267, 0.0129595, 0.0119947, 0.0111206, 0.0103268, 0.00960445, 0.00894579,
         0.00834405],
        [0.0499883, 0.0445845, 0.0399246, 0.0358844, 0.032364, 0.0292822, 0.0265729, 0.0241815, 0.0220629, 0.0201794,
         0.0184994, 0.0169964, 0.0156479, 0.0144345, 0.0133401, 0.0123504, 0.0114534, 0.0106386, 0.00989687, 0.00922037,
         0.00860216],
        [0.0511272, 0.0456219, 0.040872, 0.0367518, 0.0331599, 0.0300142, 0.0272475, 0.0248045, 0.0226392, 0.0207135,
         0.0189952, 0.0174574, 0.0160771, 0.0148348, 0.0137138, 0.0126998, 0.0117805, 0.0109452, 0.0101846, 0.00949067,
         0.00885636],
        [0.0522348, 0.0466315, 0.0417948, 0.0375973, 0.0339365, 0.030729, 0.0279067, 0.0254136, 0.023203, 0.0212363,
         0.0194809, 0.0179092, 0.016498, 0.0152275, 0.0140807, 0.013043, 0.012102, 0.0112466, 0.0104676, 0.00975668,
         0.00910664],
        [0.0533123, 0.0476145, 0.042694, 0.0384218, 0.0346942, 0.0314268, 0.0285507, 0.026009, 0.0237547, 0.0217482,
         0.0199566, 0.018352, 0.0169108, 0.0156128, 0.0144408, 0.0133801, 0.0124179, 0.011543, 0.010746, 0.0100184,
         0.00935302],
        [0.0543606, 0.0485716, 0.04357, 0.0392257, 0.0354335, 0.0321082, 0.02918, 0.0265913, 0.0242943, 0.0222492,
         0.0204225, 0.0187859, 0.0173155, 0.0159908, 0.0147943, 0.0137111, 0.0127282, 0.0118343, 0.0110197, 0.0102759,
         0.00959549],
        [0.0553807, 0.0495037, 0.0444239, 0.0400097, 0.0361551, 0.0327736, 0.0297949, 0.0271605, 0.0248222, 0.0227396,
         0.0208788, 0.0192111, 0.0177122, 0.0163615, 0.0151413, 0.0140361, 0.013033, 0.0121206, 0.0112888, 0.0105292,
         0.00983409],
        [0.0563738, 0.0504116, 0.0452562, 0.0407745, 0.0368593, 0.0334235, 0.0303958, 0.0277171, 0.0253387, 0.0232197,
         0.0213257, 0.0196277, 0.0181013, 0.0167252, 0.0154817, 0.0143552, 0.0133325, 0.0124019, 0.0115534, 0.0107783,
         0.0100688],
        [0.0573406, 0.0512963, 0.0460676, 0.0415206, 0.0375468, 0.0340583, 0.030983, 0.0282614, 0.0258441, 0.0236896,
         0.0217634, 0.020036, 0.0184826, 0.017082, 0.0158158, 0.0146685, 0.0136266, 0.0126783, 0.0118135, 0.0110232,
         0.0102998],
        [0.0582822, 0.0521584, 0.0468589, 0.0422486, 0.038218, 0.0346784, 0.0315571, 0.0287938, 0.0263386, 0.0241497,
         0.0221922, 0.0204362, 0.0188566, 0.0174319, 0.0161437, 0.0149761, 0.0139154, 0.0129499, 0.0120691, 0.0112641,
         0.0105269],
        [0.0591994, 0.0529987, 0.0476307, 0.042959, 0.0388734, 0.0352843, 0.0321182, 0.0293144, 0.0268225, 0.0246002,
         0.0226121, 0.0208283, 0.0192232, 0.0177751, 0.0164654, 0.015278, 0.0141991, 0.0132167, 0.0123204, 0.0115009,
         0.0107504],
        [0.0600932, 0.053818, 0.0483836, 0.0436525, 0.0395136, 0.0358764, 0.0326669, 0.0298237, 0.0272961, 0.0250413,
         0.0230236, 0.0212126, 0.0195826, 0.0181118, 0.0167811, 0.0155744, 0.0144778, 0.0134789, 0.0125673, 0.0117338,
         0.0109702],
        [0.0609642, 0.0546169, 0.0491183, 0.0443295, 0.0401388, 0.036455, 0.0332033, 0.030322, 0.0277596, 0.0254732,
         0.0234266, 0.0215892, 0.0199351, 0.018442, 0.0170909, 0.0158654, 0.0147514, 0.0137365, 0.0128101, 0.0119627,
         0.0111863]
    ], device='cuda', dtype=torch.float32)

    neckMaxPos = torch.t(neckMaxPos)
    neckMaxVal = torch.t(neckMaxVal)

    batch_size = 1024

    # ____________________________________________________________________________________________
    # Loop to calculate effective radii of all atoms of all frames. This loop
    # implements equations 9, 10 and 11 from J. Phys. Chem. 1996, 100, 51, 19824–19839.
    # to analytically find the de-screening contribution using a pair-wise approximation
    # Note that sander does not use the same method to find the de-screening contribution
    # and instead uses the equations found in the appendix of J. Mol. Biol. (1990) 216, 1045-1066.
    # However, they note that "Aside from the scaling of the radii, this is the same approach
    # developed in Schaefer and Froemmel, JMB 216:1045 (1990)."
    # _____________________________________________________________________________________________

    r_eff = torch.zeros((n_frames, n_atoms), device='cuda', dtype=torch.float32)
    for batch_start in range(0, n_frames, batch_size):
        batch_end = min(batch_start + batch_size, n_frames)
        batch_size_current = batch_end - batch_start

        rij_squared = torch.zeros((batch_size_current, n_atoms, n_atoms), device='cuda', dtype=torch.float32)
        inv_rborn_offset = 1.0 / rborn_offset
        rborn_offset_batch = rborn_offset.unsqueeze(0).expand(batch_size_current, -1)
        rborn_batch = rborn_offset.unsqueeze(0).expand(batch_size_current, -1)
        sjpj_batch = fs.unsqueeze(0).expand(n_atoms, -1)
        ri_batch = rborn_offset.unsqueeze(1).expand(-1, n_atoms)
        sjpj_batch_squared = sjpj_batch ** 2

        # _________________________________________________________________________________
        # Step 1: Calculate all pairwise atomic distances
        # _________________________________________________________________________________

        coords_batch = coords[batch_start:batch_end]

        coords_i = coords_batch.unsqueeze(2)
        coords_j = coords_batch.unsqueeze(1)

        diff = coords_i - coords_j
        rij_squared = torch.sum(diff ** 2, dim=-1)
        rij_squared = torch.clamp(rij_squared, min=MIN ** 2, max=MAX ** 2)
        rij = torch.sqrt(rij_squared)
        rij = torch.clamp(rij, min=MIN, max=MAX)

        # _________________________________________________________________________________
        # Step 2: Find upper and lower bound Ukk and Lkk as described by equation 10 and 11
        # of J. Phys. Chem. 1996, 100, 51, 19824–19839.
        # _________________________________________________________________________________

        sjpj = sjpj_batch.unsqueeze(0).expand(batch_size_current, -1, -1)
        pi = ri_batch.unsqueeze(0).expand(batch_size_current, -1, -1)

        condition1 = (pi < rij + sjpj)
        condition2 = (rij + sjpj <= pi)

        u_ij = torch.where(condition1, rij + sjpj,
                           torch.where(condition2, torch.ones_like(rij), torch.zeros_like(rij)))

        u_ij = torch.where(rij != MIN, u_ij, torch.zeros_like(rij))
        u_ij = torch.clamp(u_ij, min=MIN, max=MAX)

        condition1 = (pi <= rij - sjpj)
        condition2 = (rij - sjpj <= pi) & (pi < rij + sjpj)
        condition3 = (pi < rij - sjpj)

        l_ij = torch.where(condition1, rij - sjpj,
                           torch.where(condition2, pi,
                                       torch.where(condition3, torch.ones_like(rij), torch.zeros_like(rij))))
        l_ij = torch.where(rij != 0, l_ij, torch.zeros_like(rij))
        l_ij = torch.clamp(l_ij, min=MIN, max=MAX)

        # _________________________________________________________________________________
        # Step 3: Calculate de-screening contribution for each atom as described in
        # equation 9 of J. Phys. Chem. 1996, 100, 51, 19824–19839.
        # _________________________________________________________________________________

        inv_uij = torch.where(u_ij != MIN, 1.0 / u_ij, torch.zeros_like(u_ij))
        inv_uij = torch.clamp(inv_uij, min=MIN, max=MAX)
        inv_lij = torch.where(l_ij != MIN, 1.0 / l_ij, torch.zeros_like(l_ij))
        inv_lij = torch.clamp(inv_lij, min=MIN, max=MAX)
        inv_rij = torch.where(rij != MIN, 1.0 / rij, torch.zeros_like(rij))
        inv_rij = torch.clamp(inv_rij, min=MIN, max=MAX)
        inv_uij_squared = inv_uij ** 2
        inv_uij_squared = torch.clamp(inv_uij_squared, min=MIN, max=MAX)
        inv_lij_squared = inv_lij ** 2
        inv_lij_squared = torch.clamp(inv_lij_squared, min=MIN, max=MAX)
        ln_lij_over_uij = torch.where(rij != MIN, torch.log(l_ij * inv_uij), torch.zeros_like(rij))
        sjpj_over_4rij = 0.25 * sjpj * inv_rij
        sjpj_over_4rij = torch.clamp(sjpj_over_4rij, min=MIN, max=MAX)

        vdw_integral_to_sum = 0.5 * torch.where(rij != MIN, (inv_lij - inv_uij
                                                     + 0.25 * rij * (inv_uij_squared - inv_lij_squared)
                                                     + 0.5 * inv_rij * ln_lij_over_uij
                                                     + sjpj_over_4rij * (inv_lij_squared - inv_uij_squared)),
                                          torch.zeros_like(rij))
        # _________________________________________________________________________________
        # Neck Correction implemented using equation 5 and 6 of J. Chem. Theory Comput.
        # 2007, 3, 1, 156–169
        # _________________________________________________________________________________

        gbneck_cutoff = (rij < rborn_offset.unsqueeze(1) + rborn_offset.unsqueeze(0) + gbneckcut)

        # Calculate mdist, mdist2, mdist3, and mdist6
        mdist = rij - neckMaxPos[neck_index.unsqueeze(1), neck_index.unsqueeze(0)]
        mdist2 = mdist * mdist
        mdist3 = mdist2 * mdist
        mdist6 = mdist3 * mdist3

        neck = neckMaxVal[neck_index.unsqueeze(1), neck_index.unsqueeze(0)] / (1 + mdist2 + 0.3 * mdist6)

        vdw_integral_to_sum = torch.where(rij != MIN, vdw_integral_to_sum, torch.zeros_like(vdw_integral_to_sum))
        vdw_integral_to_sum = torch.where(gbneck_cutoff,
                                          vdw_integral_to_sum + gbneckscale * neck,
                                          vdw_integral_to_sum)

        vdw_integral = torch.sum(vdw_integral_to_sum, dim=-1)

        # _______________________________________________________________________________________________
        # Step 4: Calculate effective radius, using equation 4 of J. Chem. Theory Comput. 2007, 3, 156-169
        # The parameters for alpha, beta, and gamma are taken from sander
        # _______________________________________________________________________________________________

        gb_alpha_batch = gb_alpha.unsqueeze(0).expand(batch_size_current, -1)
        gb_beta_batch = gb_beta.unsqueeze(0).expand(batch_size_current, -1)
        gb_gamma_batch = gb_gamma.unsqueeze(0).expand(batch_size_current, -1)
        PSI = vdw_integral * rborn_offset_batch
        PSI = torch.clamp(PSI, min=MIN, max=MAX)
        PSI_2 = PSI ** 2
        PSI_2 = torch.clamp(PSI_2, min=MIN, max=MAX)
        PSI_3 = PSI ** 3
        PSI_3 = torch.clamp(PSI_3, min=MIN, max=MAX)
        rborn_offset_batch_inv = 1.0 / rborn_offset_batch
        rborn_batch_inv = 1.0 / rborn_batch

        r_eff_inv = rborn_offset_batch_inv - torch.tanh(
            gb_alpha_batch * PSI - gb_beta_batch * PSI_2 + gb_gamma_batch * PSI_3) * rborn_batch_inv
        r_eff_inv = torch.clamp(r_eff_inv, min=MIN, max=MAX)
        r_eff[batch_start:batch_end] = 1 / r_eff_inv

    # _________________________________________________________________________________
    # Ending the loop to calculate the effective radius, now calculate solvation energy
    # re-calculate pairwise distances and use equation 4.2 and 4.3 from AMBER manual
    # to find electrostatic solvation energy
    # _________________________________________________________________________________
    batch_size = 1024  # You can adjust this based on your GPU memory
    solvation_energy = torch.zeros(n_frames, device='cuda', dtype=torch.float32)

    for batch_start in range(0, n_frames, batch_size):
        batch_end = min(batch_start + batch_size, n_frames)
        batch_size_current = batch_end - batch_start

        coords_batch = coords[batch_start:batch_end]
        r_eff_batch = r_eff[batch_start:batch_end]

        coords_i = coords_batch.unsqueeze(2)
        coords_j = coords_batch.unsqueeze(1)

        diff = coords_i - coords_j
        rij_squared = torch.sum(diff ** 2, dim=-1)

        Ri = r_eff_batch.unsqueeze(2).expand(-1, -1, n_atoms)
        Ri = torch.clamp(Ri, min=MIN, max=MAX)
        Rj = r_eff_batch.unsqueeze(1).expand(-1, n_atoms, -1)
        Rj = torch.clamp(Rj, min=MIN, max=MAX)
        RiRj = Ri * Rj
        fGB = torch.sqrt(rij_squared + RiRj * torch.exp(-rij_squared / (4 * RiRj)))
        fGB = torch.clamp(fGB, min=MIN, max=MAX)
        qiqj = (charges.unsqueeze(0).unsqueeze(2) * charges.unsqueeze(0).unsqueeze(1))
        exp_term = torch.exp(-kappa * fGB)
        solvation_energy_batch = -0.5 * torch.sum(qiqj / fGB * (1 - (exp_term / extdiel)), dim=(1, 2))
        solvation_energy[batch_start:batch_end] = solvation_energy_batch

    return solvation_energy


if __name__ == "__main__":
    dcd_file = "../traj/cyclicpeptide/TVGGVG/TVGGVG_all_rmrp_min.dcd"
#    dcd_file2 = '../CyclicPeptideOutput/spherical_expansion.dcd'
#    dcd_file3 = '../CyclicPeptideOutput/Driver_Cyclic_Peptide_Prediction.dcd'
#    dcd_file4 = '../../Pentane/PentaneCode/pentane_r1_10ps_100ns_nowat.dcd'
#    dcd_file5 = '../../Cyclohexane/CyclohexaneCode/CHX_r1_100ns_nowat.dcd'

    prmtop_file = "../traj/cyclicpeptide/TVGGVG/protein.prmtop"
#    prmtop_file2 = '../../Pentane/PentaneCode/pentane.prmtop'
#    prmtop_file3 = '../../Cyclohexane/CyclohexaneCode/CHX.prmtop'

    energy = calculate_energy_test(dcd_file, prmtop_file)
    print("\nMean Total Energy:", energy)
