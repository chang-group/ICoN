import numpy as np

from Atom import *

def getTorsionList(path_in, top_path_in, multifragment=False):

    path = path_in
    top_path = top_path_in

    parser = PRMTOPParser(top_path)

    atoms = parser.get_atoms()
    n_atoms = len(atoms)
    n_torsions = n_atoms - 3
    bond_list = parser.bonds()
    adjacency_list = create_adjacency_list(bond_list)

    #Find Root

    initial_atom = None

    terminal_atoms = sort_atoms_by_mass_a2( \
        [a for a in atoms if len(adjacency_list[a.index]) == 1], reverse=True)
    if (initial_atom is None):
        # Select the heaviest root atoms from the heaviest terminal atom
        initial_atom = terminal_atoms[0]
    elif (not initial_atom in terminal_atoms):
        raise ValueError('Initial atom is not a terminal atom')

    second_atom = atoms[adjacency_list[initial_atom.index][0]]

    if atoms != 3:
        third_atom = sort_atoms_by_mass_a2( \
            [a for a in atoms \
            if (a.index in adjacency_list[second_atom.index]) and (a != initial_atom) and (a not in terminal_atoms)], \
            )[0]
    else:
        third_atom = sort_atoms_by_mass_a2( \
            [a for a in atoms \
            if (a.index in adjacency_list[second_atom.index]) and (a != initial_atom)], \
            )[0]

    forth_atoms = []
    for a in second_atom.bonded_atoms(bond_list):
        if a in terminal_atoms and a != initial_atom:
            forth_atoms.append(a)

    forth_atoms = sort_atoms_by_mass_a2(forth_atoms)

    root = [initial_atom, second_atom, third_atom] + forth_atoms
    #root = [atoms[0], atoms[2], atoms[3], atoms[12], atoms[13]] 
    print('Root atoms:', root)
    
    pair_to_remove = [0, 65]
    bond_list.remove(pair_to_remove)



    def check_overlap(list1, list2):
        for elem1 in list1:
            for elem2 in list2:
                if elem1 == elem2:
                    return True
        return False

    
    torsion_XYZ_inds = find_torsions(root,atoms,bond_list)
    torsion_XYZ_inds = np.array(torsion_XYZ_inds)
    root = np.array(root)
    return torsion_XYZ_inds, n_atoms, root
    #np.save(out_path + "torsion_XYZ_inds.npy",torsion_XYZ_inds)
