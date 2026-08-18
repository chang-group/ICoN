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
    #root = [atoms[155], atoms[156], atoms[157], atoms[160], atoms[161]]
    #root





    def check_overlap(list1, list2):
        for elem1 in list1:
            for elem2 in list2:
                if elem1 == elem2:
                    return True
        return False

    # HIV Protease stuff

    # This is multifragment section. Only use when needed
    # new_pair = [1563, root[2].index]   # or use [1563, 1564]
    # new_pair = [1560, 1561]   # or use [1563, 1564]
    # Create a new list that includes all original pairs plus the new pair
    # bond_list_new = bond_list + [new_pair]

    # pair_to_add = [3154,root[2].index]
    # bond_list_f = bond_list_new
    # bond_list_f = bond_list_f + [pair_to_add]

    if multifragment:
        pair_to_remove = [86, 88]
        bond_list.remove(pair_to_remove)
        pair_to_remove = [154, 155]
        bond_list.remove(pair_to_remove)


        pair_to_remove = [139, 142]
        bond_list.remove(pair_to_remove)
        bond_list = bond_list + [[139, root[2].index]]
        pair_to_remove = [93, 96]
        bond_list.remove(pair_to_remove)
        bond_list = bond_list + [[96, root[2].index]]
        pair_to_remove = [4, 6]
        bond_list.remove(pair_to_remove)
        bond_list = bond_list + [[4, root[2].index]]

    #pair_to_remove2 = [154, 155]
    # # # Remove the pair
    #bond_list.remove(pair_to_remove)
    # bond_list.remove(pair_to_remove2)
    # #bond_list_f = bond_list_f + [pair_to_add]
    # #bond_list_f
    #
    # print("origonal bond list size", len(bond_list))
    # root_indicies = [root[0].index, root[1].index, root[2].index, root[3].index, root[4].index]
    # k = 0
    # for i in bond_list:
    #
    #     if k % 15 == 0 and not check_overlap(i, root_indicies):
    #         print("i in loop", i)
    #         print("k in loop", k)
    #         bond_list.remove([i[0], i[1]])
    #         print("new bond", [i[1], root[2].index])
    #         bond_list = bond_list + [[i[1], root[2].index]]
    #     k += 1
    #
    # print("new bond list", bond_list)
    # print("new bond list size", len(bond_list))
    #find torsions

    #root = [atoms[155], atoms[156], atoms[157], atoms[160], atoms[161]]
    #root = [atoms[177], atoms[178], atoms[179], atoms[182], atoms[184]]
    root = [atoms[162], atoms[163], atoms[164], atoms[165], atoms[172]]
    torsion_XYZ_inds = find_torsions(root,atoms,bond_list)
    torsion_XYZ_inds = np.array(torsion_XYZ_inds)
    root = np.array(root)
    return torsion_XYZ_inds, n_atoms, root
    #np.save(out_path + "torsion_XYZ_inds.npy",torsion_XYZ_inds)
