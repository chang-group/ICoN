import torch


class Atom:
    """
    Represents an atom with a name, mass, and index.

    Attributes:
        name (str): The name of the atom.
        mass (float): The mass of the atom.
        index (int): The index of the atom in the molecule.
    """

    def __init__(self, name, mass, index):
        self.name = name
        self.mass = mass
        self.index = index
        self.bond_force_constants = []
        self.bond_equil_values = []
        self.angle_force_constants = []
        self.angle_equil_values = []
        self.dihedral_force_constants = []
        self.dihedral_periodicities = []
        self.dihedral_phases = []
        self.vdw_a_coef = []
        self.vdw_b_coef = []
        self.parser = None

    def __repr__(self):
        return f"Atom(name={self.name}, mass={self.mass}, index={self.index})"

    def bonded_atoms(self, bond_list=None):
        """
        Return the list of atoms bonded to this atom.

        Parameters:
            bond_list (list of tuple): Optional. A list of bonds to use instead of the default parser bonds.

        Returns:
            list of Atom: The atoms bonded to this atom.
        """
        if bond_list is None:
            # Use the default bonds from the parser if no custom bond list is provided
            bond_list = self.parser.bonds()
        adjacency_list = create_adjacency_list(bond_list)
        return [self.parser.atoms[i] for i in adjacency_list.get(self.index, [])]


class PRMTOPParser:
    """
    Parses a PRMTOP file to extract atom and bond information.

    Attributes:
        filename (str): The path to the PRMTOP file.
        atoms (list of Atom): List of Atom objects.
        atom_names (list of str): List of atom names.
        masses (list of float): List of atom masses.
        indices (list of int): List of atom indices.
        bond_indices (list of int): List of bond indices.
    """

    def __init__(self, filename):
        self.filename = filename
        self.pointers = []
        self.atoms = []
        self.atom_names = []
        self.residue_label = []
        self.atomic_number = []
        self.charges = []
        self.masses = []
        self.indices = []
        self.bond_indices = []
        self.angles_inc_hydrogen = []
        self.angles_without_hydrogen = []
        self.dihedrals_inc_hydrogen = []
        self.dihedrals_without_hydrogen = []
        self.residue_pointer = []
        self.bond_force_constants = []
        self.bond_equil_values = []
        self.angle_force_constants = []
        self.angle_equil_values = []
        self.dihedral_force_constants = []
        self.dihedral_periodicities = []
        self.dihedral_phases = []
        self.lennard_jones_acoef = []
        self.lennard_jones_bcoef = []
        self.atom_type_index = []
        self.nonbonded_parm_index = []
        self.number_excluded_atoms = []
        self.excluded_atoms_list = []
        self.scee_scale_factor = []
        self.scnb_scale_factor = []
        self.rborn = []
        self.fs = []
        self.sections = {
            'POINTERS': False,
            'ATOM_NAME': False,
            'RESIDUE_LABEL': False,
            'CHARGE': False,
            'MASS': False,
            'ATOMIC_NUMBER': False,
            'BONDS_INC_HYDROGEN': False,
            'BONDS_WITHOUT_HYDROGEN': False,
            'ANGLES_INC_HYDROGEN': False,
            'ANGLES_WITHOUT_HYDROGEN': False,
            'DIHEDRALS_WITHOUT_HYDROGEN': False,
            'DIHEDRALS_INC_HYDROGEN': False,
            'RESIDUE_POINTER': False,
            'BOND_FORCE_CONSTANT': False,
            'BOND_EQUIL_VALUE': False,
            'ANGLE_FORCE_CONSTANT': False,
            'ANGLE_EQUIL_VALUE': False,
            'DIHEDRAL_FORCE_CONSTANT': False,
            'DIHEDRAL_PERIODICITY': False,
            'DIHEDRAL_PHASE': False,
            'LENNARD_JONES_ACOEF' : False,
            'LENNARD_JONES_BCOEF' : False,
            'ATOM_TYPE_INDEX': False,
            'NONBONDED_PARM_INDEX': False,
            'NUMBER_EXCLUDED_ATOMS': False,
            'EXCLUDED_ATOMS_LIST': False,
            'SCEE_SCALE_FACTOR': False,
            'SCNB_SCALE_FACTOR': False,
            'RADII': False,
            'SCREEN': False

        }

        self.parse_prmtop()
        self.initialize_atoms()

    def get_prmtop_params(self):

        return {
            'pointers': torch.tensor(self.pointers),
            'atom_names': self.atom_names,
            'residue_label': self.residue_label,
            'bond_indices': self.get_bond_indices(),
            'atomic_number': torch.tensor(self.atomic_number),
            'residue_pointer': torch.tensor(self.residue_pointer),
            'bond_k': torch.tensor(self.bond_force_constants),
            'bond_r0': torch.tensor(self.bond_equil_values),
            'angle_indices': self.get_angle_indices(),
            'angle_k': torch.tensor(self.angle_force_constants),
            'angle_theta0': torch.tensor(self.angle_equil_values),
            'torsion_v': torch.tensor(self.dihedral_force_constants),
            'torsion_n': torch.tensor(self.dihedral_periodicities),
            'torsion_gamma': torch.tensor(self.dihedral_phases),
            'torsion_indices': self.get_torsion_indices(),
            'torsion_indices_with_sign': self.get_torsion_indices_with_sign(),
            'atom_type_index': torch.tensor(self.atom_type_index),
            'nonbonded_parm_index': torch.tensor(self.nonbonded_parm_index),
            'number_excluded_atoms': torch.tensor(self.number_excluded_atoms),
            'excluded_atoms_list': torch.tensor(self.excluded_atoms_list),
            'lennard_jones_acoef': torch.tensor(self.lennard_jones_acoef),
            'lennard_jones_bcoef': torch.tensor(self.lennard_jones_bcoef),
            'scee_scale_factor': torch.tensor(self.scee_scale_factor),
            'scnb_scale_factor': torch.tensor(self.scnb_scale_factor),
            'charges': torch.tensor(self.charges),
            'rborn': torch.tensor(self.rborn),
            'fs': torch.tensor(self.fs),
        }

    def get_bond_indices(self):
        return [[self.bond_indices[i] // 3, self.bond_indices[i + 1] // 3, self.bond_indices[i + 2]]
                for i in range(0, len(self.bond_indices), 3)]

    def get_angle_indices(self):
        angles = self.angles_inc_hydrogen + self.angles_without_hydrogen
        return [[angles[i] // 3, angles[i + 1] // 3, angles[i + 2] // 3, angles[i + 3]] for i in range(0, len(angles), 4)]

    def get_torsion_indices(self):
        dihedrals = self.dihedrals_inc_hydrogen + self.dihedrals_without_hydrogen
        return [
            [abs(dihedrals[i]) // 3, abs(dihedrals[i + 1]) // 3, abs(dihedrals[i + 2]) // 3, abs(dihedrals[i + 3]) // 3, abs(dihedrals[i + 4])]
            for i in range(0, len(dihedrals), 5)]

    def get_torsion_indices_with_sign(self):
        dihedrals_with_sign = self.dihedrals_inc_hydrogen + self.dihedrals_without_hydrogen
        return [
            [(dihedrals_with_sign[i]) // 3, (dihedrals_with_sign[i + 1]) // 3, (dihedrals_with_sign[i + 2]) // 3, (dihedrals_with_sign[i + 3]) // 3,
             (dihedrals_with_sign[i + 4])]
            for i in range(0, len(dihedrals_with_sign), 5)]

    def reset_sections(self):
        """Reset section flags to False."""
        for key in self.sections:
            self.sections[key] = False

    def parse_prmtop(self):
        """Parse the PRMTOP file and extract relevant data."""
        with open(self.filename, 'r') as file:
            lines = file.readlines()

        atom_name_width = 4

        for line in lines:
            if line.startswith('%FLAG'):
                self.reset_sections()
                section_name = line.split()[1]
                if section_name in self.sections:
                    self.sections[section_name] = True
                continue

            if line.startswith('%FORMAT'):
                continue

            elif self.sections['POINTERS']:
                self.pointers.extend(map(float, line.split()))
            if self.sections['ATOM_NAME']:
                names = [line[i:i + 4].strip() for i in range(0, len(line), 4)]
                self.atom_names.extend([name for name in names if name])
            if self.sections['RESIDUE_LABEL']:
                labels = [line[i:i + 4].strip() for i in range(0, len(line), 4)]
                self.residue_label.extend([labels for labels in labels if labels])
            elif self.sections['CHARGE']:
                self.charges.extend(map(float, line.split()))
            elif self.sections['MASS']:
                self.masses.extend(map(float, line.split()))
            elif self.sections['ATOMIC_NUMBER']:
                self.atomic_number.extend(map(float, line.split()))
            elif self.sections['BONDS_INC_HYDROGEN'] or self.sections['BONDS_WITHOUT_HYDROGEN']:
                self.bond_indices.extend(map(int, line.split()))
            elif self.sections['RESIDUE_POINTER']:
                self.residue_pointer.extend(map(float, line.split()))
            elif self.sections['BOND_FORCE_CONSTANT']:
                self.bond_force_constants.extend(map(float, line.split()))
            elif self.sections['BOND_EQUIL_VALUE']:
                self.bond_equil_values.extend(map(float, line.split()))
            elif self.sections['ANGLE_FORCE_CONSTANT']:
                self.angle_force_constants.extend(map(float, line.split()))
            elif self.sections['ANGLE_EQUIL_VALUE']:
                self.angle_equil_values.extend(map(float, line.split()))
            elif self.sections['DIHEDRAL_FORCE_CONSTANT']:
                self.dihedral_force_constants.extend(map(float, line.split()))
            elif self.sections['DIHEDRAL_PERIODICITY']:
                self.dihedral_periodicities.extend(map(float, line.split()))
            elif self.sections['DIHEDRAL_PHASE']:
                self.dihedral_phases.extend(map(float, line.split()))
            elif self.sections['LENNARD_JONES_ACOEF']:
                self.lennard_jones_acoef.extend(map(float, line.split()))
            elif self.sections['LENNARD_JONES_BCOEF']:
                self.lennard_jones_bcoef.extend(map(float, line.split()))
            elif self.sections['ANGLES_INC_HYDROGEN']:
                self.angles_inc_hydrogen.extend(map(int, line.split()))
            elif self.sections['ANGLES_WITHOUT_HYDROGEN']:
                self.angles_without_hydrogen.extend(map(int, line.split()))
            elif self.sections['DIHEDRALS_INC_HYDROGEN']:
                self.dihedrals_inc_hydrogen.extend(map(int, line.split()))
            elif self.sections['DIHEDRALS_WITHOUT_HYDROGEN']:
                self.dihedrals_without_hydrogen.extend(map(int, line.split()))
            elif self.sections['ATOM_TYPE_INDEX']:
                self.atom_type_index.extend(map(int, line.split()))
            elif self.sections['NONBONDED_PARM_INDEX']:
                self.nonbonded_parm_index.extend(map(int, line.split()))
            elif self.sections['NUMBER_EXCLUDED_ATOMS']:
                self.number_excluded_atoms.extend(map(int, line.split()))
            elif self.sections['EXCLUDED_ATOMS_LIST']:
                self.excluded_atoms_list.extend(map(int, line.split())),
            elif self.sections['SCNB_SCALE_FACTOR']:
                self.scnb_scale_factor.extend(map(float, line.split()))
            elif self.sections['SCEE_SCALE_FACTOR']:
                self.scee_scale_factor.extend(map(float, line.split()))
            elif self.sections['SCNB_SCALE_FACTOR']:
                self.scnb_scale_factor.extend(map(float, line.split()))
            elif self.sections['RADII']:
                self.rborn.extend(map(float, line.split()))
            elif self.sections['SCREEN']:
                self.fs.extend(map(float, line.split()))

    def initialize_atoms(self):
        """Initialize Atom objects from parsed data."""
        num_atoms = len(self.atom_names)
        if num_atoms != len(self.masses):
            raise ValueError(
                f"Number of atom names ({num_atoms}) does not match number of masses ({len(self.masses)}).")

        for i in range(num_atoms):
            atom = Atom(self.atom_names[i], self.masses[i], i)
            atom.parser = self
            self.indices.append(i)
            self.atoms.append(atom)

    def get_atoms(self):
        """Return the list of Atom objects."""
        return self.atoms

    def bonds(self):
        """Return a list of bonds as pairs of atom indices."""
        return [[self.bond_indices[i] // 3, self.bond_indices[i + 1] // 3] for i in range(0, len(self.bond_indices), 3)]


def create_adjacency_list(bond_list):
    """Create an adjacency list representing bonded atoms."""
    # bond_list = self.bonds()
    adjacency_list = {}
    for a, b in bond_list:
        if a not in adjacency_list:
            adjacency_list[a] = []
        if b not in adjacency_list:
            adjacency_list[b] = []
        adjacency_list[a].append(b)
        adjacency_list[b].append(a)
    return adjacency_list


def sort_atoms_by_mass_a2(atoms, reverse=False):
    """
    Sort Atom in the direction from C to N, and thus find C, CA first, then other carbon, then N,O backbone, then by mass.
    This ensures the side chain atom sort correctly

    Parameters
    ----------
    atoms : list of Atom
        List to sort
    reverse : bool
        Atoms will be in descending order if True

    Returns
    -------
    list of Atom
        Sorted list
    """

    def sorting_key(atom):
        is_priority_carbon = atom.name in ['CA']
        is_secondary_carbon = atom.name in ['C']
        is_other_carbon = 11.5 < atom.mass < 12.5 and not is_priority_carbon
        is_n_or_o = atom.name in ['N', 'O']
        # Priority order: carbon atoms with names 'C' or 'CA', other carbon atoms, then N and O, then by mass and index
        return (
        not is_priority_carbon, not is_secondary_carbon, not is_n_or_o, not is_other_carbon, -atom.mass, -atom.index)

    sorted_atoms = sorted(atoms, key=sorting_key)

    if reverse:
        sorted_atoms.reverse()

    return sorted_atoms


def find_shortest_path_to_root(specific_atom, root_atoms, atom_bond_list):
    # Build the graph as an adjacency list
    graph = {}
    for bond in atom_bond_list:
        a, b = bond
        if a not in graph:
            graph[a] = []
        if b not in graph:
            graph[b] = []
        graph[a].append(b)
        graph[b].append(a)

    # Initialize the queue and visited set
    queue = [(specific_atom, 0)]
    visited = {specific_atom}
    distances = {specific_atom: 0}

    # BFS traversal
    while queue:
        current_atom, current_distance = queue.pop(0)

        # Explore neighbors
        for neighbor in graph[current_atom]:
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, current_distance + 1))
                distances[neighbor] = current_distance + 1

    # Find the shortest distance to any root atom
    shortest_distance = float('inf')
    for root in root_atoms:
        if root in distances and distances[root] < shortest_distance:
            shortest_distance = distances[root]

    return shortest_distance


def find_torsions(root, atoms, bond_list_new):
    """Constructs a list of torsion angles

    Parameters
    ----------
    root : AtomGroup
        First three atoms in the coordinate system
    atoms : AtomGroup
        Atoms that are allowed to be part of the torsion angle

    Returns
    -------
    torsions : list of AtomGroup
        list of AtomGroup objects that define torsion angles
    """
    torsions = []
    selected_atoms = list(root)
    root_indices = [a.index for a in root]

    while len(selected_atoms) < len(atoms):
        torsionAdded = False

        for a1 in selected_atoms:
            # Find a0, which is a new atom connected to the selected atom
            a0_list = []
            for a in a1.bonded_atoms(bond_list_new):
                if (a in atoms) and (a not in selected_atoms):
                    a0_list.append(a)
            a0_list = sort_atoms_by_mass_a2(a0_list)

            for a0 in a0_list:
                # Find a2 list, which is connected to a1, is not a terminal atom, and has been selected
                a2_list = [a for a in a1.bonded_atoms(bond_list_new)
                           if a != a0 and len(a.bonded_atoms(bond_list_new)) > 1 and a in atoms and a in selected_atoms]
                if len(a2_list) > 0:
                    a2_list = sorted(a2_list,
                                     key=lambda x: find_shortest_path_to_root(x.index, root_indices, bond_list_new))

                for a2 in a2_list:
                    # Find a3 list, which is connected to a2, has been selected, and is not a1
                    a3_list = [a for a in a2.bonded_atoms(bond_list_new)
                               if a != a1 and a in atoms and a in selected_atoms]
                    if len(a3_list) > 0:
                        a3_list = sorted(a3_list,
                                         key=lambda x: find_shortest_path_to_root(x.index, root_indices, bond_list_new))

                        a3 = a3_list[0]

                        torsions.append([a0.index, a1.index, a2.index, a3.index])
                        selected_atoms.append(a0)
                        torsionAdded = True
                        break  # Exit a3 loop

                    if torsionAdded:
                        break

        if not torsionAdded:
            print('Selected atoms:')
            print([a.index + 1 for a in selected_atoms])
            print('Torsions found:')
            print(torsions)
            print("selected atoms size", len(selected_atoms))
            raise ValueError('Additional torsions not found.')

    print(torsions)
    return torsions


def find_torsion_indices_with_atoms(torsion_lists, atom_indices):
    matching_indices = []
    for i, torsion in enumerate(torsion_lists):
        if any(atom in torsion for atom in atom_indices):
            matching_indices.append(i)
    return matching_indices




