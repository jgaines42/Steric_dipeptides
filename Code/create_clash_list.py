##################################################################################
# all_clash = create_clash_list(n_atoms, bonds, angles)
#
# Determines all possible non-bonded clashes for set of atoms, bond, and angles
#
# Input:
#   n_atoms: The number of atoms in the system
#   bonds: N x 2 array of the index of bonded atoms
#   angles: N x 3 array of the index of atoms with bond angles between them
#
# Note: inputs should be 0 indexed and should be in numerical order
#       i.e a bond between atom 3 and 5 should be [3, 5] not [5, 3]
#
# Output:
#   all_clash: N x 2 array of indexes of atoms that could clash
#
##################################################################################

def create_clash_list(n_atoms, bonds, angles):
    import numpy as np

    # Loop over all first atoms
    for atom1 in range(0, n_atoms - 1):

        # Get list other atoms later in the list
        atom2 = np.arange(atom1 + 1, n_atoms)

        # Find all bonds that start with this atom
        ind0 = bonds[:, 0] == atom1
        sub_bonds = bonds[ind0, :]

        # remove atoms that are bonded to this atom
        if (sub_bonds.shape[0] > 0):
            ind1 = np.isin(atom2, sub_bonds[:, 1])
            atom2 = atom2[~ind1]

        # Find all angles that start with this atom
        ind0 = angles[:, 0] == atom1
        sub_angles = angles[ind0, :]

        # remove atoms that are in an angle with this atom
        if (sub_angles.shape[0] > 0):
            ind1 = np.isin(atom2, sub_angles[:, 2])
            atom2 = atom2[~ind1]

        # Store data in the clash list
        if (atom1 == 0):
            clash_list = np.tile(atom1, [atom2.shape[0], 1])
            clash_list = np.append(clash_list, atom2)
            clash_list = clash_list.reshape((atom2.shape[0], 2), order='F')
            all_clash = clash_list.copy()
        else:
            clash_list = np.tile(atom1, [atom2.shape[0], 1])
            clash_list = np.append(clash_list, atom2)
            clash_list = clash_list.reshape((atom2.shape[0], 2), order='F')
            all_clash = np.append(all_clash, clash_list, axis=0)
    return all_clash
