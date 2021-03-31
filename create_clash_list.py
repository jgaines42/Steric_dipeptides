def create_clash_list(n_atoms, bonds, angles, radii):
    import numpy as np

    for atom1 in range(0, n_atoms - 1):  # Loop over all first atoms

        # Get list of all atoms to check
        atom2 = np.arange(atom1 + 1, n_atoms)
        ind0 = bonds[:, 0] == atom1
        sub_bonds = bonds[ind0, :]

        # remove atoms that are bonded to this atom or in a bond angle
        if (sub_bonds.shape[0] > 0):
            ind1 = np.isin(atom2, sub_bonds[:, 1])
            atom2 = atom2[~ind1]

        ind0 = angles[:, 0] == atom1
        sub_angles = angles[ind0, :]

        if (sub_angles.shape[0] > 0):
            ind1 = np.isin(atom2, sub_angles[:, 2])
            atom2 = atom2[~ind1]

        if (atom1 == 0):
            clash_list = np.tile(atom1, [atom2.shape[0], 1])
            clash_list = np.append(clash_list, atom2)
            clash_list = clash_list.reshape((atom2.shape[0], 2), order='F')
            #print(clash_list)
            all_clash = clash_list.copy()
        else:
            clash_list = np.tile(atom1, [atom2.shape[0], 1])
            clash_list = np.append(clash_list, atom2)
            clash_list = clash_list.reshape((atom2.shape[0], 2), order='F')
            all_clash = np.append(all_clash, clash_list, axis=0)
    return all_clash