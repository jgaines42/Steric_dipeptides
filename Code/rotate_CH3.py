##################################################################################
# [energy, positioon] = rotate_CH3(Position, delta_term_CH3_1, delta_term_CH3_2, CH3_1_index, CH3_2_index,
#                           moveAtomID_CH3_1, moveAtomID_CH3_2, clash_list, radii_2, min_E, CH3_clash_list, CH3_radii_list)
#
# Input:
#   Position:           coordinates of the atoms
#   delta_term_CH3_1    term used for rotation of first CH3 group
#   delta_term_CH3_2    term used for rotation of the second CH3 group
#   CH3_1_index         index of the atoms of the dihedral of the first CH3 group
#   CH3_2_index         index of the atoms of the dihedral of the second CH3 group
#   moveAtomID_CH3_1    atoms that move when first group is rotated
#   moveAtomID_CH3_2    atoms that move with the second group is rotated
#   clash_list          all possible clashes to check
#   radii_2             (radii + radii)^2 for all pairs of atoms
#   min_E               current energy of the system (we want to beat this)
#   CH3_clash_list      all possible clashes involving the CH3 hydrogens
#   CH3_radii_list      (radii + radii)^2 for all pairs of atoms involving the CH3 hydrogens
#
# Returns:
#   total_E: the new total energy of the system
#   min_Positions: the coordinates of the system at this minimum energy
##################################################################################
def rotate_CH3(Position, delta_term_CH3_1, delta_term_CH3_2, CH3_1_index, CH3_2_index, moveAtomID_CH3_1, moveAtomID_CH3_2, clash_list, radii_2, min_E, CH3_clash_list, CH3_radii_list):
    import numpy as np
    from rotate_DA import rotate_DA

    Pos_b4_CH3 = Position.copy()

    # Get initial energy due to CH3 groups
    diff_pos = Position[CH3_clash_list[:, 0], :] - Position[CH3_clash_list[:, 1], :]
    sum_2 = np.sum(np.square(diff_pos), 1)
    ind0 = sum_2 < CH3_radii_list
    s_r_6 = np.power(CH3_radii_list[ind0] / sum_2[ind0], 3)
    E = np.power(1 - s_r_6, 2)
    CH3_E = np.sum(E)

    # Make all CH3 1 positions
    all_CH3_1 = [0] * 72
    for c3_1_loop in range(0, 72):
        setCH3_1 = c3_1_loop * 5.0
        Position = rotate_DA(Pos_b4_CH3, setCH3_1, delta_term_CH3_1, CH3_1_index, moveAtomID_CH3_1)
        all_CH3_1[c3_1_loop] = Position[moveAtomID_CH3_1, :]

    # Make all CH3 2 positions
    all_CH3_2 = [0] * 72
    for c3_2_loop in range(0, 72):
        setCH3_2 = c3_2_loop * 5.0
        Position = rotate_DA(Pos_b4_CH3, setCH3_2, delta_term_CH3_2, CH3_2_index, moveAtomID_CH3_2)
        all_CH3_2[c3_2_loop] = Position[moveAtomID_CH3_2, :]

    found_0 = 0

    total_E = min_E
    min_Position = Pos_b4_CH3.copy()

    # now do combinations
    for c3_1_loop in range(0, 72):
        Position[moveAtomID_CH3_1, :] = all_CH3_1[c3_1_loop]

        for c3_2_loop in range(0, 72):

            if (found_0 == 0):
                Position[moveAtomID_CH3_2, :] = all_CH3_2[c3_2_loop]

                # Check energy due to overlaps of H3 atoms
                diff_pos = Position[CH3_clash_list[:, 0], :] - Position[CH3_clash_list[:, 1], :]
                sum_2 = np.sum(np.square(diff_pos), 1)
                ind0 = sum_2 < CH3_radii_list
                s_r_6 = np.power(CH3_radii_list[ind0] / sum_2[ind0], 3)
                E = np.power(1 - s_r_6, 2)
                this_E = np.sum(E)

                # If this energy is less than the lowest CH3 energy found so far
                if (this_E < CH3_E):

                    # Set new lowest CH3 energy
                    CH3_E = this_E

                    # Get the total energy of the system and store in total_E
                    diff_pos = Position[clash_list[:, 0], :] - Position[clash_list[:, 1], :]
                    sum_2 = np.sum(np.square(diff_pos), 1)
                    ind0 = sum_2 < radii_2
                    s_r_6 = np.power(radii_2[ind0] / sum_2[ind0], 3)
                    E = np.power(1 - s_r_6, 2)
                    total_E = np.sum(E)

                    # Save this position
                    min_Position = Position.copy()

                # If there is no energy due to overlaps of the CH3 hydrogens, end the loop
                if (CH3_E == 0):

                    # Get the total energy of the system and store in total E
                    diff_pos = Position[clash_list[:, 0], :] - Position[clash_list[:, 1], :]
                    sum_2 = np.sum(np.square(diff_pos), 1)
                    ind0 = sum_2 < radii_2
                    s_r_6 = np.power(radii_2[ind0] / sum_2[ind0], 3)
                    E = np.power(1 - s_r_6, 2)
                    total_E = np.sum(E)
                    found_0 = 1

                    return total_E, Position

    return total_E, min_Position
