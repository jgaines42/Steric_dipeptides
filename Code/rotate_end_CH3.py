##################################################################################
# [energy, position] = rotate_end_CH3(Position, delta_term_CH3_1, CH3_1_index, moveAtomID_CH3_1,
#                                       clash_list, radii_2, min_E, CH3_clash_list, CH3_radii_list):
#
# Input:
#   Position:           coordinates of the atoms
#   delta_term_CH3_1    term used for rotation of the CH3 group
#   CH3_1_index         index of the atoms of the dihedral of the CH3 group
#   moveAtomID_CH3_1    atoms that move when the CH3 group is rotated
#   clash_list          all possible clashes to check
#   radii_2             (radii + radii)^2 for all pairs of atoms
#   min_E               current energy of the system (we want to beat this)
#   CH3_clash_list      all possible clashes involving the CH3 hydrogens
#   CH3_radii_list      (radii + radii)^2 for all pairs of atoms involving the CH3 hydrogens
#
# Returns:
#   total_E: the new total energy of the system
#   min_Positions: the coordinates of the system at this minimum energy
#
##################################################################################

def rotate_end_CH3(Position, delta_term_CH3_1, CH3_1_index, moveAtomID_CH3_1, clash_list, radii_2, min_E, CH3_clash_list, CH3_radii_list):
    import numpy as np
    from rotate_DA import rotate_DA

    Pos_b4_CH3 = Position.copy()

    # Get initial energy due to CH3 group
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
        Position = rotate_DA(Pos_b4_CH3.copy(), setCH3_1, delta_term_CH3_1, CH3_1_index, moveAtomID_CH3_1)
        all_CH3_1[c3_1_loop] = Position[moveAtomID_CH3_1, :]

    found_0 = 0

    total_E = min_E
    min_Position = Pos_b4_CH3.copy()

    # Loop over all possible CH3 positions
    for c3_1_loop in range(0, 72):

        if (found_0 == 0):
            Position[moveAtomID_CH3_1, :] = all_CH3_1[c3_1_loop]

            # Just check energy due to CH3 clashes
            diff_pos = Position[CH3_clash_list[:, 0], :] - Position[CH3_clash_list[:, 1], :]
            sum_2 = np.sum(np.square(diff_pos), 1)
            ind0 = sum_2 < CH3_radii_list
            s_r_6 = np.power(CH3_radii_list[ind0] / sum_2[ind0], 3)
            E = np.power(1 - s_r_6, 2)
            this_E = np.sum(E)

            # if this energy is less than the lowest CH3 energy so far
            if (this_E < CH3_E):
                # Store the new lowest energy
                CH3_E = this_E

                # Get full energy of the system
                diff_pos = Position[clash_list[:, 0], :] - Position[clash_list[:, 1], :]
                sum_2 = np.sum(np.square(diff_pos), 1)
                ind0 = sum_2 < radii_2
                s_r_6 = np.power(radii_2[ind0] / sum_2[ind0], 3)
                E = np.power(1 - s_r_6, 2)
                total_E = np.sum(E)

                # Store the new lowest energy positions
                min_Position = Position.copy()

            # If the energy due to CH3 was 0, store and exit
            if (CH3_E == 0):

                # Get full energy of the system
                diff_pos = Position[clash_list[:, 0], :] - Position[clash_list[:, 1], :]
                sum_2 = np.sum(np.square(diff_pos), 1)
                ind0 = sum_2 < radii_2
                s_r_6 = np.power(radii_2[ind0] / sum_2[ind0], 3)
                E = np.power(1 - s_r_6, 2)
                total_E = np.sum(E)
                found_0 = 1
                return total_E, Position

    return total_E, min_Position
