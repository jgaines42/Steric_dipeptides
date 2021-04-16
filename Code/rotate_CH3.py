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
import numpy as np
# from rotate_DA import rotate_DA
from math import sin, cos

#@profile
def rotate_CH3(Position, delta_term_CH3_1, delta_term_CH3_2, CH3_1_index, CH3_2_index, moveAtomID_CH3_1, moveAtomID_CH3_2, clash_list, radii_2, min_E, CH3_clash_list, CH3_radii_list):

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
        #Position = rotate_DA(Pos_b4_CH3, setCH3_1, delta_term_CH3_1, CH3_1_index, moveAtomID_CH3_1)
        deltaChi1_F153 = delta_term_CH3_1 - setCH3_1 * np.pi / 180.0

        # Get the coordinates of the 2nd atom in the dihedral
        subtract_atom = Pos_b4_CH3[CH3_1_index[1], :]

        # Move all atoms so 2nd atom of dihedral is at orgin
        TempPosition = Pos_b4_CH3 - subtract_atom

        # Get the vector from the origin to the 3rd atom of the dihedral.
        # This is the vector that we will rotate around
        CAtoCB_F153 = -TempPosition[CH3_1_index[2], :]
        CAtoCB_F153 = CAtoCB_F153 / np.linalg.norm(CAtoCB_F153)

        # Do complicated math
        q0 = cos(deltaChi1_F153 / 2.0)
        sindelta = sin(deltaChi1_F153 / 2.0)
        q1 = CAtoCB_F153[0] * sindelta
        q2 = CAtoCB_F153[1] * sindelta
        q3 = CAtoCB_F153[2] * sindelta
        q02 = q0 * q0
        q12 = q1 * q1
        q22 = q2 * q2
        q32 = q3 * q3

        # Q is the rotation vector
        Q = np.array([[(q02 + q12 - q22 - q32), 2.0 * (q1 * q2 - q0 * q3), 2.0 * (q0 * q2 + q1 * q3)],
                     [2.0 * (q1 * q2 + q0 * q3), (q02 - q12 + q22 - q32), 2.0 * (-q0 * q1 + q2 * q3)],
                     [2.0 * (-q0 * q2 + q1 * q3), 2 * (q0 * q1 + q2 * q3), (q02 - q12 - q22 + q32)]])

        # Extract the coordinates that should move
        move_coordinates = TempPosition[moveAtomID_CH3_1, :]

        # Use Q to rotate these atoms
        newPos = np.matmul(Q, move_coordinates.T)

        # Put the moved atoms back into the original array
        TempPosition[moveAtomID_CH3_1, :] = newPos.T

        # Translate atoms back to original location
        Position = TempPosition + subtract_atom
        all_CH3_1[c3_1_loop] = Position[moveAtomID_CH3_1, :]

    # Make all CH3 2 positions
    all_CH3_2 = [0] * 72
    for c3_2_loop in range(0, 72):
        setCH3_2 = c3_2_loop * 5.0
        #Position = rotate_DA(Pos_b4_CH3, setCH3_2, delta_term_CH3_2, CH3_2_index, moveAtomID_CH3_2)

        # Calculate how much the dihedral needs to change (in radians)
        deltaChi1_F153 = delta_term_CH3_2 - setCH3_2 * np.pi / 180.0

        # Get the coordinates of the 2nd atom in the dihedral
        subtract_atom = Pos_b4_CH3[CH3_2_index[1], :]

        # Move all atoms so 2nd atom of dihedral is at orgin
        TempPosition = Pos_b4_CH3 - subtract_atom

        # Get the vector from the origin to the 3rd atom of the dihedral.
        # This is the vector that we will rotate around
        CAtoCB_F153 = -TempPosition[CH3_2_index[2], :]
        CAtoCB_F153 = CAtoCB_F153 / np.linalg.norm(CAtoCB_F153)

        # Do complicated math
        q0 = cos(deltaChi1_F153 / 2.0)
        sindelta = sin(deltaChi1_F153 / 2.0)
        q1 = CAtoCB_F153[0] * sindelta
        q2 = CAtoCB_F153[1] * sindelta
        q3 = CAtoCB_F153[2] * sindelta
        q02 = q0 * q0
        q12 = q1 * q1
        q22 = q2 * q2
        q32 = q3 * q3

        # Q is the rotation vector
        Q = np.array([[(q02 + q12 - q22 - q32), 2.0 * (q1 * q2 - q0 * q3), 2.0 * (q0 * q2 + q1 * q3)],
                     [2.0 * (q1 * q2 + q0 * q3), (q02 - q12 + q22 - q32), 2.0 * (-q0 * q1 + q2 * q3)],
                     [2.0 * (-q0 * q2 + q1 * q3), 2 * (q0 * q1 + q2 * q3), (q02 - q12 - q22 + q32)]])

        # Extract the coordinates that should move
        move_coordinates = TempPosition[moveAtomID_CH3_2, :]

        # Use Q to rotate these atoms
        newPos = np.matmul(Q, move_coordinates.T)

        # Put the moved atoms back into the original array
        TempPosition[moveAtomID_CH3_2, :] = newPos.T

        # Translate atoms back to original location
        Position = TempPosition + subtract_atom

        all_CH3_2[c3_2_loop] = Position[moveAtomID_CH3_2, :]

    found_0 = 0

    total_E = min_E
    min_Position = Pos_b4_CH3.copy()
    diff_pos = np.zeros([CH3_clash_list.shape[0], 3])
    
    # now do combinations
    for c3_1_loop in range(0, 72):
        Position[moveAtomID_CH3_1, :] = all_CH3_1[c3_1_loop]

        for c3_2_loop in range(0, 72):

            if (found_0 == 0):
                Position[moveAtomID_CH3_2, :] = all_CH3_2[c3_2_loop]
                diff_pos = Position[CH3_clash_list[:, 0], :] - Position[CH3_clash_list[:, 1], :]
                #print(diff_pos1.shape)
                # Check energy due to overlaps of H3 atoms
                #diff_pos[:, 0] = Position[CH3_clash_list[:, 0], 0] - Position[CH3_clash_list[:, 1], 0]
                #diff_pos[:, 1] = Position[CH3_clash_list[:, 0], 1] - Position[CH3_clash_list[:, 1], 1]
                #diff_pos[:, 2] = Position[CH3_clash_list[:, 0], 2] - Position[CH3_clash_list[:, 1], 2]
                #
                #sum_2 = np.sum(np.square(diff_pos), 1)
                sum_2 = diff_pos[:, 0] * diff_pos[:, 0] + diff_pos[:, 1] * diff_pos[:, 1] + diff_pos[:, 2] * diff_pos[:, 2]
                #sum_2 = np.sum(diff_pos * diff_pos, 1)
                ind0 = sum_2 < CH3_radii_list

                s_over_r = CH3_radii_list[ind0] / sum_2[ind0]
                s_r_6 = s_over_r * s_over_r * s_over_r
                E = (1 - s_r_6) * (1 - s_r_6)
                #s_r_6 = np.power(CH3_radii_list[ind0] / sum_2[ind0], 3)
                #E = np.power(1 - s_r_6, 2)
                this_E = sum(E)

                # If this energy is less than the lowest CH3 energy found so far
                if (this_E < CH3_E):

                    # Set new lowest CH3 energy
                    CH3_E = this_E

                    # Get the total energy of the system and store in total_E
                    diff_pos2 = Position[clash_list[:, 0], :] - Position[clash_list[:, 1], :]
                    #sum_2 = np.sum(np.square(diff_pos), 1)
                    sum_2 = np.sum(diff_pos2 * diff_pos2, 1)
                    ind0 = sum_2 < radii_2
                    s_over_r = radii_2[ind0] / sum_2[ind0]
                    s_r_6 = s_over_r * s_over_r * s_over_r
                    E = (1 - s_r_6) * (1 - s_r_6)
                    # s_r_6 = np.power(radii_2[ind0] / sum_2[ind0], 3)
                    # E = np.power(1 - s_r_6, 2)
                    total_E = sum(E)

                    # Save this position
                    min_Position = Position.copy()

                # If there is no energy due to overlaps of the CH3 hydrogens, end the loop
                if (CH3_E == 0):

                    # Get the total energy of the system and store in total E
                    diff_pos2 = Position[clash_list[:, 0], :] - Position[clash_list[:, 1], :]
                    #sum_2 = np.sum(np.square(diff_pos), 1)
                    sum_2 = np.sum(diff_pos2 * diff_pos2, 1)
                    ind0 = sum_2 < radii_2
                    # s_r_6 = np.power(radii_2[ind0] / sum_2[ind0], 3)
                    # E = np.power(1 - s_r_6, 2)
                    s_over_r = radii_2[ind0] / sum_2[ind0]
                    s_r_6 = s_over_r * s_over_r * s_over_r
                    E = (1 - s_r_6) * (1 - s_r_6)
                    total_E = sum(E)
                    found_0 = 1

                    return total_E, Position

    return total_E, min_Position
