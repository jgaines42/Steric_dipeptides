
def rotate_end_CH3(Position, delta_term_CH3_1, CH3_1_index, moveAtomID_CH3_1,  clash_list, radii_2, min_E, CH3_clash_list, CH3_radii_list):
    import numpy as np
    from rotate_DA import rotate_DA
    Pos_b4_CH3 = Position.copy()

    #print(CH3_clash_list)
    #print(CH3_radii_list)
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
        #Position = Pos_b4_CH3.copy()
        setCH3_1 = c3_1_loop * 5.0
        Position = rotate_DA(Pos_b4_CH3.copy(), setCH3_1, delta_term_CH3_1, CH3_1_index, moveAtomID_CH3_1)
        all_CH3_1[c3_1_loop] = Position[moveAtomID_CH3_1, :]


    found_0 = 0

    # Target E is the energy if there were no CH3 clashes 
    target_E = min_E - CH3_E
    total_E = min_E
    min_Position = Pos_b4_CH3.copy()
    # now do combinations
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
            if (this_E < CH3_E):
                CH3_E = this_E
                diff_pos = Position[clash_list[:, 0], :] - Position[clash_list[:, 1], :]
                sum_2 = np.sum(np.square(diff_pos), 1)
                ind0 = sum_2 < radii_2
                s_r_6 = np.power(radii_2[ind0] / sum_2[ind0], 3)
                E = np.power(1 - s_r_6, 2)
                total_E = np.sum(E)
                min_Position = Position.copy()
            if (CH3_E == 0):
                diff_pos = Position[clash_list[:, 0], :] - Position[clash_list[:, 1], :]
                sum_2 = np.sum(np.square(diff_pos), 1)
                ind0 = sum_2 < radii_2
                s_r_6 = np.power(radii_2[ind0] / sum_2[ind0], 3)
                E = np.power(1 - s_r_6, 2)
                total_E = np.sum(E)
                found_0 = 1
                return total_E, Position
    return total_E, min_Position