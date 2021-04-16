############################################################################
# run_dihedral.py
#
# run_dihedral.py folder coord_file AA_name coord_num save_folder save_name
#
# Loads the gro file of a Valine residue and rotates the full phi, psi, chi space
# Calculates the energy due to atomic overlaps at each position and writes to a file
# Angles are sampled every 5 degree
#
# Input:
#   folder: path to folder with the coordinate file
#   coord_file: file name containing the coordinates. Coordinates must be in gromacs order but in Angstroms
#   AA_name: amino acid 3 letter abbreviation
#   coord_num: integer of which set of corodinates this is
#   save_folder: folder to save data
#   save_name: additional string to use in save file names
#
# Required files
#   coord_file: coordinates in Gro order, columns = x, y, z, in Angstroms
#   Val_bonds.txt: atom indexes (1 indexed) of all atomic bonds
#   Val_angles.txt: atom indexes of all bond angles
#   Val_radii.txt: Atomic radii (in Angstroms) for each atom
#   Val_Hbonds.txt: 1 indexed list of potential hydrogen bonding atoms
#
# Saves:
#   save_folder + 'Val_energy_' + save_name + '_' + coord_num + '.txt'
#       Columns:
#           1: phi
#           2: psi
#           3: chi 1
#           4: Energy
#
############################################################################


import sys
from calcDihedral import calcDihedral
#from rotate_DA import rotate_DA
from create_clash_list import create_clash_list
from rotate_CH3 import rotate_CH3
from rotate_end_CH3 import rotate_end_CH3
import time
import numpy as np
from math import sin, cos


def rotate_DA(Position, setChi, delta_term, iChiArray, moveAtomID):

    #import numpy as np
    
    # Calculate how much the dihedral needs to change (in radians)
    deltaChi1_F153 = delta_term - setChi * np.pi / 180.0

    # Get the coordinates of the 2nd atom in the dihedral
    subtract_atom = Position[iChiArray[1], :]

    # Move all atoms so 2nd atom of dihedral is at orgin
    TempPosition = Position - subtract_atom

    # Get the vector from the origin to the 3rd atom of the dihedral.
    # This is the vector that we will rotate around
    CAtoCB_F153 = -TempPosition[iChiArray[2], :]
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
    move_coordinates = TempPosition[moveAtomID, :]

    # Use Q to rotate these atoms
    newPos = np.matmul(Q, move_coordinates.T)

    # Put the moved atoms back into the original array
    TempPosition[moveAtomID, :] = newPos.T

    # Translate atoms back to original location
    new_Pos1 = TempPosition + subtract_atom
    return new_Pos1


time1 = time.perf_counter()

# Get arguments
folder = sys.argv[1]
coord_file = sys.argv[2]
AA_name = sys.argv[3]
coord_num = sys.argv[4]
save_folder = sys.argv[5]
save_name = sys.argv[6]

# Load coordinates
Coord = np.loadtxt(folder + coord_file)        # in  Angstroms

# Load list of bonds, angles and atomic radii
bonds = np.loadtxt('Val_bonds.txt')                     # 1 indexed list of bond indexes
angles = np.loadtxt('Val_angles.txt')                   # 1 indexed list of bond angle indexes
radii = np.loadtxt('Val_radii.txt')                     # Atomic radii in Angstroms
hbonds = np.loadtxt('Val_Hbond.txt', dtype=int)

n_atoms = Coord.shape[0]

# fix bonds and angles to 0 indexing
bonds = bonds - 1
angles = angles - 1
hbonds = hbonds - 1

# Create clash list
clash_list = create_clash_list(n_atoms, bonds, angles)
radii_sum = radii[clash_list[:, 0]] + radii[clash_list[:, 1]]

# modify H-bond radii
for i in range(0, hbonds.shape[0]):
    ind1 = np.isin(clash_list[:, 0], hbonds[i, 0])
    ind2 = np.isin(clash_list[:, 1], hbonds[i, 1])
    radii_sum[ind1 & ind2] = 1.5

# Get sum of radii for all possible atom pairs in clash list
radii_2 = radii_sum * radii_sum                                 # Leave as (radii_sum)^2 to improve math later

# Get index of peptide caps
prior_res_index = [0, 1, 2, 3, 4, 5]
next_res_index = np.arange(n_atoms - 6, n_atoms)


# Initialize dihedral arrays
phi_index = np.array([4, 6, 8, n_atoms - 8])
psi_index = np.array([6, 8, n_atoms - 8, next_res_index[0]])
chi1_index = np.array([6, 8, 10, 12])
CH3_1_index = np.array([8, 10, 12, 13])
CH3_2_index = np.array([8, 10, 16, 17])
end_CH3_1_index = np.array([6, 4, 1, 0])
end_CH3_2_index = np.array([20, 22, 24, 25])

# Get initial dihedral angle values
phi_init = calcDihedral(phi_index, Coord)
psi_init = calcDihedral(psi_index, Coord)
chi1_init = calcDihedral(chi1_index, Coord)
CH3_1_init = calcDihedral(CH3_1_index, Coord)
CH3_2_init = calcDihedral(CH3_2_index, Coord)
end_CH3_1_init = calcDihedral(end_CH3_1_index, Coord)
end_CH3_2_init = calcDihedral(end_CH3_2_index, Coord)
omega_1 = calcDihedral([1, 4, 6, 8], Coord)
omega_2 = calcDihedral([8, 20, 22, 24], Coord)

# initial values of angles being rotated have to be > 0? not sure why
if (phi_init < 0):
    phi_init = phi_init + 360
if (psi_init < 0):
    psi_init = psi_init + 360
if (chi1_init < 0):
    chi1_init = chi1_init + 360

# Create list of atoms to move
moveAtomID_phi = np.arange(phi_index[2] + 1, n_atoms)
moveAtomID_psi = np.arange(psi_index[2] + 1, n_atoms)
moveAtomID_chi1 = np.arange(chi1_index[2] + 1, n_atoms - 8)
moveAtomID_CH3_1 = np.arange(13, 16)
moveAtomID_CH3_2 = np.arange(17, 20)
moveAtomID_end_CH3_1 = np.array([0, 2, 3])
moveAtomID_end_CH3_2 = np.array([25, 26, 27])

delta_term_phi = np.pi * np.sign(phi_init) * phi_init / 180.0
delta_term_psi = np.pi * np.sign(psi_init) * psi_init / 180.0
delta_term_chi1 = np.pi * np.sign(chi1_init) * chi1_init / 180.0
delta_term_CH3_1 = np.pi * np.sign(CH3_1_init) * CH3_1_init / 180.0
delta_term_CH3_2 = np.pi * np.sign(CH3_2_init) * CH3_2_init / 180.0
delta_term_end_CH3_1 = np.pi * np.sign(end_CH3_1_init) * end_CH3_1_init / 180.0
delta_term_end_CH3_2 = np.pi * np.sign(end_CH3_2_init) * end_CH3_2_init / 180.0

# get clash list that include side chain CH3 atoms
ind0 = np.isin(clash_list[:, 0], moveAtomID_CH3_1)
ind1 = np.isin(clash_list[:, 1], moveAtomID_CH3_1)
ind2 = np.isin(clash_list[:, 0], moveAtomID_CH3_2)
ind3 = np.isin(clash_list[:, 1], moveAtomID_CH3_2)
CH3_clash_list = clash_list[ind0 | ind1 | ind2 | ind3, :]
CH3_radii_list = radii_2[ind0 | ind1 | ind2 | ind3]


# get clash list that includes end H3 atoms
ind0 = np.isin(clash_list[:, 0], moveAtomID_end_CH3_1)
ind1 = np.isin(clash_list[:, 1], moveAtomID_end_CH3_1)
CH3_end1_clash_list = clash_list[ind0 | ind1, :]
CH3_end1_radii_list = radii_2[ind0 | ind1]

ind0 = np.isin(clash_list[:, 0], moveAtomID_end_CH3_2)
ind1 = np.isin(clash_list[:, 1], moveAtomID_end_CH3_2)
CH3_end2_clash_list = clash_list[ind0 | ind1, :]
CH3_end2_radii_list = radii_2[ind0 | ind1]

# Get initial E
# Check clashes
diff_pos = Coord[clash_list[:, 0], :] - Coord[clash_list[:, 1], :]      # distance vector
sum_2 = np.sum(np.square(diff_pos), 1)                                  # distance^2
ind0 = sum_2 < radii_2                                                  # Compare the distance to the (sum_radii)^2
s_r_6 = np.power(radii_2[ind0] / sum_2[ind0], 3)
E = np.power(1 - s_r_6, 2)
total_E = np.sum(E)

min_E = 1000000


# All energy data
all_energy = np.zeros((36 * 36 * 36, 4))
counter = 0

# Store initial position so you can keep moving back
Initial_Position = Coord

# Loop over all phi
for phi_loop in range(0, 1):
    setPhi = phi_loop * 10.0
    Pos_phi = rotate_DA(Initial_Position, setPhi, delta_term_phi, phi_index, moveAtomID_phi)
    print(setPhi)
    
    # After setting phi, CB won't move when we rotate psi, so we can pre-make all chi structures
    all_Chi_pos = [0] * 36
    for chi1_preloop in range(0, 36):
        setChi = chi1_preloop * 10.0
        Position = rotate_DA(Pos_phi, setChi, delta_term_chi1, chi1_index, moveAtomID_chi1)
        all_Chi_pos[chi1_preloop] = Position[moveAtomID_chi1, :]


    for psi_loop in range(0, 36):
        setPsi = psi_loop * 10.0
        Pos_phi_psi = rotate_DA(Pos_phi, setPsi, delta_term_psi, psi_index, moveAtomID_psi)
        print(setPsi)

        for chi1_loop in range(0, 36):
            Position_ppc = Pos_phi_psi.copy()
            Position_ppc[moveAtomID_chi1, :] = all_Chi_pos[chi1_loop]

            setChi = chi1_loop * 10.0
            total_E = 0

            # Now that we've rotated everything, check for clashes
            diff_pos = Position_ppc[clash_list[:, 0], :] - Position_ppc[clash_list[:, 1], :]
            #sum_2 = np.sum(np.square(diff_pos), 1)
            sum_2 = np.sum(diff_pos * diff_pos, 1)
            ind0 = sum_2 < radii_2

            s_over_r = radii_2[ind0] / sum_2[ind0]
            s_r_6 = s_over_r * s_over_r * s_over_r
            E = (1 - s_r_6) * (1 - s_r_6)
            #s_r_6 = np.power(radii_2[ind0] / sum_2[ind0], 3)
            #E = np.power(1 - s_r_6, 2)
            total_E = sum(E)

            # See if any of the clashes included the SC CH3 group
            clash_used = clash_list[ind0, :]
            ind1 = np.isin(moveAtomID_CH3_1, clash_used)
            ind2 = np.isin(moveAtomID_CH3_2, clash_used)

            # If so, rotate the SC CH3 groups
            if (total_E > 0 and (np.any(ind1) or np.any(ind2))):

                min_E = total_E
                # Rotate CH3 groups
                [min_E, Position_CH3] = rotate_CH3(Position_ppc.copy(), delta_term_CH3_1, delta_term_CH3_2, CH3_1_index,
                                                   CH3_2_index, moveAtomID_CH3_1, moveAtomID_CH3_2, clash_list, radii_2,
                                                   min_E, CH3_clash_list, CH3_radii_list)
                total_E = min_E
                Position_ppc = Position_CH3

            # See if any clashes include N-terminal end CH3

            ind1 = np.isin(moveAtomID_end_CH3_1, clash_used)
            if (total_E > 0 and np.any(ind1)):
                min_E = total_E
                [min_E, new_Pos] = rotate_end_CH3(Position_ppc.copy(), delta_term_end_CH3_1, end_CH3_1_index,
                                                  moveAtomID_end_CH3_1, clash_list, radii_2, min_E,
                                                  CH3_end1_clash_list, CH3_end1_radii_list)
                total_E = min_E
                Position_ppc = new_Pos
                print('rotate end')
            # See if any clashes include Cterminal end CH3
            ind2 = np.isin(moveAtomID_end_CH3_2, clash_used)
            
            if (total_E > 0 and np.any(ind2)):
                min_E = total_E
                [min_E, new_Pos] = rotate_end_CH3(Position_ppc.copy(), delta_term_end_CH3_2, end_CH3_2_index,
                                                  moveAtomID_end_CH3_2, clash_list, radii_2, min_E,
                                                  CH3_end2_clash_list, CH3_end2_radii_list)
                total_E = min_E
                Position_ppc = new_Pos
                print('rotate end')
            all_energy[counter, 0] = setPhi
            all_energy[counter, 1] = setPsi
            all_energy[counter, 2] = setChi
            all_energy[counter, 3] = total_E
            counter = counter + 1

        #time2 = time.perf_counter()
        #print('one loop', time2-time1)
np.savetxt(save_folder + 'Val_energy_' + save_name + '_' + coord_num + '.txt', all_energy, fmt='%d %d %d %7.4f')
time2 = time.perf_counter()
print('time: ', time2 - time1)