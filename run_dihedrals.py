############################################################################
# run_dihedral.py
#
# Loads the gro file of a Valine residue and rotates the full phi, psi, chi space
# Calculates the energy due to atomic overlaps at each position and writes to a file
# Angles are sampled every 5 degree
#
# Required files
#   Val_start_1.txt: coordinates in Gro order
#   Val_bonds.txt: atom indexes (1 indexed) of all atomic bonds
#   Val_angles.txt: atom indexes of all bond angles
#   Val_radii.txt: Atomic radii (in Angstroms) for each atom
#
# Calls:
#   create_clash_list
#   
############################################################################

import numpy as np
from calcDihedral import calcDihedral
from rotate_DA import rotate_DA
from create_clash_list import create_clash_list
from rotate_CH3 import rotate_CH3


# Load coordinates
coord_file = 'Val_start.txt'

Atoms = np.loadtxt(coord_file, usecols=[1], dtype='str')
Coord = np.loadtxt(coord_file, usecols=[3, 4, 5])        # in nm
Coord = Coord * 10.0                                     # Put in Angstroms

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
clash_list = create_clash_list(n_atoms, bonds, angles, radii)
np.savetxt('clash_list.txt', clash_list, fmt='%d')
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
end_CH3_1_index = np.array([0, 1, 4, 6])
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

#print(phi_init, psi_init, chi1_init, CH3_1_init, CH3_2_init)
#print(omega_1, omega_2)

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

delta_term_phi = np.pi * np.sign(phi_init) * phi_init / 180.0
delta_term_psi = np.pi * np.sign(psi_init) * psi_init / 180.0
delta_term_chi1 = np.pi * np.sign(chi1_init) * chi1_init / 180.0
delta_term_CH3_1 = np.pi * np.sign(CH3_1_init) * CH3_1_init / 180.0
delta_term_CH3_2 = np.pi * np.sign(CH3_2_init) * CH3_2_init / 180.0

# get clash list that include side chain CH3 atoms
ind0 = np.isin(clash_list[:, 0], moveAtomID_CH3_1)
ind1 = np.isin(clash_list[:, 1], moveAtomID_CH3_1)
ind2 = np.isin(clash_list[:, 0], moveAtomID_CH3_2)
ind3 = np.isin(clash_list[:, 1], moveAtomID_CH3_2)
CH3_clash_list = clash_list[ind0 | ind1 | ind2 | ind3, :]
CH3_radii_list = radii_2[ind0 | ind1 | ind2 | ind3]

# Get initial E
# Check clashes
diff_pos = Coord[clash_list[:, 0], :] - Coord[clash_list[:, 1], :]      # distance vector
sum_2 = np.sum(np.square(diff_pos), 1)                                  # distance^2
ind0 = sum_2 < radii_2                                                  # Compare the distance to the (sum_radii)^2
s_r_6 = np.power(radii_2[ind0] / sum_2[ind0], 3)
E = np.power(1 - s_r_6, 2)
total_E =  np.sum(E)
print('init E', total_E)
print(clash_list[ind0, :])


min_E = 1000000

# Move phi/psi dihedrals to 0
# Position = rotate_DA(Initial_Position, 0, delta_term_phi, phi_index, moveAtomID_phi)
# Position = rotate_DA(Position, 0, delta_term_psi, psi_index, moveAtomID_psi)
# phi_init = 0
# psi_init = 0
# delta_term_phi = np.pi * np.sign(phi_init) * phi_init / 180.0
# delta_term_psi = np.pi * np.sign(psi_init) * psi_init / 180.0
# Initial_Position = Position.copy()
# print(phi_init, psi_init)

# All energy data
all_energy = np.zeros((72 * 72 * 72, 4))
counter = 0

# Store initial position so you can keep moving back
Initial_Position = Coord
# Loop over all phi
for phi_loop in range(0, 72):
    #Position = Initial_Position.copy()
    setPhi = phi_loop * 5.0
    Pos_phi = rotate_DA(Initial_Position.copy(), setPhi, delta_term_phi, phi_index, moveAtomID_phi)
    print(setPhi)
    #Pos_phi = Position.copy()
    for psi_loop in range(0, 72):
        #Position = Pos_phi.copy()
        setPsi = psi_loop * 5.0
        Pos_phi_psi = rotate_DA(Pos_phi.copy(), setPsi, delta_term_psi, psi_index, moveAtomID_psi)
        #Pos_b4_chi1 = Pos_phi_psi.copy()
        for chi1_loop in range(12, 13):
            #Position = Pos_b4_chi1.copy()
            setChi = chi1_loop * 5.0
            Position_ppc = rotate_DA(Pos_phi_psi.copy(), setChi, delta_term_chi1, chi1_index, moveAtomID_chi1)
            total_E = 0
            # Now that we've rotated everything, check for clashes
            diff_pos = Position_ppc[clash_list[:, 0], :] - Position_ppc[clash_list[:, 1], :]
            sum_2 = np.sum(np.square(diff_pos), 1)
            ind0 = sum_2 < radii_2

            s_r_6 = np.power(radii_2[ind0] / sum_2[ind0], 3)
            E = np.power(1 - s_r_6, 2)
            total_E = np.sum(E)
           
            clash_used = clash_list[ind0, :]
            ind1 = np.isin(moveAtomID_CH3_1, clash_used)
            ind2 = np.isin(moveAtomID_CH3_2, clash_used)

            # Try rotating CH3
            if (total_E > 0 and (np.any(ind1) or np.any(ind2))):
                #print(setPhi, setPsi, setChi)
                min_E = total_E
                #Pos_b4_CH3 = Position.copy()
                # Rotate CH3 groups
                min_E = rotate_CH3(Position_ppc.copy(), delta_term_CH3_1, delta_term_CH3_2, CH3_1_index, CH3_2_index, moveAtomID_CH3_1, moveAtomID_CH3_2, clash_list, radii_2, min_E, CH3_clash_list, CH3_radii_list)
                total_E = min_E
            # Try rotating terminal CH3
            # if (total_E > 0 and (np.any(ind1) or np.any(ind2)) and total_E < 10):
            #     #print(setPhi, setPsi, setChi)
            #     min_E = total_E
            #     Pos_b4_CH3 = Position.copy()
            #     # Rotate CH3 groups
            #     min_E = rotate_CH3(Position, delta_term_CH3_1, delta_term_CH3_2, CH3_1_index, CH3_2_index, moveAtomID_CH3_1, moveAtomID_CH3_2, clash_list, radii_2, min_E)
            #     total_E = min_E

            all_energy[counter, 0] = setPhi
            all_energy[counter, 1] = setPsi
            all_energy[counter, 2] = setChi
            all_energy[counter, 3] = total_E
            counter = counter + 1
            # Figure out if you should add CH3 rotation

print('min', min_E)
np.savetxt('Val_energy_180_new_0330.txt', all_energy, fmt='%d %d %d %7.4f')

# fig = plt.figure()
# ax = plt.axes(projection='3d')
# bb = [0, 1, 4, 6, 8, 20, 22, 24, 25]
# ax.plot3D(Position[bb, 0], Position[bb, 1], Position[bb, 2])
# ax.set_aspect('equal', 'box')
# plt.show()