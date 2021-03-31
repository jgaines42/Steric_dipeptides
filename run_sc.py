import numpy as np
from calcDihedral import calcDihedral
from rotate_DA import rotate_DA
from create_clash_list import create_clash_list
from rotate_CH3 import rotate_CH3
import matplotlib.pyplot as plt  # type: ignore
# Load coordinates
Val_data = 'Val_start.txt'

Atoms = np.loadtxt(Val_data, usecols=[1], dtype='str')
Coord = np.loadtxt(Val_data, usecols=[3, 4, 5])  # in nm
bonds = np.loadtxt('Val_bonds.txt')
angles = np.loadtxt('Val_angles.txt')
radii = np.loadtxt('Val_radii.txt')
n_atoms = Coord.shape[0]

# fix bonds and angles to 0 indexing
bonds = bonds - 1
angles = angles - 1

# Create clash list
clash_list = create_clash_list(n_atoms, bonds, angles, radii)
np.savetxt('clash_list.txt', clash_list, fmt='%d')
# Get radii^2
radii_sum = radii[clash_list[:, 0]] + radii[clash_list[:, 1]]
radii_2 = radii_sum * radii_sum

Coord = Coord * 10.0		# Put in Angstroms


# Get index of peptide caps

prior_res_index = [0, 1, 2, 3, 4, 5]
next_res_index = np.arange(n_atoms - 6, n_atoms)


# Initialize dihedral arrays
phi_index = np.array([4, 6, 8, n_atoms - 8])
#print(Atoms[phi_index])
psi_index = np.array([6, 8, n_atoms - 8, next_res_index[0]])
#print(Atoms[psi_index])
chi1_index = np.array([6, 8, 10, 12])
#print(Atoms[chi1_index])
CH3_1_index = np.array([8, 10, 12, 13])
#print(Atoms[CH3_1_index])
CH3_2_index = np.array([8, 10, 16, 17])
#print(Atoms[CH3_2_index])

# Get initial dihedrals
phi_init = calcDihedral(phi_index, Coord)
psi_init = calcDihedral(psi_index, Coord)
chi1_init = calcDihedral(chi1_index, Coord)
CH3_1_init = calcDihedral(CH3_1_index, Coord)
CH3_2_init = calcDihedral(CH3_2_index, Coord)
print(phi_init, psi_init, chi1_init, CH3_1_init, CH3_2_init)

omega_1 = calcDihedral([1, 4, 6, 8], Coord)
omega_2 = calcDihedral([8, 20, 22, 24], Coord)
print('omegas', omega_1, omega_2)

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
# print(Atoms[moveAtomID_phi])
# print(Atoms[moveAtomID_psi])
# print(Atoms[moveAtomID_chi1])
# print(Atoms[moveAtomID_CH3_1])
# print(Atoms[moveAtomID_CH3_2])

delta_term_phi = np.pi * np.sign(phi_init) * phi_init / 180.0
delta_term_psi = np.pi * np.sign(psi_init) * psi_init / 180.0
delta_term_chi1 = np.pi * np.sign(chi1_init) * chi1_init / 180.0
delta_term_CH3_1 = np.pi * np.sign(CH3_1_init) * CH3_1_init / 180.0
delta_term_CH3_2 = np.pi * np.sign(CH3_2_init) * CH3_2_init / 180.0
print('delta', delta_term_phi, delta_term_psi)
# Get initial E
total_E = 0

# Check clashes
diff_pos = Coord[clash_list[:, 0], :] - Coord[clash_list[:, 1], :]
sum_2 = np.sum(np.square(diff_pos), 1)
ind0 = sum_2 < radii_2
s_r_6 = np.power(radii_2[ind0] / sum_2[ind0], 3)
E = np.power(1 - s_r_6, 2)
total_E = total_E + np.sum(E)
print('init E', total_E)
print(clash_list[ind0, :])
print(sum_2[ind0])
print(radii_2[ind0])

Initial_Position = Coord
min_E = 1000000

# All energy data
all_energy = np.zeros((72, 4))
counter = 0
# Loop over all phi
low_E = 1000
for chi1_loop in range(0, 72):
    Position = Initial_Position.copy()
    setChi = chi1_loop * 5.0
    #setChi = 180
    Position = rotate_DA(Position, setChi, delta_term_chi1, chi1_index, moveAtomID_chi1)
    total_E = 0
    # Now that we've rotated everything, check for clashes
    diff_pos = Position[clash_list[:, 0], :] - Position[clash_list[:, 1], :]
    sum_2 = np.sum(np.square(diff_pos), 1)
    ind0 = sum_2 < radii_2

    s_r_6 = np.power(radii_2[ind0] / sum_2[ind0], 3)
    E = np.power(1 - s_r_6, 2)
    total_E = np.sum(E)

    clash_used = clash_list[ind0, :]
    ind1 = np.isin(moveAtomID_CH3_1, clash_used)
    ind2 = np.isin(moveAtomID_CH3_2, clash_used)

    # Try rotating CH3
    if (total_E > 0 and (np.any(ind1) or np.any(ind2)) and total_E < 10):
        #print(setPhi, setPsi, setChi)
        min_E = total_E
        Pos_b4_CH3 = Position.copy()
        # Rotate CH3 groups
        min_E = rotate_CH3(Position, delta_term_CH3_1, delta_term_CH3_2, CH3_1_index, CH3_2_index, moveAtomID_CH3_1, moveAtomID_CH3_2, clash_list, radii_2, min_E)
        total_E = min_E
    if (total_E < low_E):
        low_E = total_E
    all_energy[counter, 0] = phi_init
    all_energy[counter, 1] = psi_init
    all_energy[counter, 2] = setChi
    all_energy[counter, 3] = total_E
    counter = counter + 1
    # Figure out if you should add CH3 rotation
    #np.savetxt('temp.txt', Position/10, fmt='%2.3f')
print('min', low_E)
np.savetxt('Val_energy_rotSC_box.txt', all_energy, fmt='%d %d %d %7.4f')

# fig = plt.figure()
# ax = plt.axes(projection='3d')
# bb = [0, 1, 4, 6, 8, 20, 22, 24, 25]
# ax.plot3D(Position[bb, 0], Position[bb, 1], Position[bb, 2])
# ax.set_aspect('equal', 'box')
# plt.show()