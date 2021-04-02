# create Dunbrack coordinates
import numpy as np
from calcDihedral import calcDihedral
from create_clash_list import create_clash_list

n_sample = 1000
# Load bonds
bonds = np.loadtxt('Val_bonds.txt', dtype=int)  
angles = np.loadtxt('Val_angles.txt', dtype=int)
hbonds = np.loadtxt('Val_Hbond.txt', dtype=int)
# fix bonds and angles to 0 indexing
bonds = bonds - 1
angles = angles - 1
hbonds = hbonds - 1

all_bl = np.zeros([bonds.shape[0], n_sample])
all_angles = np.zeros([angles.shape[0], n_sample])
all_dihedral = np.zeros([9, n_sample])
file1 = '/Users/jmorte02/Documents/Projects/Dun_coordinates/Val_coordinates.txt'
counter = 0

# dihedral positions
phi_index = np.array([4, 6, 8, 20])
psi_index = np.array([6, 8, 20, 22])
chi1_index = np.array([6, 8, 10, 12])
CH3_1_index = np.array([8, 10, 12, 13])
CH3_2_index = np.array([8, 10, 16, 17])
end_CH3_1_index = np.array([0, 1, 4, 6])
end_CH3_2_index = np.array([20, 22, 24, 25])
omega_1_index = np.array([1, 4, 6, 8])
omega_2_index = np.array([8, 20, 22, 24])

radii = np.loadtxt('Val_radii.txt')
# Create clash list
clash_list = create_clash_list(28, bonds, angles, radii)
# Get radii^2
radii_sum = radii[clash_list[:, 0]] + radii[clash_list[:, 1]]


# modify H-bond radii
for i in range(0, hbonds.shape[0]):
    ind1 = np.isin(clash_list[:, 0], hbonds[i, 0])
    ind2 = np.isin(clash_list[:, 1], hbonds[i, 1])
    radii_sum[ind1&ind2] = 1.5
    print(clash_list[ind1&ind2, :])


radii_2 = radii_sum * radii_sum
all_time = np.zeros([n_sample, 1])
all_E = np.zeros([n_sample, 1])
with open(file1, 'rb') as inputfile:
    for i in range(0, n_sample):
        header = read_trr_header(inputfile)
        time1 = header['time']
        data = read_trr_data(inputfile, header)
        coord = data['x'] * 10.0

        # Check clashes
        diff_pos = coord[clash_list[:, 0], :] - coord[clash_list[:, 1], :]
        sum_2 = np.sum(np.square(diff_pos), 1)
        ind0 = sum_2 < radii_2
        s_r_6 = np.power(radii_2[ind0] / sum_2[ind0], 3)
        E = np.power(1 - s_r_6, 2)
        total_E = np.sum(E)
        all_E[i] = total_E / 72.0
        all_time[i] = time1
        np.savetxt('coord_data/Val_coordinates_' + str(i) + '.txt', coord, fmt = '%6.3f')