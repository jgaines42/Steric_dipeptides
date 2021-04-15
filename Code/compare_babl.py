
import numpy as np
import matplotlib.pyplot as plt

n_bins = 100
Atoms = np.genfromtxt('Val_atoms.txt', dtype='str')
angle_list = np.loadtxt('Val_angles.txt', dtype=int)
angle_list = angle_list - 1


all_bonds = np.loadtxt('../../Dun_coordinates/all_bonds.txt')
bond_length_means = np.zeros([all_bonds.shape[0], 4])
for i in range(0, all_bonds.shape[0]):
    bond_length_means[i, 0] = np.mean(all_bonds[i, :])
    bond_length_means[i, 1] = np.std(all_bonds[i, :])
    bond_length_means[i, 2] = np.min(all_bonds[i, :])
    bond_length_means[i, 3] = np.max(all_bonds[i, :])
print(bond_length_means[0, :])

#all_angles = np.loadtxt('MD_data/all_angles.txt')
all_angles = np.loadtxt('../../Dun_coordinates/all_angles.txt')
bond_angle_means = np.zeros([all_angles.shape[0], 4])
for i in range(0, all_angles.shape[0]):
    bond_angle_means[i, 0] = np.mean(all_angles[i, :])
    bond_angle_means[i, 1] = np.std(all_angles[i, :])
    bond_angle_means[i, 2] = np.min(all_angles[i, :])
    bond_angle_means[i, 3] = np.max(all_angles[i, :])


all_bonds = np.loadtxt('../MD_data/all_bonds.txt')
bond_length_means_MD = np.zeros([all_bonds.shape[0], 4])
for i in range(0, all_bonds.shape[0]):
    bond_length_means_MD[i, 0] = np.mean(all_bonds[i, :])
    bond_length_means_MD[i, 1] = np.std(all_bonds[i, :])
    bond_length_means_MD[i, 2] = np.min(all_bonds[i, :])
    bond_length_means_MD[i, 3] = np.max(all_bonds[i, :])

bond_length_means_MD = np.loadtxt('../MD_data/bond_length_means.txt')

#all_angles = np.loadtxt('MD_data/all_angles.txt')
bond_angle_means_MD = np.loadtxt('../MD_data/bond_angle_means.txt')
# for i in range(0, all_angles.shape[0]):
#     bond_angle_means_MD[i, 0] = np.mean(all_angles[i, :])
#     bond_angle_means_MD[i, 1] = np.std(all_angles[i, :])
#     bond_angle_means_MD[i, 2] = np.min(all_angles[i, :])
#     bond_angle_means_MD[i, 3] = np.max(all_angles[i, :])


ind0 = bond_length_means[:, 0] < 3
sub_means = bond_length_means[ind0, :]
sub_MD = bond_length_means_MD[ind0, :]
print(sub_means.shape[0])
plt.errorbar(np.arange(0, sub_means.shape[0]), sub_means[:, 0], yerr=sub_means[:, 1],  fmt='o')
plt.errorbar(np.arange(0, sub_MD.shape[0])+.33, sub_MD[:, 0], yerr=sub_MD[:, 1],  fmt='x')
plt.show()


ind0 = bond_angle_means[:, 0] > 100
sub_means = bond_angle_means[ind0, :]
sub_MD = bond_angle_means_MD[ind0, :]
plt.errorbar(np.arange(1, sub_means.shape[0] + 1), sub_means[:, 0], yerr=sub_means[:, 1],  fmt='o', capsize=2, elinewidth=2)
plt.errorbar(np.arange(1, sub_MD.shape[0]+1) + 0.33, sub_MD[:, 0], yerr=sub_MD[:, 1],  fmt='o',capsize=2, elinewidth=2)
plt.plot([5.7, 5.7], [100, 130], 'k')
plt.plot([10.7, 10.7], [100, 130], 'k')
plt.plot([15.7, 15.7], [100, 130], 'k')
plt.plot([20.7, 20.7], [100, 130], 'k')
plt.plot([25.7, 25.7], [100, 130], 'k')
plt.plot([30.7, 30.7], [100, 130], 'k')
plt.plot([35.7, 35.7], [100, 130], 'k')

plt.show()
print(angle_list[ind0, :])
sub_list = angle_list[ind0, :]
for i in range(0, sub_list.shape[0]):
	print(Atoms[sub_list[i, 0]], Atoms[sub_list[i, 1]], Atoms[sub_list[i, 2]])


