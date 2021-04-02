import numpy as np
import matplotlib.pyplot as plt


all_bonds = np.loadtxt('MD_Data/all_bonds.txt')
all_angles = np.loadtxt('MD_Data/all_angles.txt')
bond_angle_means = np.loadtxt('bond_angle_means.txt')
bond_length_means = np.loadtxt('bond_length_means.txt')
bond_dihedrals = np.loadtxt('MD_Data/all_dihedral.txt')
data = np.loadtxt('MD_Data/MD_energy_1.txt')

ind0 = bond_dihedrals[7, :] < 0
bond_dihedrals[7, ind0] = bond_dihedrals[7, ind0] + 360
ind0 = bond_dihedrals[8, :] < 0
bond_dihedrals[8, ind0] = bond_dihedrals[8, ind0] + 360
plt.hist(bond_dihedrals[7, :])
plt.show()
plt.hist(bond_dihedrals[8, :])
plt.show()
omega1_mean = np.mean(bond_dihedrals[7, :])
omega1_std = np.std(bond_dihedrals[7, :])
omega2_mean = np.mean(bond_dihedrals[8, :])
omega2_std = np.std(bond_dihedrals[8, :])
print(omega1_mean, omega1_std, omega2_mean, omega2_std)

bond_min_use = bond_length_means[:, 0] - 3.0 * bond_length_means[:, 1]
bond_max_use = bond_length_means[:, 0] + 3.0 * bond_length_means[:, 1]

angle_min_use = bond_angle_means[:, 0] - 3.0 * bond_angle_means[:, 1]
angle_max_use = bond_angle_means[:, 0] + 3.0 * bond_angle_means[:, 1]

print('number of samples', data.shape[0])

counter = 0
to_use_array = np.zeros(data.shape[0])
for i in range(0, data.shape[0]):
    # ind0 = bond_dihedrals[7, i] < omega1_mean - omega1_std * 3.0
    # ind1 = bond_dihedrals[7, i] > omega1_mean + omega1_std * 3.0
    # ind2 = bond_dihedrals[8, i] < omega2_mean - omega2_std * 3.0
    # ind3 = bond_dihedrals[8, i] > omega2_mean + omega2_std * 3.0
    # if (ind0 or ind1 or ind2 or ind3):
    #     to_use_array[i] = 0
    #     counter = counter + 1
    # else:
    this_bonds = all_bonds[:, i]
    this_angles = all_angles[:, i]
    ind0 = this_bonds <= bond_min_use
    ind1 = this_bonds >= bond_max_use
    ind2 = this_angles <= angle_min_use
    ind3 = this_angles >= angle_max_use
    if (ind0.any() or ind1.any() or ind2.any() or ind3.any()):
        to_use_array[i] = 0
        counter = counter + 1
    else:
        to_use_array[i] = 1

print(counter)

ind0 = to_use_array == 1
used_energy = data[ind0]
discarded_energy = data[~ind0]
print(used_energy.shape)
print(discarded_energy.shape)
fig, axs = plt.subplots(1, 2, tight_layout=True)
# We can set the number of bins with the `bins` kwarg
axs[0].hist(used_energy)
axs[1].hist(discarded_energy)
plt.show()

fig, axs = plt.subplots(2, 1, tight_layout=True)
ind0 = data > 20
sub_angles = all_angles[:, 0:data.shape[0]]
axs[0].boxplot(sub_angles[:, ind0].T)

ind0 = data > 20

axs[1].boxplot(sub_angles[:, ~ind0].T)
plt.show()