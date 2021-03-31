import numpy as np
import matplotlib.pyplot as plt

n_bins = 100
Atoms = np.genfromtxt('Val_atoms.txt', dtype='str')
angle_list = np.loadtxt('Val_angles.txt', dtype=int)
angle_list = angle_list - 1
# all_dihedral = np.loadtxt('MD_data/all_dihedral.txt')

# ind0 = all_dihedral[5, :] < 0
# all_dihedral[5, ind0] = all_dihedral[5, ind0] + 360
# ind0 = all_dihedral[6, :] < 0
# all_dihedral[6, ind0] = all_dihedral[6, ind0] + 360
# ind0 = all_dihedral[7, :] < 0
# all_dihedral[7, ind0] = all_dihedral[7, ind0] + 360
# ind0 = all_dihedral[8, :] < 0
# all_dihedral[8, ind0] = all_dihedral[8, ind0] + 360

# fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)
# bins1 = np.arange(0, 360, 5)
# # We can set the number of bins with the `bins` kwarg
# axs[0].hist(all_dihedral[5, :], bins=bins1)
# axs[1].hist(all_dihedral[6, :], bins=bins1)
# plt.show()
# print(np.mean(all_dihedral[8, :]), np.std(all_dihedral[8, :]))

# print(np.mean(all_dihedral[7, :]), np.std(all_dihedral[7, :]))

all_bonds = np.loadtxt('MD_data/all_bonds.txt')
# bins1 = np.arange(0.9, 1.7, 0.01)
# bins_H = np.arange(0.99, 1.15, 0.001)

# bond_length_means = np.zeros([all_bonds.shape[0], 4])
# for i in range(0, all_bonds.shape[0]):
#     bond_length_means[i, 0] = np.mean(all_bonds[i, :])
#     bond_length_means[i, 1] = np.std(all_bonds[i, :])
#     bond_length_means[i, 2] = np.min(all_bonds[i, :])
#     bond_length_means[i, 3] = np.max(all_bonds[i, :])
# np.savetxt('MD_data/bond_length_means.txt', bond_length_means, fmt='%6.4f')

all_angles = np.loadtxt('MD_data/all_angles.txt')

bond_angle_means = np.zeros([all_angles.shape[0], 4])
for i in range(0, all_angles.shape[0]):
    bond_angle_means[i, 0] = np.mean(all_angles[i, :])
    bond_angle_means[i, 1] = np.std(all_angles[i, :])
    bond_angle_means[i, 2] = np.min(all_angles[i, :])
    bond_angle_means[i, 3] = np.max(all_angles[i, :])
# np.savetxt('MD_data/bond_angle_means.txt', bond_angle_means, fmt='%5.2f')

# fig1, ax1 = plt.subplots()
# ax1.boxplot(all_angles.T)
# plt.show()

# for i in range(4, 5):
#     plt.hist(all_angles[i, :], 30)
# plt.show()
# print(np.mean(all_angles[4, :]))

# fig = plt.figure()

# # Create an axes instance
# ax = fig.add_axes([0,0,1,1])

# # Create the boxplot
# bp = ax.violinplot(all_angles.T)
# plt.show()

# f, axs = plt.subplots(1, 8, sharey=True)
# for i in range(0, 8):
#     #vert_hist = np.histogram(all_angles[i,:], bins=100)
#     axs[i].hist(all_angles[i, :], bins=100, orientation="horizontal")
#     axs[i].set_title(Atoms[angle_list[i, 0]] + '-' + Atoms[angle_list[i, 1]] + '-' + Atoms[angle_list[i, 2]])
#     for j in range(0, 11):
#         axs[i].plot([0, 40000], [all_angles[i, j], all_angles[i, j]], label=str(j))
# axs[7].legend()
# plt.show()

# f, axs = plt.subplots(1, 8, sharey=True)
# for i in range(0, 8):
#     #vert_hist = np.histogram(all_angles[i,:], bins=100)
#     axs[i].hist(all_angles[i+8, :], bins=100, orientation="horizontal")
#     axs[i].set_title(Atoms[angle_list[i+8, 0]] + '-' + Atoms[angle_list[i+8, 1]] + '-' + Atoms[angle_list[i+8, 2]])
#     for j in range(0, 11):
#         axs[i].plot([0, 40000], [all_angles[i+8, j], all_angles[i+8, j]], label=str(j))
# axs[7].legend()
# plt.show()

min_allowed = bond_angle_means[:, 0] - 3 * bond_angle_means[:, 1]
max_allowed = bond_angle_means[:, 0] + 3 * bond_angle_means[:, 1]

for loop_angle in range(0, 6):
    f, axs = plt.subplots(1, 8, sharey=True)
    for i in range(0, 8):
        #vert_hist = np.histogram(all_angles[i,:], bins=100)
        #axs[i].hist(all_angles[i + (loop_angle * 8), :], bins=100, orientation="horizontal")
        #axs[i].set_title(Atoms[angle_list[i + (loop_angle * 8), 0]] + '-' + Atoms[angle_list[i + (loop_angle * 8), 1]] + '-' + Atoms[angle_list[i + (loop_angle * 8), 2]])
        for j in range(0, 100):
            #axs[i].plot([0, 40000], [all_angles[i + (loop_angle * 8), j], all_angles[i + (loop_angle * 8), j]], label=str(j))
            if (all_angles[i + (loop_angle * 8), j] < min_allowed[i + loop_angle * 8] or all_angles[i + (loop_angle * 8), j] > max_allowed[i + loop_angle * 8]):
                print('sample ' + str(j) + ' angle ' + str(i + (loop_angle * 8)))
    #axs[7].legend()
    #plt.show()

# h_bonding = [0, 1, 2, 6, 8, 11, 14, 15, 16, 15, 18, 19, 22, 24, 25, 26]
# fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
# for i in range(0, 16):
#     a, bin_edges = np.histogram(all_bonds[h_bonding[i], :], bins=bins_H)
#     axs.plot(bins_H[0:bins_H.shape[0]-1], a/all_bonds.shape[1]/0.01)

# for i in range(7, 14):
#     a, bin_edges = np.histogram(all_bonds[i, :], bins=bins1)
#     axs[0][1].plot(bins1[0:bins1.shape[0]-1], a/all_bonds.shape[1]/0.01)

# for i in range(14, 21):
#     a, bin_edges = np.histogram(all_bonds[i, :], bins=bins1)
#     axs[1][0].plot(bins1[0:bins1.shape[0]-1], a/all_bonds.shape[1]/0.01)

# for i in range(21, 27):
#     a, bin_edges = np.histogram(all_bonds[i, :], bins=bins1)
#     axs[1][1].plot(bins1[0:bins1.shape[0]-1], a/all_bonds.shape[1]/0.01)
plt.show()