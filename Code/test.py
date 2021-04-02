import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('../Dun_coordinates/Val_coordinates.txt')
print(data[0, 0:5])

# data = np.loadtxt('MD_Data/MD_energy_1.txt')
# dihedral = np.loadtxt('MD_Data/all_dihedral.txt')
# MD_energy = np.loadtxt('../Dipeptides/Val_dipeptide/s1/prod_energy.xvg', skiprows=37)
# print(data.shape)
# print(dihedral.shape)
# # plt.plot(MD_energy[:, 0], MD_energy[:, 1])
# # plt.show()

# # plt.hist(MD_energy[:, 1])
# # plt.show()
# # plt.hist(MD_energy[:, 8])
# # plt.show()
# # plt.hist(MD_energy[:, 9])
# # plt.show()
# ind0 = dihedral[2, :] < 0
# dihedral[2, ind0] = dihedral[2, ind0] + 360
# plt.plot(dihedral[2, 0:100000], data[0:100000], '.')
# plt.show()
# ind0 = data > .5
# x = data[ind0]
# sub_dih = dihedral[2, 0:100000]

# print(x.shape[0]/data.shape[0])
# ind0 = data > 1
# x = data[ind0]
# print(x.shape[0]/data.shape[0])
# print(data.shape)
# print(data.max())
# #bins = np.arange(0, 40, 2)
# plt.hist(data)
# plt.show()

# P = np.exp(-1.0 * data / 0.01)
# plt.hist(P)
# plt.show()

# plt.plot(P, '.')
# plt.show()
# print(np.exp(-1.0 * 1 / 0.01))
# print(np.exp(-1.0 *  data.max() / 0.01))

# plt.plot(MD_energy[0:50000,9], data[0:100000:2], '.')
# plt.show()

# plt.plot(MD_energy[0:50000,8], data[0:100000:2], '.')
# plt.show()

# plt.plot(MD_energy[0:50000,3], data[0:100000:2], '.')
# plt.show()

# fig, axs = plt.subplots(nrows=1, ncols=1)

# axs.hist2d(MD_energy[0:50000,9], data[0:100000:2], bins=100)
# plt.show()


# bond_lengths = np.loadtxt('MD_data/all_bonds.txt')
# bond_angles = np.loadtxt('MD_data/all_angles.txt')
# dihedrals = np.loadtxt('MD_data/all_dihedral.txt')
# ind0 = data > 5.0
# sub_angles = bond_angles[:, 0:data.shape[0]]
# sub_lengths = bond_lengths[:, 0:data.shape[0]]
# sub_dih = dihedrals[:, 0:data.shape[0]]
# sub_angles = sub_angles[:, ind0].T
# sub_lengths = sub_lengths[:, ind0].T
# sub_dih = sub_dih[:, ind0].T
# np.savetxt('highE_lengths.txt', sub_lengths)
# np.savetxt('highE_angles.txt', sub_angles)
# np.savetxt('highE_dihedral.txt', sub_dih)

# for i in range(0, bond_lengths.shape[0]):
# 	plt.plot(bond_lengths[i,0:100000], data, '.')

# 	plt.show()