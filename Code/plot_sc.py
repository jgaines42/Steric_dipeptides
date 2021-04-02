import numpy as np
import matplotlib.pyplot as plt  # type: ignore
Val_data = np.loadtxt('Val_energy_rotSC_box.txt')

plt.plot(Val_data[0:72, 2], Val_data[0:72, 3])
plt.show()
P = np.exp(-1.0 * Val_data[0:72, 3] / 0.01)
plt.plot(Val_data[0:72, 2], P)
plt.show()
# print(np.sum(P))
# ind0 = Val_data[:, 2] == 60

# data_60 = Val_data[ind0, :]
# print(data_60.shape)
# phi_psi_60 = np.zeros([72, 72])
# np.savetxt('Val_60.txt', data_60)
# #data_60[0, 3] = 0
# #data_60[1000, 3] = 0
# print(data_60[1000, :])
# for i in range(0, data_60.shape[0]):
#     phi = data_60[i, 0]
#     psi = data_60[i, 1]
#     if (phi >= 180):
#         phi = phi - 360
#     if (psi >= 180):
#         psi = psi - 360
#     phi = int(phi / 5) + 36
#     psi = int(psi / 5) + 36
#     # if (phi >= 72):
#     #     phi = phi - 72
#     # if (psi >= 72):
#     #     psi = psi - 72
#     phi_psi_60[psi, phi] = data_60[i, 3]
#     if (data_60[i, 3] == 0):
#         print(phi, psi)

# #phi_psi_60[10, 20] = 0
# fig = plt.figure(figsize=(10, 8))
# ax = fig.add_subplot(111)
# ax.set_aspect('equal', adjustable='box')

# i2 = plt.imshow(phi_psi_60, origin='lower', cmap='viridis_r')
# plt.colorbar()

# plt.clim(0, 1)
# plt.show()


# P = np.exp(-1.0 * phi_psi_60 / 0.01)
# print(np.sum(P))
# fig = plt.figure(figsize=(10, 8))
# ax = fig.add_subplot(111)
# ax.set_aspect('equal', adjustable='box')

# i2 = plt.imshow(P / np.sum(P), origin='lower', cmap='viridis')

# plt.colorbar()
# plt.show()