# Use this script to plot the results of phi/psi rotation with chi sampling of 60, 180, 300

import numpy as np
import matplotlib.pyplot as plt  # type: ignore

all_60 = np.zeros([36, 36])
all_180 = np.zeros([36, 36])
all_300 = np.zeros([36, 36])
runs = 100
save_tag = '040721'
for run_loop in range(0, runs):
    Val_data = np.loadtxt('../rotation_results/full_rotation/Val_energy_' + save_tag + '_' + str(run_loop) + '.txt')


    # For aech phi/psi, get the full energy array

    # Find the lowest energy value, use this to convert everything to probability

    # Check that the lowest Prob between 0-60, 60-180, 180-300 and 300-360 is zero

    # Find the highest prob in each region (0-120, 120-240, 240-360)

    # Get rotamer probabilities by summing in each region

    # Store in Val_60, Val_180 and Val_300

    for rot_loop in range(0, 3):

        this_chi = (rot_loop * 120) + 60
        ind0 = Val_data[:, 2] == this_chi

        data_60 = Val_data[ind0, :]
        phi_psi_60 = np.zeros([36, 36])
        phi_psi_60 = phi_psi_60 - 100

        for i in range(0, data_60.shape[0]):
            phi = data_60[i, 0]
            psi = data_60[i, 1]
            if (phi >= 180):
                phi = phi - 360
            if (psi >= 180):
                psi = psi - 360
            phi = int(phi / 10) + 18
            psi = int(psi / 10) + 18
            phi_psi_60[psi, phi] = data_60[i, 3]

        if (rot_loop == 0):
            Val_60 = phi_psi_60.copy()
        elif (rot_loop == 1):
            Val_180 = phi_psi_60.copy()
        else:
            Val_300 = phi_psi_60.copy()

    for i in range(0, 36):
        for j in range(0, 36):

            if (Val_60[i, j] == Val_180[i, j] and Val_60[i, j] < Val_300[i, j]):
                A = np.exp(-1.0 * (Val_180[i, j] - Val_60[i, j]) / 0.01)
                B = np.exp(-1.0 * (Val_300[i, j] - Val_60[i, j]) / 0.01)
                sum1 = A + B + 1
                Val_60[i, j] = 1 / sum1
                Val_180[i, j] = A / sum1
                Val_300[i, j] = B / sum1

            elif (Val_60[i, j] < Val_180[i, j] and Val_60[i, j] < Val_300[i, j]):
                A = np.exp(-1.0 * (Val_180[i, j] - Val_60[i, j]) / 0.01)
                B = np.exp(-1.0 * (Val_300[i, j] - Val_60[i, j]) / 0.01)

                sum1 = A + B + 1
                Val_60[i, j] = 1 / sum1
                Val_180[i, j] = A / sum1
                Val_300[i, j] = B / sum1

            elif (Val_180[i, j] < Val_60[i, j] and Val_180[i, j] < Val_300[i, j]):
                A = np.exp(-1.0 * (Val_60[i, j] - Val_180[i, j]) / 0.01)
                B = np.exp(-1.0 * (Val_300[i, j] - Val_180[i, j]) / 0.01)

                sum1 = A + B + 1

                Val_60[i, j] = A / sum1
                Val_180[i, j] = 1 / sum1
                Val_300[i, j] = B / sum1
            else:
                A = np.exp(-1.0 * (Val_60[i, j] - Val_300[i, j]) / 0.01)
                B = np.exp(-1.0 * (Val_180[i, j] - Val_300[i, j]) / 0.01)

                sum1 = A + B + 1
                Val_60[i, j] = A / sum1
                Val_180[i, j] = B / sum1
                Val_300[i, j] = 1 / sum1
    all_60 = all_60 + Val_60
    all_180 = all_180 + Val_180
    all_300 = all_300 + Val_300

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)
ax.set_aspect('equal', adjustable='box')
ind0 = all_60 == 0
all_60[ind0] = 'nan'
i2 = plt.imshow(all_60 / runs, origin='lower', cmap='Blues')
plt.colorbar()
plt.title('P($\phi, \psi$)', fontsize=20)
plt.xlabel('$\phi$', fontsize=20)
plt.ylabel('$\psi$', fontsize=20)
locs, labels = plt.xticks()            # Get locations and labels
plt.xticks(np.arange(0, 37, 6), np.arange(-180, 181, 60), fontsize=20)
locs, labels = plt.yticks()            # Get locations and labels
plt.yticks(np.arange(0, 37, 6), np.arange(-180, 181, 60), fontsize=20)
plt.axis([-0.5, 35.5, -0.5, 35.5])
plt.savefig('../rotation_results/full_rotationVal_steric_60_' + save_tag + '_all' + str(runs) + '.png', bbox_inches='tight')
plt.show()

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)
ax.set_aspect('equal', adjustable='box')
ind0 = all_180 == 0
all_180[ind0] = 'nan'
i2 = plt.imshow(all_180 / runs, origin='lower', cmap='Blues')
plt.colorbar()
plt.title('P($\phi, \psi$)', fontsize=20)
plt.xlabel('$\phi$', fontsize=20)
plt.ylabel('$\psi$', fontsize=20)
locs, labels = plt.xticks()            # Get locations and labels
plt.xticks(np.arange(0, 37, 6), np.arange(-180, 181, 60), fontsize=20)
locs, labels = plt.yticks()            # Get locations and labels
plt.yticks(np.arange(0, 37, 6), np.arange(-180, 181, 60), fontsize=20)
plt.axis([-0.5, 35.5, -0.5, 35.5])
plt.savefig('../rotation_results/full_rotation/Val_steric_180_' + save_tag + '_all' + str(runs) + '.png', bbox_inches='tight')
plt.show()

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)
ax.set_aspect('equal', adjustable='box')
ind0 = all_300 == 0
all_300[ind0] = 'nan'
i2 = plt.imshow(all_300 / runs, origin='lower', cmap='Blues')
plt.colorbar()
plt.title('P($\phi, \psi$)', fontsize=20)
plt.xlabel('$\phi$', fontsize=20)
plt.ylabel('$\psi$', fontsize=20)
locs, labels = plt.xticks()            # Get locations and labels
plt.xticks(np.arange(0, 37, 6), np.arange(-180, 181, 60), fontsize=20)
locs, labels = plt.yticks()            # Get locations and labels
plt.yticks(np.arange(0, 37, 6), np.arange(-180, 181, 60), fontsize=20)
plt.axis([-0.5, 35.5, -0.5, 35.5])
plt.savefig('../rotation_results/full_rotation/Val_steric_300_' + save_tag + '_all' + str(runs) + '.png', bbox_inches='tight')
plt.show()