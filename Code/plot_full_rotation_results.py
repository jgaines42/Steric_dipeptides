# Use this script to plot the results of phi/psi rotation with chi sampling of 60, 180, 300

import numpy as np
import matplotlib.pyplot as plt  # type: ignore

all_60 = np.zeros([36, 36])
all_180 = np.zeros([36, 36])
all_300 = np.zeros([36, 36])

all_max_60 = np.zeros([36, 36])
all_max_180 = np.zeros([36, 36])
all_max_300 = np.zeros([36, 36])
all_Chi_pos = [0] * 36


folder = '../rotation_results/full_bemeta/Val/'
runs = 100
Te = 0.01
save_tag = '040821'
prob_spots = 0

for run_loop in range(0, runs):
    Val_data = np.loadtxt(folder + 'Val_energy_' + save_tag + '_' + str(run_loop) + '.txt')

    run_max_60 = np.zeros([36, 36])
    run_max_180 = np.zeros([36, 36])
    run_max_300 = np.zeros([36, 36])

    this_60 = np.zeros([36, 36])
    this_180 = np.zeros([36, 36])
    this_300 = np.zeros([36, 36])

    counter = 0
    # For each phi/psi, get the full energy array
    for phi_loop in range(0, 36):
        for psi_loop in range(0, 36):
            this_data = Val_data[counter:counter + 36, :]
            counter = counter + 36
            # Find the lowest energy value, use this to convert everything to probability
            min_E = this_data[:, 3].min()
            this_E = this_data[:, 3]

            A = np.exp(-1.0 * (this_E - min_E) / Te)
            sum1 = np.sum(A)
            this_prob = A / sum1


            # Check that the lowest Prob between 0-60, 60-180, 180-300 and 300-360 is zero
            min_1 = this_prob[0:6].min()
            min_2 = this_prob[6:18].min()
            min_3 = this_prob[18:30].min()
            min_4 = this_prob[30:36].min()

            if (min_1 > 0.0001 or min_2 > 0.0001 or min_3 > 0.0001 or min_4 > 0.0001):
                # print("problem ", run_loop, phi_loop, psi_loop)
                # f, (ax1, ax2) = plt.subplots(2,1)
                # ax1.plot(this_E)
                # ax2.plot(this_prob)
                # plt.show()
                # asfd
                prob_spots = prob_spots + 1
            # Find the highest prob in each region (0-120, 120-240, 240-360)

            max1 = np.max(this_prob[0:12])
            result = np.where(this_prob[0:12] == max1)
            result = result[0]
            ind_max = result[int(len(result) / 2)]
            if (np.max(result) == 11):
                f, (ax1, ax2) = plt.subplots(2, 1)
                ax1.plot(this_E)
                ax2.plot(this_prob)
                plt.show()
                asfd

            max1 = np.max(this_prob[12:24])
            result = np.where(this_prob[12:24] == max1)
            result = result[0] + 12
            ind_max = result[int(len(result) / 2)]

            max1 = np.max(this_prob[24:36])
            result = np.where(this_prob[24:36] == max1)
            result = result[0] + 24
            ind_max = result[int(len(result) / 2)]

            max_60 = np.argmax(this_prob[0:12]) * 10
            max_180 = (np.argmax(this_prob[12:24]) + 12) * 10
            max_300 = (np.argmax(this_prob[24:36]) + 24) * 10

            phi = this_data[0, 0]
            psi = this_data[0, 1]
            if (phi >= 180):
                phi = phi - 360
            if (psi >= 180):
                psi = psi - 360
            phi = int(phi / 10) + 18
            psi = int(psi / 10) + 18

            run_max_60[psi, phi] = max_60
            run_max_180[psi, phi] = max_180
            run_max_300[psi, phi] = max_300

            # Get rotamer probabilities by summing in each region
            this_60[psi, phi] = np.sum(this_prob[0:12])
            this_180[psi, phi] = np.sum(this_prob[12:24])
            this_300[psi, phi] = np.sum(this_prob[24:36])

            # Store in Val_60, Val_180 and Val_300
    all_60 = all_60 + this_60
    all_180 = all_180 + this_180
    all_300 = all_300 + this_300

    all_max_60 = all_max_60 + run_max_60
    all_max_180 = all_max_180 + run_max_180
    all_max_300 = all_max_300 + run_max_300

    # for rot_loop in range(0, 3):

    #     this_chi = (rot_loop * 120) + 60
    #     ind0 = Val_data[:, 2] == this_chi

    #     data_60 = Val_data[ind0, :]
    #     phi_psi_60 = np.zeros([36, 36])
    #     phi_psi_60 = phi_psi_60 - 100

    #     for i in range(0, data_60.shape[0]):
    #         phi = data_60[i, 0]
    #         psi = data_60[i, 1]
    #         if (phi >= 180):
    #             phi = phi - 360
    #         if (psi >= 180):
    #             psi = psi - 360
    #         phi = int(phi / 10) + 18
    #         psi = int(psi / 10) + 18
    #         phi_psi_60[psi, phi] = data_60[i, 3]

    #     if (rot_loop == 0):
    #         Val_60 = phi_psi_60.copy()
    #     elif (rot_loop == 1):
    #         Val_180 = phi_psi_60.copy()
    #     else:
    #         Val_300 = phi_psi_60.copy()

    # for i in range(0, 36):
    #     for j in range(0, 36):

    #         if (Val_60[i, j] == Val_180[i, j] and Val_60[i, j] < Val_300[i, j]):
    #             A = np.exp(-1.0 * (Val_180[i, j] - Val_60[i, j]) / 0.01)
    #             B = np.exp(-1.0 * (Val_300[i, j] - Val_60[i, j]) / 0.01)
    #             sum1 = A + B + 1
    #             Val_60[i, j] = 1 / sum1
    #             Val_180[i, j] = A / sum1
    #             Val_300[i, j] = B / sum1

    #         elif (Val_60[i, j] < Val_180[i, j] and Val_60[i, j] < Val_300[i, j]):
    #             A = np.exp(-1.0 * (Val_180[i, j] - Val_60[i, j]) / 0.01)
    #             B = np.exp(-1.0 * (Val_300[i, j] - Val_60[i, j]) / 0.01)

    #             sum1 = A + B + 1
    #             Val_60[i, j] = 1 / sum1
    #             Val_180[i, j] = A / sum1
    #             Val_300[i, j] = B / sum1

    #         elif (Val_180[i, j] < Val_60[i, j] and Val_180[i, j] < Val_300[i, j]):
    #             A = np.exp(-1.0 * (Val_60[i, j] - Val_180[i, j]) / 0.01)
    #             B = np.exp(-1.0 * (Val_300[i, j] - Val_180[i, j]) / 0.01)

    #             sum1 = A + B + 1

    #             Val_60[i, j] = A / sum1
    #             Val_180[i, j] = 1 / sum1
    #             Val_300[i, j] = B / sum1
    #         else:
    #             A = np.exp(-1.0 * (Val_60[i, j] - Val_300[i, j]) / 0.01)
    #             B = np.exp(-1.0 * (Val_180[i, j] - Val_300[i, j]) / 0.01)

    #             sum1 = A + B + 1
    #             Val_60[i, j] = A / sum1
    #             Val_180[i, j] = B / sum1
    #             Val_300[i, j] = 1 / sum1
    # all_60 = all_60 + Val_60
    # all_180 = all_180 + Val_180
    # all_300 = all_300 + Val_300

print(prob_spots)

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
plt.savefig(folder + 'full_rotationVal_steric_60_' + save_tag + '_all' + str(runs) + '_new.png', bbox_inches='tight')
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
plt.savefig(folder + 'Val_steric_180_' + save_tag + '_all' + str(runs) + '_new.png', bbox_inches='tight')
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
plt.savefig(folder + 'Val_steric_300_' + save_tag + '_all' + str(runs) + '_new.png', bbox_inches='tight')
plt.show()


fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)
ax.set_aspect('equal', adjustable='box')
ind0 = all_60 == 0
all_60[ind0] = 'nan'
i2 = plt.imshow(all_max_60 / runs, origin='lower', cmap='Blues')
plt.colorbar()
plt.title('P($\phi, \psi$)', fontsize=20)
plt.xlabel('$\phi$', fontsize=20)
plt.ylabel('$\psi$', fontsize=20)
locs, labels = plt.xticks()            # Get locations and labels
plt.xticks(np.arange(0, 37, 6), np.arange(-180, 181, 60), fontsize=20)
locs, labels = plt.yticks()            # Get locations and labels
plt.yticks(np.arange(0, 37, 6), np.arange(-180, 181, 60), fontsize=20)
plt.show()

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)
ax.set_aspect('equal', adjustable='box')
ind0 = all_60 == 0
all_60[ind0] = 'nan'
i2 = plt.imshow(all_max_180 / runs, origin='lower', cmap='Blues')
plt.colorbar()
plt.title('P($\phi, \psi$)', fontsize=20)
plt.xlabel('$\phi$', fontsize=20)
plt.ylabel('$\psi$', fontsize=20)
locs, labels = plt.xticks()            # Get locations and labels
plt.xticks(np.arange(0, 37, 6), np.arange(-180, 181, 60), fontsize=20)
locs, labels = plt.yticks()            # Get locations and labels
plt.yticks(np.arange(0, 37, 6), np.arange(-180, 181, 60), fontsize=20)
plt.show()

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)
ax.set_aspect('equal', adjustable='box')
ind0 = all_60 == 0
all_60[ind0] = 'nan'
i2 = plt.imshow(all_max_300 / runs, origin='lower', cmap='Blues')
plt.colorbar()
plt.title('P($\phi, \psi$)', fontsize=20)
plt.xlabel('$\phi$', fontsize=20)
plt.ylabel('$\psi$', fontsize=20)
locs, labels = plt.xticks()            # Get locations and labels
plt.xticks(np.arange(0, 37, 6), np.arange(-180, 181, 60), fontsize=20)
locs, labels = plt.yticks()            # Get locations and labels
plt.yticks(np.arange(0, 37, 6), np.arange(-180, 181, 60), fontsize=20)
plt.show()


