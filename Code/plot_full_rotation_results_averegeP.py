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
all_probs = [[0] * 36 for i in range(36)]

folder = '../rotation_results_bemeta/Val/'
runs = 100
Te = 0.01
save_tag = '042221'
prob_spots = 0
save_this_chi = np.zeros([36*36, 36])
for run_loop in range(0, runs):
    Val_data = np.loadtxt(folder + 'Val_energy_' + save_tag + '_' + str(run_loop) + '.txt')

    run_max_60 = np.zeros([36, 36])
    run_max_180 = np.zeros([36, 36])
    run_max_300 = np.zeros([36, 36])

    this_60 = np.zeros([36, 36])
    this_180 = np.zeros([36, 36])
    this_300 = np.zeros([36, 36])


    counter = 0
    counter1 = 0
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
            all_probs[phi_loop][psi_loop] = all_probs[phi_loop][psi_loop] + this_prob
            # plt.plot(this_prob)
            save_this_chi[counter1, :] = this_prob/runs + save_this_chi[counter1, :]
            counter1 = counter1 + 1
        # plt.close()
#save_this_chi = save_this_chi / runs
ind0 = save_this_chi < 10E-20
save_this_chi[ind0] = 0
np.savetxt('this_chi.txt', save_this_chi)

# Now extrac tto desired arrays
for phi_loop in range(0, 36):
    for psi_loop in range(0, 36):
        this_data = all_probs[phi_loop][psi_loop]
        psi = psi_loop
        phi = phi_loop
        if (phi >= 18):
            phi = phi - 36
        if (psi >= 18):
            psi = psi - 36
        phi = phi + 18
        psi = psi + 18
        all_60[psi, phi] = np.sum(this_data[0:12])
        all_180[psi, phi] = np.sum(this_data[12:24])
        all_300[psi, phi] = np.sum(this_data[24:36])

        max1 = np.max(this_data[0:12])
        result = np.where(this_data[0:12] == max1)
        result = result[0]
        ind_max = result[int(len(result) / 2)]

        if (all_60[psi, phi] / runs > 0.1):
            all_max_60[psi, phi] = ind_max * 10
            if (ind_max * 10 == 0):
                print('60 ', phi_loop, psi_loop)
                print('60 ', phi, psi)
                print(np.sum(this_data / runs))
                plt.plot(this_data / runs)
                plt.close()
        else:
            all_max_60[psi, phi] = -1

        max1 = np.max(this_data[12:24])
        result = np.where(this_data[12:24] == max1)
        result = result[0]
        ind_max = result[int(len(result) / 2)]

        if (all_180[psi, phi] / runs > 0.1):
            all_max_180[psi, phi] = (ind_max + 12) * 10
            if (ind_max * 10 == 0):
                print('180 ', phi_loop, psi_loop)
                print('180 ', phi, psi)
                print(np.sum(this_data / runs))
                plt.plot(this_data/ runs)
                plt.close()
        else:
            all_max_180[psi, phi] = -1

        max1 = np.max(this_data[24:36])
        result = np.where(this_data[24:36] == max1)
        result = result[0]
        ind_max = result[int(len(result) / 2)]
        if (all_300[psi, phi] / runs > 0.1):
            all_max_300[psi, phi] = (ind_max + 24) * 10
            if (all_max_300[psi, phi] == 0):
                print(ind_max)
                print(result)
                wesdfa
            if (ind_max * 10 == 0):
                print('300 ', phi_loop, psi_loop)
                print('300 ', phi, psi)
                print(np.sum(this_data / runs))
                plt.plot(this_data/ runs)
                plt.close()
        else:
            all_max_300[psi, phi] = -1
        if (all_max_300[psi, phi] == 0):
            print(ind_max)
            print(result)
            wesdfa
        min_1 = this_data[0:6].min() / runs
        min_2 = this_data[6:18].min() / runs
        min_3 = this_data[18:30].min() / runs
        min_4 = this_data[30:36].min() / runs
        # if (min_1 > 0.001 or min_2 > 0.001 or min_3 > 0.001 or min_4 > 0.001):
        #     print("problem ", run_loop, phi_loop, psi_loop, min_1, min_2, min_3, min_4)
            
        #     plt.plot(this_data / runs)
        #     plt.plot([0, 36], [0, 0])
        #     plt.plot([11, 11], [0, .1])
        #     plt.plot([23, 23], [0, .1])
        #     plt.close()
#asfd

a = all_60 + all_180 + all_300


fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)
ax.set_aspect('equal', adjustable='box')
ind0 = all_60 == 0
all_60[ind0] = 'nan'
i2 = plt.imshow(all_60 / runs, origin='lower', cmap='Blues')
plt.colorbar()
plt.clim([0, 1])
plt.title('P($\phi, \psi$)', fontsize=20)
plt.xlabel('$\phi$', fontsize=20)
plt.ylabel('$\psi$', fontsize=20)
locs, labels = plt.xticks()            # Get locations and labels
plt.xticks(np.arange(0, 37, 6), np.arange(-180, 181, 60), fontsize=20)
locs, labels = plt.yticks()            # Get locations and labels
plt.yticks(np.arange(0, 37, 6), np.arange(-180, 181, 60), fontsize=20)
plt.axis([-0.5, 35.5, -0.5, 35.5])
plt.savefig(folder + 'full_rotationVal_steric_60_' + save_tag + '_all' + str(runs) + '_new.png', bbox_inches='tight')
plt.close()

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
plt.clim([0, 1])
plt.savefig(folder + 'Val_steric_180_' + save_tag + '_all' + str(runs) + '_new.png', bbox_inches='tight')
plt.close()

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
plt.clim([0, 1])
plt.savefig(folder + 'Val_steric_300_' + save_tag + '_all' + str(runs) + '_new.png', bbox_inches='tight')
plt.close()



cmap1 = plt.get_cmap('seismic', 13)
# set limits .5 outside true range
cmap1.set_bad(color='gray')


fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)
ax.set_aspect('equal', adjustable='box')
ind0 = all_max_60 == -1
all_max_60[ind0] = 'nan'


i2 = plt.imshow(all_max_60 - 60, origin='lower', cmap=cmap1)

plt.title('P($\phi, \psi$)', fontsize=20)
plt.xlabel('$\phi$', fontsize=20)
plt.ylabel('$\psi$', fontsize=20)
locs, labels = plt.xticks()            # Get locations and labels
plt.xticks(np.arange(0, 37, 6), np.arange(-180, 181, 60), fontsize=20)
locs, labels = plt.yticks()            # Get locations and labels
plt.yticks(np.arange(0, 37, 6), np.arange(-180, 181, 60), fontsize=20)
plt.clim([-65, 65])
plt.colorbar(ticks=np.arange(-60, 70, 10))
plt.savefig(folder + 'Val_steric_60_max_' + save_tag + '_all' + str(runs) + '_new.png', bbox_inches='tight')
plt.close()

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)
ax.set_aspect('equal', adjustable='box')
ind0 = all_max_180 == -1
all_max_180[ind0] = 'nan'
i2 = plt.imshow(all_max_180 - 180, origin='lower', cmap=cmap1)
plt.title('P($\phi, \psi$)', fontsize=20)
plt.xlabel('$\phi$', fontsize=20)
plt.ylabel('$\psi$', fontsize=20)
locs, labels = plt.xticks()            # Get locations and labels
plt.xticks(np.arange(0, 37, 6), np.arange(-180, 181, 60), fontsize=20)
locs, labels = plt.yticks()            # Get locations and labels
plt.yticks(np.arange(0, 37, 6), np.arange(-180, 181, 60), fontsize=20)
plt.clim([-65, 65])
plt.colorbar(ticks=np.arange(-60, 70, 10))
plt.savefig(folder + 'Val_steric_180_max_' + save_tag + '_all' + str(runs) + '_new.png', bbox_inches='tight')
plt.close()

np.savetxt('all_max_300.txt', all_max_300)

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)
ax.set_aspect('equal', adjustable='box')
ind0 = all_max_300 == -1
all_max_300 = all_max_300 - 300
all_max_300[ind0] = 'nan'
i2 = plt.imshow(all_max_300, origin='lower', cmap=cmap1)
plt.title('P($\phi, \psi$)', fontsize=20)
plt.xlabel('$\phi$', fontsize=20)
plt.ylabel('$\psi$', fontsize=20)
locs, labels = plt.xticks()            # Get locations and labels
plt.xticks(np.arange(0, 37, 6), np.arange(-180, 181, 60), fontsize=20)
locs, labels = plt.yticks()            # Get locations and labels
plt.yticks(np.arange(0, 37, 6), np.arange(-180, 181, 60), fontsize=20)
plt.clim([-65, 65])
plt.colorbar(ticks=np.arange(-60, 70, 10))
plt.savefig(folder + 'Val_steric_300_max_' + save_tag + '_all' + str(runs) + '_new.png', bbox_inches='tight')
plt.close()


