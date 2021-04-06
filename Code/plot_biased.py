import numpy as np  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import math

folder = '/Users/jmorte02/Documents/Projects/Dipeptides/Val_dipeptide/s1/bemeta/biased/'
data = np.loadtxt(folder + 'prod_angle_s1.xvg', skiprows=17)
nchi = 1
AA = 'Val'

print(data.shape[0])

# Loop over all residues in the cyclic peptide
for val_loop in range(0, 1):
    if (val_loop == 0): 
        Val_chi1 = data[:, 4].copy()
        if (nchi >= 2):
            Val_chi2 = data[:, 5].copy()
        else:
            Val_chi2 = np.zeros([1,1])
        Val_phi = data[:, 2].copy()
        Val_psi = data[:, 3].copy()
    plt.plot(Val_psi, '.')
    plt.show()

    # Make chi values between 0 and 360
    ind0 = Val_chi1 < 0
    Val_chi1[ind0] = Val_chi1[ind0] + 360

    if (nchi >= 2):
        ind0 = Val_chi2 < 0
        Val_chi2[ind0] = Val_chi2[ind0] + 360
    # if (nchi >= 3):
    #     ind0 = Val_chi3 < 0
    #     Val_chi3[ind0] = Val_chi3[ind0] + 360
    # if (nchi >= 4):
    #     ind0 = Val_chi4 < 0
    #     Val_chi4[ind0] = Val_chi4[ind0] + 360

    num_steps = Val_chi1.shape[0]

    # Plot phi/psi heatmap of full trajectory
    square_chi1 = np.zeros([36, 36])
    counts = np.zeros([36, 36])
    max_phi = 0
    max_psi = 0
    for i in range(0, Val_phi.shape[0]):
        phi = Val_phi[i]
        psi = Val_psi[i]
        phi_index = (math.floor(phi / 10) + 17)
        psi_index = (math.floor(psi / 10) + 17)
        square_chi1[psi_index, phi_index] = square_chi1[psi_index, phi_index] + 1
    square_chi1 = square_chi1 / Val_phi.shape[0] / 100
    np.savetxt(folder + AA + '_phipsiAll.txt', square_chi1)
    print(np.sum(square_chi1))

    ind0 = square_chi1 == 0
    square_chi1[ind0] = 'nan'

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)
    ax.set_aspect('equal', adjustable='box')
    plt.imshow(square_chi1, origin='lower', cmap='Blues')
    # plt.plot(beta_x + 17.5, (beta_y + 17.5), 'r')
    # plt.plot(beta_x1 + 17.5, (beta_y1 + 17.5), 'r')
    # plt.plot(pII_x + 17.5, (pII_y + 17.5), 'r')
    # plt.plot(pII_x1 + 17.5, (pII_y1 + 17.5), 'r')
    # plt.plot(alphaL_x + 17.5, (alphaL_y + 17.5), 'r')
    # plt.plot(alphaR_x + 17.5, (alphaR_y + 17.5), 'r')
    plt.colorbar()
    plt.xlabel('$\phi$', fontsize=20)
    plt.ylabel('$\psi$', fontsize=20)

    locs, labels = plt.xticks()            # Get locations and labels
    plt.xticks(np.arange(0, 36, 6), np.arange(-180, 180, 60), fontsize=20)
    locs, labels = plt.yticks()            # Get locations and labels
    plt.yticks(np.arange(0, 36, 6), np.arange(-180, 180, 60), fontsize=20)
    #plt.clim(0, 0.0001)
    plt.axis([-.5, 35.5, -.5, 35.5])
    if (val_loop == 0):
        plt.title(AA, fontsize=20)
        plt.savefig(folder + AA + '_phipsiAll.png', bbox_inches='tight')
    plt.show()

    # Bin phi/psi by chi1 rotamer bins
    ind0 = Val_chi1 < 120
    ind1 = Val_chi1 >= 0
    Val_60_phi = Val_phi[ind0 & ind1]
    Val_60_psi = Val_psi[ind0 & ind1]
    print(Val_60_phi.shape[0] / Val_phi.shape[0])

    ind0 = Val_chi1 < 240
    ind1 = Val_chi1 >= 120
    Val_180_phi = Val_phi[ind0 & ind1]
    Val_180_psi = Val_psi[ind0 & ind1]
    print(Val_180_phi.shape[0] / Val_phi.shape[0])

    ind0 = Val_chi1 < 360
    ind1 = Val_chi1 >= 240
    Val_300_phi = Val_phi[ind0 & ind1]
    Val_300_psi = Val_psi[ind0 & ind1]
    print(Val_300_phi.shape[0] / Val_phi.shape[0])


    # Check for values in the chi1 minimums
    ind0 = Val_chi1 < 10
    ind1 = Val_chi1 >= 0
    val_bad = Val_chi1[ind0 & ind1]
    print(val_bad.shape[0] / Val_phi.shape[0])

    ind0 = Val_chi1 < 130
    ind1 = Val_chi1 >= 110
    val_bad = Val_chi1[ind0 & ind1]
    print(val_bad.shape[0] / Val_phi.shape[0])

    ind0 = Val_chi1 < 250
    ind1 = Val_chi1 >= 230
    val_bad = Val_chi1[ind0 & ind1]
    print(val_bad.shape[0] / Val_phi.shape[0])

    ind0 = Val_chi1 < 370
    ind1 = Val_chi1 >= 350
    val_bad = Val_chi1[ind0 & ind1]
    print(val_bad.shape[0] / Val_phi.shape[0])

    # Plot each 1D chi distribution
    bins = np.arange(0, 360, 5) + 2.5

    for loop_chi in range(0, nchi):
        if (loop_chi == 0):
            this_chi = Val_chi1
            save_str = 'chi1'
        elif (loop_chi == 1):
            this_chi = Val_chi2
            save_str = 'chi2'
        elif (loop_chi == 2):
            this_chi = Val_chi3
            save_str = 'chi3'
        elif (loop_chi == 3):
            this_chi = Val_chi4
            save_str = 'chi4'

        n = np.zeros([72, 1])
        for i in range(0, this_chi.shape[0]):
            chi = this_chi[i]
            chi_index = math.floor(chi / 5)
            n[chi_index] = n[chi_index] + 1
        n = n / this_chi.shape[0]/5

        plt.plot(bins[0:72], n, 'k')
        plt.xticks(np.arange(0, 360, step=60))

        plt.axis([0, 360, 0, 0.03])
        if (loop_chi == 0):
            plt.xlabel('$\chi_1$')
            plt.ylabel('P($\chi_1$)')
        elif (loop_chi == 1):
            plt.xlabel('$\chi_2$')
            plt.ylabel('P($\chi_2$)')
        elif (loop_chi == 2):
            plt.xlabel('$\chi_3$')
            plt.ylabel('P($\chi_3$)')
        elif (loop_chi == 3):
            plt.xlabel('$\chi_4$')
            plt.ylabel('P($\chi_4$)')
        if (val_loop == 0):
            plt.title(AA, fontsize=20)
            plt.savefig(folder + AA + '_P' + save_str + '.png', bbox_inches='tight')
        plt.close()
        plt.close()

    # Plot chi1/chi2 if relevant
    if (nchi >= 2):
        n = np.zeros([36, 36])
        for i in range(0, this_chi.shape[0]):
            chi1 = Val_chi1[i]
            chi2 = Val_chi2[i]
            chi1_index = math.floor(chi1 / 10)
            chi2_index = math.floor(chi2 / 10)
            n[chi2_index, chi1_index] = n[chi2_index, chi1_index] + 1
        n = n / this_chi.shape[0]
        square_chi1 = n.copy()
        ind0 = square_chi1 == 0
        square_chi1[ind0] = 'nan'
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111)
        ax.set_aspect('equal', adjustable='box')
        plt.imshow(square_chi1, origin='lower', cmap='jet')
        plt.colorbar()
        plt.xlabel('$\chi_1$', fontsize=20)
        plt.ylabel('$\chi_2$', fontsize=20)

        locs, labels = plt.xticks()            # Get locations and labels
        plt.xticks(np.arange(0, 36, 6), np.arange(0, 360, 60), fontsize=20)
        locs, labels = plt.yticks()            # Get locations and labels
        plt.yticks(np.arange(0, 36, 6), np.arange(0, 360, 60), fontsize=20)
        plt.axis([-.5, 35.5, -.5, 35.5])
        plt.clim(0, 0.08)

        plt.title(AA, fontsize=20)
        plt.savefig(folder + AA + '_chi1_chi2.png', bbox_inches='tight')

        plt.close()

    # # Plot 60

    # Loop over chi1 rotamers and plot phi/psi

    for chi_rotamers in range(0, 3):
        if (chi_rotamers == 0):
            this_phi = Val_60_phi
            this_psi = Val_60_psi
            save_str = '60'
        elif (chi_rotamers == 1):
            this_phi = Val_180_phi
            this_psi = Val_180_psi
            save_str = '180'
        elif (chi_rotamers == 2):
            this_phi = Val_300_phi
            this_psi = Val_300_psi
            save_str = '300'

        square_chi1 = np.zeros([36, 36])
        counts = np.zeros([36, 36])
        for i in range(0, this_phi.shape[0]):
            phi = this_phi[i]
            phi_index = (math.floor(phi / 10) + 17)
            psi = this_psi[i]
            psi_index = (math.floor(psi / 10) + 17)
            square_chi1[psi_index, phi_index] = square_chi1[psi_index, phi_index] + 1


        if (chi_rotamers == 0):
            square_chi1_60 = square_chi1.copy()
        elif (chi_rotamers == 1):
            square_chi1_180 = square_chi1.copy()
        elif (chi_rotamers == 2):
            square_chi1_300 = square_chi1.copy()
        square_chi1 = square_chi1
        ind0 = square_chi1 == 0
        square_chi1[ind0] = 'nan'


        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111)
        ax.set_aspect('equal', adjustable='box')
        to_plot = square_chi1 / num_steps
        plt.imshow(square_chi1 / num_steps / 100, origin='lower', cmap='jet')
        # plt.plot(beta_x + 17.5, (beta_y + 17.5), 'r')
        # plt.plot(beta_x1 + 17.5, (beta_y1 + 17.5), 'r')
        # plt.plot(pII_x + 17.5, (pII_y + 17.5), 'r')
        # plt.plot(pII_x1 + 17.5, (pII_y1 + 17.5), 'r')
        # plt.plot(alphaL_x + 17.5, (alphaL_y + 17.5), 'r')
        # plt.plot(alphaR_x + 17.5, (alphaR_y + 17.5), 'r')
        plt.colorbar()
        plt.xlabel('$\phi$', fontsize=20)
        plt.ylabel('$\psi$', fontsize=20)

        locs, labels = plt.xticks()            # Get locations and labels
        plt.xticks(np.arange(0, 36, 6), np.arange(-180, 180, 60), fontsize=20)
        locs, labels = plt.yticks()            # Get locations and labels
        plt.yticks(np.arange(0, 36, 6), np.arange(-180, 180, 60), fontsize=20)
        plt.clim(0, 0.0001)
        plt.axis([-.5, 35.5, -.5, 35.5])

        plt.title('$\chi_1$ = ' + save_str + '$^{\circ}$', fontsize=20)

        if (val_loop == 0):
            plt.title(AA, fontsize=20)
            plt.savefig(folder + AA + '_' + save_str + '.png', bbox_inches='tight')
            np.savetxt(folder + AA + '_' + save_str + '.txt', to_plot)
