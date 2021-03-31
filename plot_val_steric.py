import numpy as np
import matplotlib.pyplot as plt  # type: ignore
Val_data = np.loadtxt('Val_energy_180_new_0330.txt')

save_tag = '_180_new_0330.txt'
ind0 = Val_data[:, 2] == 180

data_60 = Val_data[ind0, :]
print(data_60.shape)
phi_psi_60 = np.zeros([72, 72])
phi_psi_60 = phi_psi_60 - 100
np.savetxt('Val' + save_tag, data_60)

for i in range(0, data_60.shape[0]):
    phi = data_60[i, 0]
    psi = data_60[i, 1]
    if (phi >= 180):
        phi = phi - 360
    if (psi >= 180):
        psi = psi - 360
    phi = int(phi / 5) + 36
    psi = int(psi / 5) + 36
    # if (phi >= 72):
    #     phi = phi - 72
    # if (psi >= 72):
    #     psi = psi - 72
    phi_psi_60[psi, phi] = data_60[i, 3]
    # if (data_60[i, 3] == 0):
    #     print(phi, psi)

#phi_psi_60[10, 20] = 0
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)
ax.set_aspect('equal', adjustable='box')

i2 = plt.imshow(phi_psi_60, origin='lower', cmap='viridis_r')
plt.colorbar()

plt.clim(0, 1)
plt.show()


P = np.exp(-1.0 * phi_psi_60 / 0.01)
print(np.sum(P))
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)
ax.set_aspect('equal', adjustable='box')

i2 = plt.imshow(P / np.sum(P), origin='lower', cmap='viridis')

plt.colorbar()
plt.show()
np.savetxt('Val_steric' + save_tag, P)
np.savetxt('Val_steric_energy' + save_tag, phi_psi_60)

Val_60 = np.loadtxt('Val_steric_energy_60_new_0330.txt')
Val_180 = np.loadtxt('Val_steric_energy_180_new_0330.txt')
Val_300 = np.loadtxt('Val_steric_energy_300_new_0330.txt')
print(Val_60[60, 25])
print(Val_180[60, 25])
print(Val_300[60, 25])

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)
ax.set_aspect('equal', adjustable='box')
# ind0 = Val_60 == 0
# Val_60[ind0] = 'nan'

i2 = plt.imshow(Val_60, origin='lower', cmap='Blues')
plt.colorbar()
plt.title('P($\phi, \psi$)', fontsize=20)
plt.xlabel('$\phi$', fontsize=20)
plt.ylabel('$\psi$', fontsize=20)
locs, labels = plt.xticks()            # Get locations and labels
plt.xticks(np.arange(0, 73, 12), np.arange(-180, 181, 60), fontsize=20)
locs, labels = plt.yticks()            # Get locations and labels
plt.yticks(np.arange(0, 73, 12), np.arange(-180, 181, 60), fontsize=20)
plt.axis([-0.5, 71.5, -0.5, 71.5])
plt.clim([0, 1])
plt.savefig('Val_steric_60_raw.png', bbox_inches='tight')
plt.show()

for i in range(0, 72):
    for j in range(0, 72):
        if (Val_60[i, j] < Val_180[i, j] and Val_60[i, j] < Val_300[i, j]):
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

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)
ax.set_aspect('equal', adjustable='box')
ind0 = Val_60 == 0
Val_60[ind0] = 'nan'
i2 = plt.imshow(Val_60, origin='lower', cmap='Blues')
plt.colorbar()
plt.title('P($\phi, \psi$)', fontsize=20)
plt.xlabel('$\phi$', fontsize=20)
plt.ylabel('$\psi$', fontsize=20)
locs, labels = plt.xticks()            # Get locations and labels
plt.xticks(np.arange(0, 73, 12), np.arange(-180, 181, 60), fontsize=20)
locs, labels = plt.yticks()            # Get locations and labels
plt.yticks(np.arange(0, 73, 12), np.arange(-180, 181, 60), fontsize=20)
plt.axis([-0.5, 71.5, -0.5, 71.5])
plt.savefig('Val_steric_60_new.png', bbox_inches='tight')
plt.show()

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)
ax.set_aspect('equal', adjustable='box')
ind0 = Val_180 == 0
Val_180[ind0] = 'nan'
i2 = plt.imshow(Val_180, origin='lower', cmap='Blues')
plt.colorbar()
plt.title('P($\phi, \psi$)', fontsize=20)
plt.xlabel('$\phi$', fontsize=20)
plt.ylabel('$\psi$', fontsize=20)
locs, labels = plt.xticks()            # Get locations and labels
plt.xticks(np.arange(0, 73, 12), np.arange(-180, 181, 60), fontsize=20)
locs, labels = plt.yticks()            # Get locations and labels
plt.yticks(np.arange(0, 73, 12), np.arange(-180, 181, 60), fontsize=20)
plt.axis([-0.5, 71.5, -0.5, 71.5])
plt.savefig('Val_steric_180_new.png', bbox_inches='tight')
plt.show()

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)
ax.set_aspect('equal', adjustable='box')
ind0 = Val_300 == 0
Val_300[ind0] = 'nan'
i2 = plt.imshow(Val_300, origin='lower', cmap='Blues')
plt.colorbar()
plt.title('P($\phi, \psi$)', fontsize=20)
plt.xlabel('$\phi$', fontsize=20)
plt.ylabel('$\psi$', fontsize=20)
locs, labels = plt.xticks()            # Get locations and labels
plt.xticks(np.arange(0, 73, 12), np.arange(-180, 181, 60), fontsize=20)
locs, labels = plt.yticks()            # Get locations and labels
plt.yticks(np.arange(0, 73, 12), np.arange(-180, 181, 60), fontsize=20)
plt.axis([-0.5, 71.5, -0.5, 71.5])
plt.savefig('Val_steric_300_new.png', bbox_inches='tight')
plt.show()