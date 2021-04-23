import numpy as np
import matplotlib.pyplot as plt  # type: ignore
folder = '../rotation_results/full_bemeta/Val/'
runs = 1000
Te = 0.01
save_tag = '042221'
all_60 = np.loadtxt(folder + 'Val_steric_60_' + save_tag + '_all' + str(runs) + '_new.txt')
all_180 = np.loadtxt(folder + 'Val_steric_180_' + save_tag + '_all' + str(runs) + '_new.txt')
all_300 = np.loadtxt(folder + 'Val_steric_300_' + save_tag + '_all' + str(runs) + '_new.txt')
all_max_60 = np.loadtxt(folder + 'Val_steric_60_max_' + save_tag + '_all' + str(runs) + '_new.txt')
all_max_180 = np.loadtxt(folder + 'Val_steric_180_max_' + save_tag + '_all' + str(runs) + '_new.txt')
all_max_300 = np.loadtxt(folder + 'Val_steric_300_max_' + save_tag + '_all' + str(runs) + '_new.txt')
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
plt.xticks(np.arange(0, 73, 12), np.arange(-180, 181, 60), fontsize=20)
locs, labels = plt.yticks()            # Get locations and labels
plt.yticks(np.arange(0, 73, 12), np.arange(-180, 181, 60), fontsize=20)
plt.axis([-0.5, 71.5, -0.5, 71.5])
plt.savefig(folder + 'Val_steric_60_' + save_tag + '_all' + str(runs) + '_new.png', bbox_inches='tight')
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
plt.xticks(np.arange(0, 73, 12), np.arange(-180, 181, 60), fontsize=20)
locs, labels = plt.yticks()            # Get locations and labels
plt.yticks(np.arange(0, 73, 12), np.arange(-180, 181, 60), fontsize=20)
plt.axis([-0.5, 71.5, -0.5, 71.5])
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
plt.xticks(np.arange(0, 73, 12), np.arange(-180, 181, 60), fontsize=20)
locs, labels = plt.yticks()            # Get locations and labels
plt.yticks(np.arange(0, 73, 12), np.arange(-180, 181, 60), fontsize=20)
plt.axis([-0.5, 71.5, -0.5, 71.5])
plt.clim([0, 1])
plt.savefig(folder + 'Val_steric_300_' + save_tag + '_all' + str(runs) + '_new.png', bbox_inches='tight')
plt.close()



cmap1 = plt.get_cmap('seismic', 25)
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
plt.xticks(np.arange(0, 73, 12), np.arange(-180, 181, 60), fontsize=20)
locs, labels = plt.yticks()            # Get locations and labels
plt.yticks(np.arange(0, 73, 12), np.arange(-180, 181, 60), fontsize=20)
plt.clim([-62.5, 62.5])
plt.colorbar(ticks=np.arange(-60, 70, 5))
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
plt.xticks(np.arange(0, 73, 12), np.arange(-180, 181, 60), fontsize=20)
locs, labels = plt.yticks()            # Get locations and labels
plt.yticks(np.arange(0, 73, 12), np.arange(-180, 181, 60), fontsize=20)
plt.clim([-62.5, 62.5])
plt.colorbar(ticks=np.arange(-60, 70, 5))
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
plt.xticks(np.arange(0, 73, 12), np.arange(-180, 181, 60), fontsize=20)
locs, labels = plt.yticks()            # Get locations and labels
plt.yticks(np.arange(0, 73, 12), np.arange(-180, 181, 60), fontsize=20)
plt.clim([-62.5, 62.5])
plt.colorbar(ticks=np.arange(-60, 70, 5))
plt.savefig(folder + 'Val_steric_300_max_' + save_tag + '_all' + str(runs) + '_new.png', bbox_inches='tight')
plt.close()