# create Dunbrack coordinates
import numpy as np
from calcDihedral import calcDihedral
from create_clash_list import create_clash_list
from Bio.PDB import*
import Bio.PDB as pdb

parser=PDBParser()
file1 = '/Users/jmorte02/Documents/Projects/Dun_coordinates/Val_coordinates.txt'
s = parser.get_structure("my_pdb", "Val_example.pdb")
io=PDBIO()
io.set_structure(s)
model = s[0]
chain = model
atoms = [a for a in chain.get_atoms() if pdb.is_aa(a.parent)]


n_sample = 10
# Load bonds
bonds = np.loadtxt('Val_bonds.txt', dtype=int)  
angles = np.loadtxt('Val_angles.txt', dtype=int)
hbonds = np.loadtxt('Val_Hbond.txt', dtype=int)
# fix bonds and angles to 0 indexing
bonds = bonds - 1
angles = angles - 1
hbonds = hbonds - 1

all_bl = np.zeros([bonds.shape[0], n_sample])
all_angles = np.zeros([angles.shape[0], n_sample])
all_dihedral = np.zeros([9, n_sample])
counter = 0

# dihedral positions
phi_index = np.array([4, 6, 8, 20])
psi_index = np.array([6, 8, 20, 22])
chi1_index = np.array([6, 8, 10, 12])
CH3_1_index = np.array([8, 10, 12, 13])
CH3_2_index = np.array([8, 10, 16, 17])
end_CH3_1_index = np.array([0, 1, 4, 6])
end_CH3_2_index = np.array([20, 22, 24, 25])
omega_1_index = np.array([1, 4, 6, 8])
omega_2_index = np.array([8, 20, 22, 24])

radii = np.loadtxt('Val_radii.txt')
# Create clash list
clash_list = create_clash_list(28, bonds, angles, radii)
# Get radii^2
radii_sum = radii[clash_list[:, 0]] + radii[clash_list[:, 1]]


# modify H-bond radii
for i in range(0, hbonds.shape[0]):
    ind1 = np.isin(clash_list[:, 0], hbonds[i, 0])
    ind2 = np.isin(clash_list[:, 1], hbonds[i, 1])
    radii_sum[ind1&ind2] = 1.5
    print(clash_list[ind1&ind2, :])


radii_2 = radii_sum * radii_sum
all_time = np.zeros([n_sample, 1])
all_E = np.zeros([n_sample, 1])

dun_data = np.loadtxt(file1, usecols=np.arange(0, n_sample*3))

for i in range(0, n_sample):
    this_coord = dun_data[:, i * 3 : (i+1) * 3]
    #print(this_coord.shape)
    # Fix atom ordering
    coord = np.zeros([28, 3])
    coord[1:22, :] = this_coord[1:22, :]
    coord[1, :] = this_coord[0, :]
    coord[2:4, :] = coord[0, :]
    coord[4, :] = this_coord[1, :]
    coord[5, :] = this_coord[2, :]
    coord[6, :] = this_coord[3, :]
    coord[7, :] = this_coord[10, :]
    coord[8, :] = this_coord[4, :]
    coord[9, :] = this_coord[11, :]
    coord[10, :] = this_coord[7, :]
    coord[11, :] = this_coord[12, :]
    coord[12, :] = this_coord[8, :]
    coord[13:16, :] = this_coord[13:16, :]
    coord[16, :] = this_coord[9, :]
    coord[17:20, :] = this_coord[16:19, :]
    coord[20:22, :] = this_coord[5:7, :]
    #print(this_coord[19:22, :])
    coord[22, :] = this_coord[19, :]
    coord[23, :] = this_coord[21, :]
    coord[24, :] = this_coord[20, :]
    #coord = coord.T
    #print(coord[0,:])
    # Check clashes
    diff_pos = coord[clash_list[:, 0], :] - coord[clash_list[:, 1], :]
    sum_2 = np.sum(np.square(diff_pos), 1)
    ind0 = sum_2 < radii_2
    ind1 = sum_2 > 0
    s_r_6 = np.power(radii_2[ind0&ind1] / sum_2[ind0&ind1], 3)
    E = np.power(1 - s_r_6, 2)
    total_E = np.sum(E)
    all_E[i] = total_E / 72.0
    #print(clash_list[ind0&ind1,:])
    np.savetxt('../Dun_coordinates/Val_coordinates_' + str(i) + '.txt', coord, fmt = '%6.3f')

    parents = []
    counter = 0;
    for  a in chain.get_atoms() :
        if pdb.is_aa(a.parent):
            parents.append(a.parent)
            counter = counter + 1;
    #xyzs = [(a.coord) for a in atoms]
    #xyzarr = np.array(xyzs)
    xyzarr = this_coord
    id_counter = 1
    # Write to PDB file
    f = open('../Dun_coordinates/Val_coordinates_' + str(i) + '.pdb', 'w')
    for i in range(0, len(atoms)):
        new_res = parents[i].get_id()[1];
        if atoms[i].get_name() == 'N':
            id_counter = id_counter+1
        if len(atoms[i].get_name())<4:
            f.write('{:6s}{:5d}  {:<4}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s} \n'.format('ATOM', i, atoms[i].get_name(), parents[i].get_resname(),atoms[i].get_full_id()[2],  id_counter, '',xyzarr[i][0], xyzarr[i][1], xyzarr[i][2], atoms[i].get_occupancy(), atoms[i].get_bfactor(),atoms[i].get_name()[0] ))
        else:
            f.write('{:6s}{:5d} {:<4} {:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s} \n'.format('ATOM', i, atoms[i].get_name(), parents[i].get_resname(),atoms[i].get_full_id()[2],  id_counter, '',xyzarr[i][0], xyzarr[i][1], xyzarr[i][2], atoms[i].get_occupancy(), atoms[i].get_bfactor(),atoms[i].get_name()[0] ))				
    f.close()

    # # Call reduce to add hydrogen atoms
    # reduce1 = reduce_folder + "./reduce -Trim -quiet " +  folder_name + file_name + "_ordered_l_u.pdb>"+ folder_name + file_name + "_noH.pdb"
    # reduce2 = reduce_folder + "./reduce -quiet " +  folder_name + file_name + "_noH.pdb>" +folder_name + file_name + "_H_l_u.pdb"
    # os.system(reduce1)
    # os.system(reduce2)


#     # calculate bonds
#     diff_pos = coord[bonds[:, 0], :] - coord[bonds[:, 1], :]
#     all_bl[:, i] = np.sqrt(np.sum(np.square(diff_pos), 1))

#     # calculate angles
#     for angle_loop in range(0, angles.shape[0]):
#         ba = coord[angles[angle_loop, 0], :] - coord[angles[angle_loop, 1], :]
#         bc = coord[angles[angle_loop, 2], :] - coord[angles[angle_loop, 1], :]

#         cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
#         all_angles[angle_loop, i] = np.arccos(cosine_angle) * 180 / np.pi

#     # calcualte dihedrals
#     all_dihedral[0, i] = calcDihedral(phi_index, coord)
#     all_dihedral[1, i] = calcDihedral(psi_index, coord)
#     all_dihedral[2, i] = calcDihedral(chi1_index, coord)
#     all_dihedral[3, i] = calcDihedral(CH3_1_index, coord)
#     all_dihedral[4, i] = calcDihedral(CH3_2_index, coord)
#     all_dihedral[5, i] = calcDihedral(end_CH3_1_index, coord)
#     all_dihedral[6, i] = calcDihedral(end_CH3_2_index, coord)
#     all_dihedral[7, i] = calcDihedral(omega_1_index, coord)
#     all_dihedral[8, i] = calcDihedral(omega_2_index, coord)

# np.savetxt('../Dun_coordinates/Energy.txt', all_E)
# np.savetxt('../Dun_coordinates/all_bonds.txt', all_bl, fmt='%5.2f')
# np.savetxt('../Dun_coordinates/all_angles.txt', all_angles, fmt='%6.2f')
# np.savetxt('../Dun_coordinates/all_dihedral.txt', all_dihedral, fmt='%7.2f')