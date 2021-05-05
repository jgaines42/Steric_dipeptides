/******************************************************************************************************************************
 rotate_amino_acid.cxx

 Samples full phi/psi/chi space of an amino acid. Can change the step size of sampling

 Dependencies: dihedral_functions.cxx, dihedral_functions.h

 Compile: g++ rotate_amino_acid.cxx dihedral_functions.cxx

 Arguements:
 	1. folder:		String of the full path to the folder containing the input coordiantes file
 	2. input_file:	String of the file name of the input coordinate file
 	3. AA_name:		3 letter abbreviation of the amino acid
 	4. coord_number: integer of which set of coordiantes is being used (used in name of saved files)
 	5. save_folder:	full path to folder to save data
 	6. save_name:	additional string to used in the name of saved files (ex. 042221 to specify the date)
******************************************************************************************************************************/

#include "dihedral_functions.h" 
#include "rotate_CH3_group.h"
#include "aa_setup.h" 

#include <iostream>
#include <fstream>
#include <ostream>
#include <cstdlib>
#include <cmath>
#include <string.h>
#include <stdlib.h>     /* srand, rand */


using namespace std;
#define PI 3.14159265


int main(int argc, char **argv){

	/*******************************
	 * Process initial arguments
	 *******************************/
	if (argc < 6){
		cerr << " Invalid input " << endl ;
		cerr << endl;
		return 0;
	}
	
	
	string folder, input_file, AA_name, save_folder, save_name, coord_num;	

	folder = argv[1];		// Get folder name	
	input_file = argv[2];	// Get file name
	AA_name = argv[3];		// Get AA name
	coord_num = argv[4];	// Get which coord
	save_folder= argv[5];	// Get save folder
	save_name= argv[6];		// get same name

	// Initialize variabales
	int return_value = 0;
	int i, j, k;

	

	// Get dihedral angle indexes
	int index_array1[9][4];
	get_dihedral_indexes(AA_name, index_array1);

	int phi_index[4] = {index_array1[0][0], index_array1[0][1], index_array1[0][2], index_array1[0][3]};
	int psi_index[4] = {index_array1[1][0], index_array1[1][1], index_array1[1][2], index_array1[1][3]};
	int chi1_index[4] = {index_array1[2][0], index_array1[2][1], index_array1[2][2], index_array1[2][3]};
	int chi2_index[4] = {index_array1[3][0], index_array1[3][1], index_array1[3][2], index_array1[3][3]};
	int CH3_1_index[4] = {index_array1[5][0], index_array1[5][1], index_array1[5][2], index_array1[5][3]};
	int CH3_2_index[4] = {index_array1[6][0], index_array1[6][1], index_array1[6][2], index_array1[6][3]};
	int end_CH3_1_index[4] = {index_array1[7][0], index_array1[7][1], index_array1[7][2], index_array1[7][3]};
	int end_CH3_2_index[4] = {index_array1[8][0], index_array1[8][1], index_array1[8][2], index_array1[8][3]};


	// Get move atom ID indexes. First value is the actual lenght of the array
	int moveAtomID_phi[30], moveAtomID_psi[30], moveAtomID_chi1[30], moveAtomID_chi2[30], moveAtomID_chi3[30], moveAtomID_CH3_1[30], moveAtomID_CH3_2[30], moveAtomID_end_CH3_1[30], moveAtomID_end_CH3_2[30] ;

	get_move_indexes(AA_name, moveAtomID_phi, moveAtomID_psi, moveAtomID_chi1, moveAtomID_chi2, moveAtomID_chi3, moveAtomID_CH3_1, moveAtomID_CH3_2, moveAtomID_end_CH3_1, moveAtomID_end_CH3_2);

	int n_atoms, clash_size, CH3_clash_size, n_chi;
	if (AA_name.compare("Val") == 0){
		clash_size = 303;		// number of clashes (length of Clash_list_val.txt)
		CH3_clash_size = 129;	// number of clashes involving terminal CH3 atoms (length of CH3_clash_list_val.txt)
		n_atoms = 28;
		n_chi = 1;
	}
	else if (AA_name.compare("Ile") == 0){
		clash_size = 381;
		CH3_clash_size = 147;
		n_atoms = 31;
		n_chi = 2;
	}

	
	
	// Load clash list file
	int clash_atom1[clash_size], clash_atom2[clash_size]; 		// variables to hold atom indexes
	float radii_sum[clash_size], radii2[clash_size];		// variables to hold sum of radii and (sum of radii)^2
	string clash_file = "Clash_list_val.txt";
	string radii_file = "Clash_list_radii_sum_val.txt";
	load_clash_array(clash_file, clash_size, clash_atom1, clash_atom2);
	load_radii_array(radii_file, clash_size, radii_sum);
	
	//square the radii
	for (i = 0; i < clash_size; i ++){
	  		radii2[i] = radii_sum[i] * radii_sum[i];
	  	}

	// load CH3 clashlist file
	int CH3_clash_atom1[CH3_clash_size], CH3_clash_atom2[CH3_clash_size];		// varables to hold atom indexes
	float CH3_radii2[CH3_clash_size];											// holds (sum of radii)^2
	clash_file = "CH3_clash_list_val.txt";
	radii_file = "CH3_clash_list_radii2_sum_val.txt";
	load_clash_array(clash_file, CH3_clash_size, CH3_clash_atom1, CH3_clash_atom2);
	load_radii_array(radii_file, CH3_clash_size, CH3_radii2); // This radii is already squared

	// Load coordinates
	float Position[n_atoms][3];
	ifstream myfile;
	myfile.open((folder + input_file).c_str());
	if (myfile.is_open())
	  {
	  for (i = 0; i < n_atoms; i ++){
	  		myfile >> Position[i][0] >> Position[i][1] >> Position[i][2];
	  	}
	  	myfile.close();
	  }
	else{
		cout << "problem with position File" << endl;
		exit(EXIT_FAILURE);
	}


	// Get initial dihedral values for each dihedral
	float phi_init, psi_init, chi1_init, chi2_init, CH3_1_init, CH3_2_init, end_CH3_1_init, end_CH3_2_init;
	phi_init = calculate_DA(phi_index, Position);
	psi_init= calculate_DA(psi_index, Position);
	chi1_init = calculate_DA(chi1_index, Position);
	chi2_init = calculate_DA(chi2_index, Position);
	CH3_1_init= calculate_DA(CH3_1_index, Position);
	CH3_2_init= calculate_DA(CH3_2_index, Position);
	end_CH3_1_init= calculate_DA(end_CH3_1_index, Position);
	end_CH3_2_init= calculate_DA(end_CH3_2_index, Position);

	// Delta term is used when rotating the dihedral angles. Since it is a constant we can calculate it once here
	float delta_term_phi, delta_term_psi, delta_term_chi1, delta_term_chi2, delta_term_CH3_1, delta_term_CH3_2, delta_term_end_CH3_1, delta_term_end_CH3_2;
	delta_term_phi = abs(PI * phi_init / 180.0);
	delta_term_psi = abs(PI *  psi_init / 180.0);
	delta_term_chi1 = abs(PI * chi1_init / 180.0);
	delta_term_chi2 = abs(PI * chi2_init / 180.0);
	delta_term_CH3_1 = abs(PI *  CH3_1_init / 180.0);
	delta_term_CH3_2 = abs(PI *  CH3_2_init / 180.0);
	delta_term_end_CH3_1 = abs(PI * end_CH3_1_init / 180.0);
	delta_term_end_CH3_2 = abs(PI * end_CH3_2_init / 180.0);

	// Declare a bunch of variables
	float newPos[n_atoms][3];
	float setPhi, setPsi, setChi1, setChi2, setCH3_1, setCH3_2, E, minE, CH3_E;
	float Pos_phi[n_atoms][3], Pos_psi[n_atoms][3], Pos_chi1[n_atoms][3], Pos_chi2[n_atoms][3], Pos_CH3_1[n_atoms][3], Pos_CH3_2[n_atoms][3];
	float CH3_min[n_atoms][3];
	float save_data[72 * 72][n_chi+3];
	int save_counter = 0;

	//Clear output file
	string n1, n2, n3, f1;
	n1 = "_energy_";
	n2 = "_";
	n3 = ".txt";
	f1 = save_folder + AA_name + n1 + save_name + n2 + coord_num + n3;
	ofstream out_file1(f1.c_str(), ios::out);
	out_file1 << fixed;

	for (int phi_loop = 0; phi_loop < 72; phi_loop++){					// loop over all phi values

		// Rotate to the new phi value
		setPhi = phi_loop * 5.0;
		rotate_DA(Position, setPhi, delta_term_phi, phi_index, moveAtomID_phi, n_atoms, Pos_phi);

		for (int psi_loop = 0; psi_loop < 72; psi_loop++){				// loop over all psi values

			// Rotate to the new psi valu
			setPsi = psi_loop * 5.0;
			rotate_DA(Pos_phi, setPsi, delta_term_psi, psi_index, moveAtomID_psi, n_atoms, Pos_psi);

			for (int chi1_loop =0; chi1_loop < 72; chi1_loop++ ){ 		// loop over all chi1 values

				// Rotate to the new chi1 valu
				setChi1 = chi1_loop * 5.0;
				rotate_DA(Pos_psi, setChi1, delta_term_chi1, chi1_index, moveAtomID_chi1, n_atoms, Pos_chi1);
				if (n_chi >= 2){

					for (int chi2_loop = 0; chi2_loop < 72; chi2_loop++){
					
						setChi2 = chi2_loop * 5.0;
						rotate_DA(Pos_chi1, setChi2, delta_term_chi2, chi2_index, moveAtomID_chi2, n_atoms, Pos_chi2);
						// Get energy
						E = get_energy(Pos_chi2, radii2, clash_atom1, clash_atom2, clash_size);

						// now we have to consider the terminal CH3 groups
						// Get energy due to CH3 atoms
						CH3_E = get_energy(Pos_chi2, CH3_radii2, CH3_clash_atom1, CH3_clash_atom2, CH3_clash_size);
						// If there is energy due to overlap of the CH3 hydrogens, rotate these extra dihedrals, finding the lowest energy configuration
						if (CH3_E > 0){	

							rotate_CH3_group(n_atoms, Pos_chi2, CH3_radii2, CH3_clash_atom1, CH3_clash_atom2, CH3_clash_size, delta_term_CH3_1, delta_term_CH3_2, CH3_1_index, CH3_2_index, moveAtomID_CH3_1, moveAtomID_CH3_2, CH3_min);
							E = get_energy(CH3_min, radii2, clash_atom1, clash_atom2, clash_size);
						}

						// Save data
						save_data[save_counter][0] = setPhi;
						save_data[save_counter][1] = setPsi;
						save_data[save_counter][2] = setChi1;
						save_data[save_counter][3] = setChi2;
						save_data[save_counter][5] = E;
						save_counter = save_counter + 1;
						if (save_counter >= 72 * 72){
							for (i=0; i < save_counter; i ++){
								out_file1 << save_data[i][0] << " " << save_data[i][1] << " " << save_data[i][2] << " " << save_data[i][3] << " " << save_data[i][4] << endl;
							}
							save_counter = 0;
						}
					} // end chi2 loop
				} // end if Chi 2
				else {			// if only chi1
					// Get energy
					E = get_energy(Pos_chi1, radii2, clash_atom1, clash_atom2, clash_size);

					// now we have to consider the terminal CH3 groups
					// Get energy due to CH3 atoms
					CH3_E = get_energy(Pos_chi1, CH3_radii2, CH3_clash_atom1, CH3_clash_atom2, CH3_clash_size);
					// If there is energy due to overlap of the CH3 hydrogens, rotate these extra dihedrals, finding the lowest energy configuration
					if (CH3_E > 0){
						rotate_CH3_group(n_atoms, Pos_chi1, CH3_radii2, CH3_clash_atom1, CH3_clash_atom2, CH3_clash_size, delta_term_CH3_1, delta_term_CH3_2, CH3_1_index, CH3_2_index, moveAtomID_CH3_1, moveAtomID_CH3_2, CH3_min);
						E = get_energy(CH3_min, radii2, clash_atom1, clash_atom2, clash_size);
					}
					// Save data
					save_data[save_counter][0] = setPhi;
					save_data[save_counter][1] = setPsi;
					save_data[save_counter][2] = setChi1;
					save_data[save_counter][3] = E;
					save_counter = save_counter + 1;
					if (save_counter >= 72 * 72){
						for (i=0; i < save_counter; i ++){
							out_file1 << save_data[i][0] << " " << save_data[i][1] << " " << save_data[i][2] << " " << save_data[i][3] << endl;
						}
						save_counter = 0;
					}


				} // end of if only Chi 1
			} // end of Chi1 loop
		} // end psi loop
	} // end phi loop

	out_file1.close();
}