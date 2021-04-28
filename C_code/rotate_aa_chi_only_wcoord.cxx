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

	int n_atoms = 28;

	// dihedral index variables for Val
	int phi_index[4] = {4, 6, 8, 20};
	int psi_index[4] = {6, 8, 20, 22};
	int chi1_index[4] = {6, 8, 10, 12};
	int CH3_1_index[4] = {8, 10, 12, 13};
	int CH3_2_index[4] = {8, 10, 16, 17};
	int end_CH3_1_index[4] = {6, 4, 1, 0};
	int end_CH3_2_index[4] = {20, 22, 24, 25};

	// Atom indexes of atoms that move for each dihedral in Val (first item is the size of the array)
	int moveAtomID_phi[20] = {20, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27};
	int moveAtomID_psi [9]= {9, 21, 22, 23, 24, 25, 26, 27, 28};
	int moveAtomID_chi1[10] = {10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
	int moveAtomID_CH3_1[4] = {4, 13, 14, 15};
	int moveAtomID_CH3_2 [4]= {4, 17, 18, 19};
	int moveAtomID_end_CH3_1[4] = {4, 0, 2, 3};
	int moveAtomID_end_CH3_2[4] = {4, 25, 26, 27};

	int clash_size = 303;		// number of clashes (length of Clash_list_val.txt)
	int CH3_clash_size = 129;	// number of clashes involving terminal CH3 atoms (length of CH3_clash_list_val.txt)
	
	// Load clash list file
	int clash_atom1[clash_size], clash_atom2[clash_size]; 		// variables to hold atom indexes
	float radii_sum[clash_size], radii2[clash_size];		// variables to hold sum of radii and (sum of radii)^2
	string clash_file = "Clash_list_Val.txt";
	string radii_file = "Clash_list_radii_sum_Val.txt";
	load_clash_array(clash_file, clash_size, clash_atom1, clash_atom2);
	load_radii_array(radii_file, clash_size, radii_sum);
	
	//square the radii
	for (i = 0; i < clash_size; i ++){
	  		radii2[i] = radii_sum[i] * radii_sum[i];
	  	}

	// load CH3 clashlist file
	int CH3_clash_atom1[CH3_clash_size], CH3_clash_atom2[CH3_clash_size];		// varables to hold atom indexes
	float CH3_radii[CH3_clash_size], CH3_radii2[CH3_clash_size];				// variables to hold sum of radii and (sum of radii)^2
	clash_file = "CH3_clash_list_Val.txt";
	radii_file = "CH3_clash_list_radii_sum_Val.txt";
	load_clash_array(clash_file, CH3_clash_size, CH3_clash_atom1, CH3_clash_atom2);
	load_radii_array(radii_file, CH3_clash_size, CH3_radii);
	//square the radii
	for (i = 0; i < clash_size; i ++){
	  		CH3_radii2[i] = CH3_radii[i] * CH3_radii[i];
	  	}


	// Load coordinates
	float all_Position[28*1000][3];
	ifstream myfile;
	myfile.open((folder + input_file).c_str());
	if (myfile.is_open())
	  {
	  for (i = 0; i < n_atoms*1000; i ++){
	  		myfile >> all_Position[i][0] >> all_Position[i][1] >> all_Position[i][2];
	  	}
	  	myfile.close();
	  }
	else{
		cout << "problem with position File" << endl;
		exit(EXIT_FAILURE);
	}

	float phi_init, psi_init, chi1_init, CH3_1_init, CH3_2_init, end_CH3_1_init, end_CH3_2_init;
	float delta_term_phi, delta_term_psi, delta_term_chi1, delta_term_CH3_1, delta_term_CH3_2, delta_term_end_CH3_1, delta_term_end_CH3_2;
		// Declare a bunch of variables
	float newPos[n_atoms][3];
	float setPhi, setPsi, setChi1, setCH3_1, setCH3_2, E, minE, CH3_E;
	float Pos_phi[n_atoms][3], Pos_psi[n_atoms][3], Pos_chi1[n_atoms][3], Pos_CH3_1[n_atoms][3], Pos_CH3_2[n_atoms][3];
	float CH3_min[n_atoms][3];
	float Position[n_atoms][3];
	float save_data[7200][4];
	int save_counter = 0;

	//Clear output file
	string n1, n2, n3, f1;
	n1 = "Val_energy_1000_";
	n2 = "_";
	n3 = ".txt";
	f1 = save_folder + n1 + save_name + n2 + coord_num + n3;
	ofstream out_file1(f1.c_str(), ios::out);
	out_file1 << fixed;


	// Loop over all 1000 sets of initial coordinates
	for (int coord_loop = 0; coord_loop < 1000; coord_loop++){
		for (i = 0; i < n_atoms; i++){
			Position[i][0] = all_Position[i + coord_loop*n_atoms][0];
			Position[i][1] = all_Position[i + coord_loop*n_atoms][1];
			Position[i][2] = all_Position[i + coord_loop*n_atoms][2];
		}

		// Get initial dihedral values for each dihedral
		phi_init = calculate_DA(phi_index, Position);
		psi_init= calculate_DA(psi_index, Position);
		chi1_init = calculate_DA(chi1_index, Position);
		CH3_1_init= calculate_DA(CH3_1_index, Position);
		CH3_2_init= calculate_DA(CH3_2_index, Position);
		end_CH3_1_init= calculate_DA(end_CH3_1_index, Position);
		end_CH3_2_init= calculate_DA(end_CH3_2_index, Position);

		// Delta term is used when rotating the dihedral angles. Since it is a constant we can calculate it once here
		delta_term_phi = abs(PI * phi_init / 180.0);
		delta_term_psi = abs(PI *  psi_init / 180.0);
		delta_term_chi1 = abs(PI * chi1_init / 180.0);
		delta_term_CH3_1 = abs(PI *  CH3_1_init / 180.0);
		delta_term_CH3_2 = abs(PI *  CH3_2_init / 180.0);
		delta_term_end_CH3_1 = abs(PI * end_CH3_1_init / 180.0);
		delta_term_end_CH3_2 = abs(PI * end_CH3_2_init / 180.0);



	

		setPhi = phi_init;
		setPsi = psi_init;


		for (int chi1_loop =0; chi1_loop < 72; chi1_loop++ ){ 		// loop over all chi1 values
			

			// Rotate to the new chi1 valu
			setChi1 = chi1_loop * 5.0;
			rotate_DA(Position, setChi1, delta_term_chi1, chi1_index, moveAtomID_chi1, n_atoms, Pos_chi1);

			// Get energy
			E = get_energy(Pos_chi1, radii2, clash_atom1, clash_atom2, clash_size);

			// now we have to consider the terminal CH3 groups
			// Get energy due to CH3 atoms
			CH3_E = get_energy(Pos_chi1, CH3_radii2, CH3_clash_atom1, CH3_clash_atom2, CH3_clash_size);
			
			// If there is energy due to overlap of the CH3 hydrogens, rotate these extra dihedrals, finding the lowest energy configuration
			if (CH3_E > 0){

				// Store the current coordinates in CH3_min
				for (i = 0; i < n_atoms; i++){
					CH3_min[i][0] = Pos_chi1[i][0];
					CH3_min[i][1] = Pos_chi1[i][1];
					CH3_min[i][2] = Pos_chi1[i][2];
				}

				minE = CH3_E;	// Lowest CH3 energy found so far

				for (int CH3_1_loop = 0; CH3_1_loop < 72; CH3_1_loop++){		// Rotate first CH3 group
					
					if (minE == 0){
						break;
					}
					
					setCH3_1 = CH3_1_loop * 5.0;
					rotate_DA(Pos_chi1, setCH3_1, delta_term_CH3_1, CH3_1_index, moveAtomID_CH3_1, n_atoms, Pos_CH3_1);

					for (int CH3_2_loop = 0; CH3_2_loop < 72; CH3_2_loop++){	// Rotate second CH3 group
						setCH3_2 = CH3_2_loop * 5.0;
						rotate_DA(Pos_CH3_1, setCH3_2, delta_term_CH3_2, CH3_2_index, moveAtomID_CH3_2, n_atoms, Pos_CH3_2);
						
						// Calculate energy due to CH3 hydrogen atoms
						CH3_E = get_energy(Pos_CH3_2, CH3_radii2, CH3_clash_atom1, CH3_clash_atom2, CH3_clash_size);

						// If we found a new minimum, store the energy and coordiantes
						if (CH3_E < minE){
							minE = CH3_E;
							// copy array to CH3_min
							for (i = 0; i < n_atoms; i++){
								CH3_min[i][0] = Pos_CH3_2[i][0];
								CH3_min[i][1] = Pos_CH3_2[i][1];
								CH3_min[i][2] = Pos_CH3_2[i][2];
							}
							if (minE == 0){ 	// If the CH3 energy is 0, exit the loops
								CH3_1_loop = 73;
								CH3_2_loop = 73;
								break;
							}
						}
					}

				}
				// Now that we've rotated and found the lowest energy CH3 configuration, calculate the full energy of the molecule
				E = get_energy(CH3_min, radii2, clash_atom1, clash_atom2, clash_size);
			}
			
			// Save data
			save_data[save_counter][0] = setPhi;
			save_data[save_counter][1] = setPsi;
			save_data[save_counter][2] = setChi1;
			save_data[save_counter][3] = E;
			save_counter = save_counter + 1;
			if (save_counter >= 7200){
				for (i=0; i < save_counter; i ++){
					out_file1 << save_data[i][0] << " " << save_data[i][1] << " " << save_data[i][2] << " " << save_data[i][3] << endl;
				}
				save_counter = 0;
			}
			
		} // end chi 1 loop
	} // end coord loop

	out_file1.close();
}