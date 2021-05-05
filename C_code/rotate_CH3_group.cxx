#include "rotate_CH3_group.h"  
#include "dihedral_functions.h" 

using namespace std;

void rotate_CH3_group(int n_atoms, float Position[][3], float CH3_radii2[], int CH3_clash_atom1[], int CH3_clash_atom2[], int CH3_clash_size, float delta_term_CH3_1, float delta_term_CH3_2, int CH3_1_index[], int CH3_2_index[], int moveAtomID_CH3_1[], int moveAtomID_CH3_2[], float CH3_min[][3]){
	float Pos_CH3_1[n_atoms][3];
	float Pos_CH3_2[n_atoms][3];
	int i, j, k;
	// Store the current coordinates in CH3_min
	for (i = 0; i < n_atoms; i++){
		CH3_min[i][0] = Position[i][0];
		CH3_min[i][1] = Position[i][1];
		CH3_min[i][2] = Position[i][2];
	}
	float CH3_E = get_energy(Position, CH3_radii2, CH3_clash_atom1, CH3_clash_atom2, CH3_clash_size);
	float minE = CH3_E;	// Lowest CH3 energy found so far
	float setCH3_1, setCH3_2;

	for (int CH3_1_loop = 0; CH3_1_loop < 72; CH3_1_loop++){		// Rotate first CH3 group
		if (minE == 0){
			break;
		}
		
		setCH3_1 = CH3_1_loop * 5.0;
		rotate_DA(Position, setCH3_1, delta_term_CH3_1, CH3_1_index, moveAtomID_CH3_1, n_atoms, Pos_CH3_1);

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
			} // end new min
		} // end for all CH3 2

	} // end for all CH3 1

}