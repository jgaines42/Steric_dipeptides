#include "rotate_OH_group.h"  
#include "dihedral_functions.h" 

using namespace std;

void rotate_OH_group(int n_atoms, float Position[][3], float OH_radii2[], int OH_clash_atom1[], int OH_clash_atom2[], int OH_clash_size, float delta_term_OH_1, int OH_1_index[], int moveAtomID_OH[], float OH_min[][3]){
	float Pos_OH_1[n_atoms][3];
	int i, j, k;
	// Store the current coordinates in CH3_min
	for (i = 0; i < n_atoms; i++){
		OH_min[i][0] = Position[i][0];
		OH_min[i][1] = Position[i][1];
		OH_min[i][2] = Position[i][2];
	}
	float OH_E = get_energy(Position, OH_radii2, OH_clash_atom1, OH_clash_atom2, OH_clash_size);
	float minE = OH_E;	// Lowest CH3 energy found so far
	float setOH;

	for (int OH_1_loop = 0; OH_1_loop < 72; OH_1_loop++){		// Rotate first CH3 group
		if (minE == 0){
			break;
		}
		
		setOH = OH_1_loop * 5.0;
		rotate_DA(Position, setOH, delta_term_OH_1, OH_1_index, moveAtomID_OH, n_atoms, Pos_OH_1);
			
		// Calculate energy due to CH3 hydrogen atoms
		OH_E = get_energy(Pos_OH_1, OH_radii2, OH_clash_atom1, OH_clash_atom2, OH_clash_size);

		// If we found a new minimum, store the energy and coordiantes
		if (OH_E < minE){
			minE = OH_E;
			// copy array to CH3_min
			for (i = 0; i < n_atoms; i++){
				OH_min[i][0] = Pos_OH_1[i][0];
				OH_min[i][1] = Pos_OH_1[i][1];
				OH_min[i][2] = Pos_OH_1[i][2];
			}
			if (minE == 0){ 	// If the CH3 energy is 0, exit the loops
				OH_1_loop = 73;
				break;
			}
		} // end new min

	} // end for all CH3 1

}