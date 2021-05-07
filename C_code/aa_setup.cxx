#include "aa_setup.h"  
#define PI 3.14159265
using namespace std;
// int phi_index[], int psi_index[], int chi1_index[], int chi2_index[], int chi3_index[], int CH3_1_index[], int CH3_2_index[], int end_CH3_1_index[], int end_CH3_2_index[]
void get_dihedral_indexes(std::string aa_name, int index_array[9][4]){
	
	if (aa_name.compare("Val") == 0){
		int this_index[10][4] = {{4, 6, 8, 20}, 		// phi
								{6, 8, 20, 22},		// psi
								{6, 8, 10, 12},		// chi1
								{-1, -1, -1, -1}, 	//chi2
								{-1, -1, -1, -1}, 	//chi3
								{8, 10, 12, 13}, 	//CH3 1
								{8, 10, 16, 17}, 	//CH3 2
								{6, 4, 1, 0},		// end CH3 1
								{20, 22, 24, 25}, 	// end CH3 2
								{-1, -1, -1, -1}, 	//OH
							};
		for (int i = 0; i < 10; i++){
			for (int j = 0; j < 4; j++){
				index_array[i][j]=this_index[i][j];
			}
		}

	}

	else if (aa_name.compare("Ile") == 0){
		int this_index[10][4] = {{4, 6, 8, 23}, 		// phi
								{6, 8, 23, 25},		// psi
								{6, 8, 10, 16},		// chi1
								{8, 10, 16, 19}, 	//chi2
								{-1, -1, -1, -1}, 	//chi3
								{8, 10, 12, 13}, 	//CH3 1
								{10, 16, 19, 20}, 	//CH3 2
								{6, 4, 1, 0},		// end CH3 1
								{23, 25, 27, 28}, 	// end CH3 2
								{-1, -1, -1, -1}, 	//OH
							};
		for (int i = 0; i < 10; i++){
			for (int j = 0; j < 4; j++){
				index_array[i][j]=this_index[i][j];
			}
		}
	}
	else if (aa_name.compare("Ser") == 0){
		int this_index[10][4] = {{4, 6, 8, 15}, 		// phi
								{6, 8, 15, 17},		// psi
								{6, 8, 10, 13},		// chi1
								{-1, -1, -1, -1}, 	//chi2
								{-1, -1, -1, -1}, 	//chi3
								{-1, -1, -1, -1}, 	//CH3 1
								{-1, -1, -1, -1}, 	//CH3 2
								{6, 4, 1, 0},		// end CH3 1
								{15, 17, 19, 21}, 	// end CH3 2
								{8, 10, 13, 14},	// OH
							};
		for (int i = 0; i < 10; i++){
			for (int j = 0; j < 4; j++){
				index_array[i][j]=this_index[i][j];
			}
		}
	}
	else {
		int this_index[10][4] = {{4, 6, 8, 20}, 		// phi
								{6, 8, 20, 22},		// psi
								{6, 8, 10, 12},		// chi1
								{-1, -1, -1, -1}, 	//chi2
								{-1, -1, -1, -1}, 	//chi3
								{8, 10, 12, 13}, 	//CH3 1
								{8, 10, 16, 17}, 	//CH3 2
								{6, 4, 1, 0},		// end CH3 1
								{20, 22, 24, 25}, 	// end CH3 2
								{-1, -1, -1, -1}, 	//OH
							};
		for (int i = 0; i < 10; i++){
			for (int j = 0; j < 4; j++){
				index_array[i][j]=this_index[i][j];
			}
		}
	}

}


void get_move_indexes(std::string aa_name, int moveAtomID_phi[], int moveAtomID_psi[], int moveAtomID_chi1[], int moveAtomID_chi2[], int moveAtomID_chi3[], int moveAtomID_CH3_1[], int moveAtomID_CH3_2[], int moveAtomID_end_CH3_1[], int moveAtomID_end_CH3_2[], int moveAtomID_OH[]){
		
	if (aa_name.compare("Val") == 0){
		int ID_phi[30] = {20, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27};
		int ID_psi [30]= {8, 21, 22, 23, 24, 25, 26, 27};
		int ID_chi1[30] = {10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
		int ID_CH3_1[30] = {4, 13, 14, 15};
		int ID_CH3_2 [30]= {4, 17, 18, 19};
		int ID_end_CH3_1[30] = {4, 0, 2, 3};
		int ID_end_CH3_2[30] = {4, 25, 26, 27};
		std::copy(ID_phi, ID_phi+ID_phi[0], moveAtomID_phi);
		std::copy(ID_psi, ID_psi+ID_psi[0], moveAtomID_psi);
		std::copy(ID_chi1, ID_chi1+ID_chi1[0], moveAtomID_chi1);
		std::copy(ID_CH3_1, ID_CH3_1+ID_CH3_1[0], moveAtomID_CH3_1);
		std::copy(ID_CH3_2, ID_CH3_2+ID_CH3_2[0], moveAtomID_CH3_2);
		std::copy(ID_end_CH3_1, ID_end_CH3_1+ID_end_CH3_1[0], moveAtomID_end_CH3_1);
		std::copy(ID_end_CH3_2, ID_end_CH3_2+ID_end_CH3_2[0], moveAtomID_end_CH3_2);
	}
	else if (aa_name.compare("Ile") == 0){
		int ID_phi[30] = {23, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};
		int ID_psi [30]= {8, 24, 25, 26, 27, 28, 29, 30};
		int ID_chi1[30] = {13, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22};
		int ID_chi2[30] = {7, 17, 18, 19, 20, 21, 22};
		int ID_CH3_1[30] = {4, 13, 14, 15};
		int ID_CH3_2 [30]= {4, 20, 21, 22};
		int ID_end_CH3_1[30] = {4, 0, 2, 3};
		int ID_end_CH3_2[30] = {4, 28, 29, 30};
		std::copy(ID_phi, ID_phi+ID_phi[0], moveAtomID_phi);
		std::copy(ID_psi, ID_psi+ID_psi[0], moveAtomID_psi);
		std::copy(ID_chi1, ID_chi1+ID_chi1[0], moveAtomID_chi1);
		std::copy(ID_CH3_1, ID_CH3_1+ID_CH3_1[0], moveAtomID_CH3_1);
		std::copy(ID_CH3_2, ID_CH3_2+ID_CH3_2[0], moveAtomID_CH3_2);
		std::copy(ID_end_CH3_1, ID_end_CH3_1+ID_end_CH3_1[0], moveAtomID_end_CH3_1);
		std::copy(ID_end_CH3_2, ID_end_CH3_2+ID_end_CH3_2[0], moveAtomID_end_CH3_2);
	}
	else if (aa_name.compare("Ser") == 0){
		int ID_phi[30] = {15, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22};
		int ID_psi [30]= {8, 16, 17, 18, 19, 20, 21, 22};
		int ID_chi1[30] = {5, 11, 12, 13, 14};
		int ID_end_CH3_1[30] = {4, 0, 2, 3};
		int ID_end_CH3_2[30] = {4, 28, 29, 30};
		int ID_end_OH[30] = {2, 14};
		std::copy(ID_phi, ID_phi+ID_phi[0], moveAtomID_phi);
		std::copy(ID_psi, ID_psi+ID_psi[0], moveAtomID_psi);
		std::copy(ID_chi1, ID_chi1+ID_chi1[0], moveAtomID_chi1);
		std::copy(ID_end_OH, ID_end_OH+ID_end_OH[0], moveAtomID_OH);
		std::copy(ID_end_CH3_1, ID_end_CH3_1+ID_end_CH3_1[0], moveAtomID_end_CH3_1);
		std::copy(ID_end_CH3_2, ID_end_CH3_2+ID_end_CH3_2[0], moveAtomID_end_CH3_2);

	}


}

