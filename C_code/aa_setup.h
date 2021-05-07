#ifndef aa_sethup_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define aa_sethup_H

#include <vector>
#include <iostream>
#include <fstream>
#include <ostream>
#include <cstdlib>
#include <cmath>
#include <string.h>
#include <stdlib.h> 

/******************************************************************************************************************************
 float calculate_DA

 Calculates the dihedral angle of a set of atoms

 Input:
 - AtomArray	: Array of 4 atomic indexes for the dihedral angle
 - Position[][3]: Nx3 array of atomic coordinates
******************************************************************************************************************************/
void get_dihedral_indexes(std::string aa_name, int index_array[9][4]); //int* return_phi, int return_psi[], int return_c1[], int return_c2[], int return_c3[], int return_CH31[], int return_CH32[], int return_eCH31[], int return_eCH32[]);

//void get_dihedral_indexes(std::string aa_name, int index_array[9][4]);
//int phi_index[], int psi_index[], int chi1_index[], int chi2_index[], int chi3_index[], int CH3_1_index[], int CH3_2_index[], int end_CH3_1_index[], int end_CH3_2_index[]
void get_move_indexes(std::string aa_name, int moveAtomID_phi[], int moveAtomID_psi[], int moveAtomID_chi1[], int moveAtomID_chi2[], int moveAtomID_chi3[], int moveAtomID_CH3_1[], int moveAtomID_CH3_2[], int moveAtomID_end_CH3_1[], int moveAtomID_end_CH3_2[], int moveAtomID_OH[]);

#endif