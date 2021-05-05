#ifndef CH3_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define CH3_H


#include <stdlib.h> 

void rotate_CH3_group(int n_atoms, float Position[][3], float CH3_radii2[], int CH3_clash_atom1[], int CH3_clash_atom2[], int CH3_clash_size, float delta_term_CH3_1, float delta_term_CH3_2, int CH3_1_index[], int CH3_2_index[], int moveAtomID_CH3_1[], int moveAtomID_CH3_2[], float CH3_min[][3]);
#endif