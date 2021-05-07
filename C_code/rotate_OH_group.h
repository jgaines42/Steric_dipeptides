#ifndef OH_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define OH_H


#include <stdlib.h> 

void rotate_OH_group(int n_atoms, float Position[][3], float OH_radii2[], int OH_clash_atom1[], int OH_clash_atom2[], int OH_clash_size, float delta_term_OH_1, int OH_1_index[], int moveAtomID_OH[], float CH3_min[][3]);
#endif