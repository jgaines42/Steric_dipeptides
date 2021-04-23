#ifndef dihedral_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define dihedral_H

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
float calculate_DA(int AtomArray[], float Position[][3]);


/******************************************************************************************************************************
 void rotate_DA

 Rotates the specified dihedral angles to the value in setChi

 Input:
 - Position[][3]: Nx3 array of atomic coordinates
 - setChi		: desired angle (in degrees) 
 - delta_term	: value based on the original angle --- abs(PI * init_angle / 180.0)
 - iChiArray	: Array of 4 indexes corresponding to the 4 atoms of the dihedral angle
 - moveAtomID	: Array of atomic indexes that will move when angle is rotated (first value is the length of the array)
 - n_atoms		: number of atoms in Position
 - TempPosition : Array that will be used to return the rotated coordinates
******************************************************************************************************************************/
void rotate_DA(float Position[][3], float setChi, float delta_term, int iChiArray[], int moveAtomID[], int n_atoms,float TempPosition[][3]);


/******************************************************************************************************************************
 float get_energy

 Gets the repulsive LJ energy due to atomic overlaps between atoms specified in clash_atom1 and clash_atom2

 Input:
 - Position[][3]: Nx3 array of atomic coordinates
 - radii2		: (r1 + r2)^2 for all atom pairs listed in clash_atom1 and clash_atom2
 - clash_atom1	: Index in Position of first atom of each pair
 - clash_atom2	: Index in Position of second atom of each pair
 - clash_size	: Length of clash_atom1 (number of pairs to check)
******************************************************************************************************************************/
float get_energy(float Position[][3], float radii2[], int clash_atom1[], int clash_atom2[], int clash_size);

/******************************************************************************************************************************
 float load_clash_array

 Loads the atomic indexes of all possible clashes in this system

 Input:
 - clash_file		: String with full path to the file to load
 - clash_size 	: Number of rows in the file
 - clash_atom1	: array that will store all indexes of the first atom in each pair
 - clash_atom2	: array that will store all indexes of the second atom in each pair
******************************************************************************************************************************/
void load_clash_array(std::string clash_file, int clash_size, int clash_atom1[], int clash_atom2[]);


/******************************************************************************************************************************
 float load_radii_array

 Loads the atomic radii file for all possible clashes in this system

 Input:
 - radii_file	: String with full path to the file to load
 - clash_size 	: Number of rows in the file
 - radii2		: array that will store the values loaded from the file
******************************************************************************************************************************/
void load_radii_array(std::string radii_file, int clash_size, float radii2[]);


/******************************************************************************************************************************
 float dot_product

 Calculates the dot product of two vectors (assumes vectors are of length 3)

 Input:
 - vector_a	: vector of length 3
 - vector_b : vector of length 3
******************************************************************************************************************************/
float dot_product(float vector_a[], float vector_b[]);

/******************************************************************************************************************************
 float cross_product

 Calculates the cross product of two vectors (assumes vectors are of length 3)

 Input:
 - vector_a	: vector of length 3
 - vector_b : vector of length 3
 - temp 	: vector of length 3 used to return the cross product
******************************************************************************************************************************/
void cross_product(float vector_a[], float vector_b[], float temp[]);

#endif