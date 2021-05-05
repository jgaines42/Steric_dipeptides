/*
AminoAcidClass.h
Jennifer Mortensen

Template class for an Amino acid
*/

#ifndef AminoA
#define AminoA

#include <iostream>
#include <fstream>
#include <ostream>
#include <cstdlib>
#include <math.h>
#include <string>
#include <sstream>
#include <vector>

class AminoAcid
{


public:

	/**********************
	Constructors
	***********************/
	
	//Default constructor
	AminoAcid();
	AminoAcid(std::string Res_name);
	int num_atoms;
	int num_dihedral;
	int moveAtomID_psi[];
	int moveAtomID_chi1[];
	int moveAtomID_CH3_1[];
	int moveAtomID_CH3_2 [];
	int moveAtomID_end_CH3_1[];
	int moveAtomID_end_CH3_2[];

	~AminoAcid(){}	

};
#endif