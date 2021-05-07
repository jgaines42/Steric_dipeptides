#include "dihedral_functions.h"  
#define PI 3.14159265
using namespace std;


/******************************************************************************************************************************
 float calculate_DA

 Calculates the dihedral angle of a set of atoms
 New approach is from https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python

 Input:
 - AtomArray    : Array of 4 atomic indexes for the dihedral angle
 - Position[][3]: Nx3 array of atomic coordinates
******************************************************************************************************************************/
float calculate_DA(int AtomArray[], float Position[][3]){

	float b0[3], b1[3], b2[3], v[3], w[3], x, y;
	float b_mag = 0.0;
	float temp[3];
	float p0[3], p1[3], p2[3], p3[3];

    // Get coordinates of the 4 atoms in the dihedral angles
	for (int i = 0; i<3; i++){
		p0[i] = Position[AtomArray[0]][i];
		p1[i] = Position[AtomArray[1]][i];
		p2[i] = Position[AtomArray[2]][i];
		p3[i] = Position[AtomArray[3]][i];
	}

    for (int i = 0;i < 3;i++){
    	b0[i] = -1.0 * (p1[i] - p0[i]);
    	b1[i] = p2[i] - p1[i];
    	b2[i] = p3[i] - p2[i];
    	b_mag = b_mag + b1[i] * b1[i];
    }
    
    // normalize b1 so that it does not influence magnitude of vector
    // reflections that come next
    b_mag = sqrt(b_mag);
   	b1[0] = b1[0] / b_mag;
   	b1[1] = b1[1] / b_mag;
   	b1[2] = b1[2] / b_mag;

    // vector projections
    // v = projection of b0 onto plane perpendicular to b1
    //   = b0 minus component that aligns with b1
    // w = projection of b2 onto plane perpendicular to b1
    //   = b2 minus component that aligns with b1
    float dot_b0_b1 = dot_product(b0, b1);
    float dot_b2_b1 = dot_product(b2, b1);
    for (int i= 0; i<3; i++){
    	v[i] = b0[i] - dot_b0_b1 * b1[i];
    	w[i] = b2[i] - dot_b2_b1 * b1[i];
    }
    
    // angle between v and w in a plane is the torsion angle
    // v and w may not be normalized but that's fine since tan is y/x
    x = dot_product(v, w);
    cross_product(b1, v, temp);
    y = dot_product(temp, w);

    float DA = atan2(y, x) * 180.0 / PI;

    // other code requires the initial dihedral to be >= 0
    if (DA < 0){
        DA = DA + 360;
    }

    return DA;

}


/******************************************************************************************************************************
 void rotate_DA

 Rotates the specified dihedral angles to the value in setChi

 Input:
 - Position[][3]: Nx3 array of atomic coordinates
 - setChi       : desired angle (in degrees) 
 - delta_term   : value based on the original angle --- abs(PI * init_angle / 180.0)
 - iChiArray    : Array of 4 indexes corresponding to the 4 atoms of the dihedral angle
 - moveAtomID   : Array of atomic indexes that will move when angle is rotated (first value is the length of the array)
 - n_atoms      : number of atoms in Position
 - TempPosition : Array that will be used to return the rotated coordinates
******************************************************************************************************************************/

void rotate_DA(float Position[][3], float setChi, float delta_term, int iChiArray[], int moveAtomID[], int n_atoms, float TempPosition[][3]){

	float q0, q1, q2, q3, q02, q12, q22, q32, sindelta;
	int i, j, k;


	// Calculate how much the dihedral needs to change (in radians)
    float deltaChi1_F153 = delta_term - setChi * PI / 180.0;

    // Get the coordinates of the 2nd atom in the dihedral
    float subtract_atom[3];
    for (j = 0; j < 3; j++){
    	subtract_atom[j] = Position[iChiArray[1]][j];
    }

  
    // Move all atoms so 2nd atom of dihedral is at orgin
    for (i = 0; i< n_atoms; i++){
    	TempPosition[i][0] = Position[i][0] - subtract_atom[0];
    	TempPosition[i][1] = Position[i][1] - subtract_atom[1];
    	TempPosition[i][2] = Position[i][2] - subtract_atom[2];
	}

    // Get the vector from the origin to the 3rd atom of the dihedral.
    // This is the vector that we will rotate around

	float CAtoCB_F153[3];
	float CACB_mag = 0.0;
    for (i = 0; i < 3; i++){
    	CAtoCB_F153[i] = -TempPosition[iChiArray[2]][i];
    	CACB_mag = CACB_mag + CAtoCB_F153[i]*CAtoCB_F153[i];
    }
    CACB_mag = sqrt(CACB_mag);

    CAtoCB_F153[0] = CAtoCB_F153[0] / CACB_mag;
    CAtoCB_F153[1] = CAtoCB_F153[1] / CACB_mag;
    CAtoCB_F153[2] = CAtoCB_F153[2] / CACB_mag;

    // Do complicated math
    q0 = cos(deltaChi1_F153 / 2.0);
    sindelta = sin(deltaChi1_F153 / 2.0);
    q1 = CAtoCB_F153[0] * sindelta;
    q2 = CAtoCB_F153[1] * sindelta;
    q3 = CAtoCB_F153[2] * sindelta;
    q02 = q0 * q0;
    q12 = q1 * q1;
    q22 = q2 * q2;
    q32 = q3 * q3;

    // Q is the rotation vector
    float Q[3][3] = {{(q02 + q12 - q22 - q32), 2.0 * (q1 * q2 - q0 * q3), 2.0 * (q0 * q2 + q1 * q3)},
         {2.0 * (q1 * q2 + q0 * q3), (q02 - q12 + q22 - q32), 2.0 * (-q0 * q1 + q2 * q3)},
         {2.0 * (-q0 * q2 + q1 * q3), 2 * (q0 * q1 + q2 * q3), (q02 - q12 - q22 + q32)}};

    // Extract the coordinates that should move
    int n_move = moveAtomID[0] - 1;
    float move_coordinates[n_move][3];
    for (i = 1; i <= n_move; i++){
	    move_coordinates[i-1][0] = TempPosition[moveAtomID[i]][0];
	    move_coordinates[i-1][1] = TempPosition[moveAtomID[i]][1];
	    move_coordinates[i-1][2] = TempPosition[moveAtomID[i]][2];
	}

    // Initialize newPos
	float newPos[n_move][3];
	for (i = 0; i < n_move; i++){
		for (j = 0; j < 3; j++){
			newPos[i][j] = 0;
		}
	}

	// Use Q to rotate these atoms
    // This code does matrix multiplication
	for(i = 0; i < 3; ++i){
        for(j = 0; j < n_move; j++){
            for(k = 0; k < 3; k++)
            {
                newPos[j][i] += Q[i][k] * move_coordinates[j][k];
            }
        }
	}
	
    // Put the moved atoms back into the original array
	for (i = 1; i <= n_move; i++){
	    TempPosition[moveAtomID[i]][0] = newPos[i-1][0];
	    TempPosition[moveAtomID[i]][1] = newPos[i-1][1];
	    TempPosition[moveAtomID[i]][2] = newPos[i-1][2];
	}

    // Translate atoms back to original location
	for (i = 0; i< n_atoms; i++){
    	for (j = 0; j < 3; j++){
	    	TempPosition[i][j] = TempPosition[i][j] + subtract_atom[j];
	    }
	}
	
}

/******************************************************************************************************************************
 float get_energy

 Gets the repulsive LJ energy due to atomic overlaps between atoms specified in clash_atom1 and clash_atom2

 Input:
 - Position[][3]: Nx3 array of atomic coordinates
 - radii2       : (r1 + r2)^2 for all atom pairs listed in clash_atom1 and clash_atom2
 - clash_atom1  : Index in Position of first atom of each pair
 - clash_atom2  : Index in Position of second atom of each pair
 - clash_size   : Length of clash_atom1 (number of pairs to check)
******************************************************************************************************************************/
float get_energy(float Position[][3], float radii2[], int clash_atom1[], int clash_atom2[], int clash_size){
	int loc1, loc2;
	float E = 0.0;
	float D1, E1, E3;
    float dif1, dif2, dif3; // X, Y, Z compents of vector between the two atoms

    // Loop over all pairs of atoms
	for (int i = 0;i < clash_size;i++){

        // Get each atom id
		loc1 = clash_atom1[i];
		loc2 = clash_atom2[i];

        // Calculate distance^2 between the two atoms
		dif1 = Position[loc1][0] - Position[loc2][0];
		dif2 = Position[loc1][1] - Position[loc2][1];
		dif3 = Position[loc1][2] - Position[loc2][2];
		D1 = dif1*dif1 + dif2*dif2 + dif3*dif3;

        // If distance^2 is less than (r1 + r2)^2, calculate the repulsive LJ energy
		if (D1 < radii2[i]){

			E3 = (radii2[i] / D1) * (radii2[i] / D1) * (radii2[i] / D1);
			E1 = (1 - E3);
            
			E = E + E1*E1;
		}
	}
	return E;
}


/******************************************************************************************************************************
 float load_clash_array

 Loads the atomic indexes of all possible clashes in this system

 Input:
 - clash_file       : String with full path to the file to load
 - clash_size   : Number of rows in the file
 - clash_atom1  : array that will store all indexes of the first atom in each pair
 - clash_atom2  : array that will store all indexes of the second atom in each pair
******************************************************************************************************************************/
void load_clash_array(std::string clash_file, int clash_size, int clash_atom1[], int clash_atom2[]){
    ifstream myfile;
    myfile.open(clash_file.c_str());
    if (myfile.is_open())
      {
      for (int i = 0; i < clash_size; i ++){
            myfile >> clash_atom1[i] >> clash_atom2[i];
        }
        myfile.close();
      }
    else{
        cout << "problem with clash File " << clash_file << endl;
        exit(EXIT_FAILURE);
    }
}


/******************************************************************************************************************************
 float load_radii_array

 Loads the atomic radii file for all possible clashes in this system

 Input:
 - radii_file   : String with full path to the file to load
 - clash_size   : Number of rows in the file
 - radii2       : array that will store the values loaded from the file
******************************************************************************************************************************/
void load_radii_array(std::string radii_file, int clash_size, float radii2[]){
    ifstream myfile;
    myfile.open(radii_file.c_str());
    if (myfile.is_open())
      {
      for (int i = 0; i < clash_size; i ++){
            myfile >> radii2[i];
        }
        myfile.close();
      }
    else{
        cout << "problem with radii File" << endl;
        exit(EXIT_FAILURE);
    }


}


/******************************************************************************************************************************
 float dot_product

 Calculates the dot product of two vectors (assumes vectors are of length 3)
 from https://www.tutorialspoint.com/cplusplus-program-for-dot-product-and-cross-product-of-two-vectors

 Input:
 - vector_a : vector of length 3
 - vector_b : vector of length 3
******************************************************************************************************************************/
float dot_product(float vector_a[], float vector_b[]) {
   float product = 0.0;
   product = product + vector_a[0] * vector_b[0] + vector_a[1] * vector_b[1] + vector_a[2] * vector_b[2];
   return product;
}

/******************************************************************************************************************************
 float cross_product

 Calculates the cross product of two vectors (assumes vectors are of length 3)
 from https://www.tutorialspoint.com/cplusplus-program-for-dot-product-and-cross-product-of-two-vectors

 Input:
 - vector_a : vector of length 3
 - vector_b : vector of length 3
 - temp     : vector of length 3 used to return the cross product
******************************************************************************************************************************/
void cross_product(float vector_a[], float vector_b[], float temp[]) {
   temp[0] = vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1];
   temp[1] = vector_a[2] * vector_b[0] - vector_a[0] * vector_b[2];
   temp[2] =  vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0];
}