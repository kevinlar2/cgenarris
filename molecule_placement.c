#include <stdio.h>
#include <math.h>
#include <string.h>
#include "spglib.h"
#include "algebra.h"
#include "read_input.h"
#include "spg_generation.h"
#include "molecule_placement.h"
#include "molecule_utils.h"
#include "crystal_utils.h"
#include "randomgen.h"
#include "combinatorics.h"


#define PI 3.141592653

extern float viewing_directions[16][3];
#pragma omp threadprivate(viewing_directions)
extern const int num_viewing_direction;

/*
//temp database
DATABASE database[530]; 



void load_database()
{
//hall number 60 is [59]	
database[59].I.number_of_positions = 4;

database[59].I.wyckoff_letter_number[0] = 1;  
database[59].I.wyckoff_letter_number[1] = 2; 
database[59].I.wyckoff_letter_number[2] = 3; 
database[59].I.wyckoff_letter_number[3] = 4; 

database[59].I.positions[0][0] = 0.5;
database[59].I.positions[0][1] = 0;
database[59].I.positions[0][2] = 0.5;

database[59].I.positions[1][0] = 0;
database[59].I.positions[1][1] = 0;
database[59].I.positions[1][2] = 0.5;

database[59].I.positions[2][0] = 0.5;
database[59].I.positions[2][1] = 0;
database[59].I.positions[2][2] = 0;

database[59].I.positions[3][0] = 0;
database[59].I.positions[3][1] = 0;
database[59].I.positions[3][2] = 0;

//mirror plane
//number of mirror planes for spg
//hallnumber 400 is [399]
database[399].M.number_of_positions = 5;
//number corresponding to wyckoff letter
database[399].M.wyckoff_letter_number[0] = 20;
database[399].M.normal[0] = 2 ;

//cordinates od the first position of first mirror site
database[399].M.first_position_trans[0][0] = 0;
database[399].M.first_position_trans[0][1] = 0.5;
database[399].M.first_position_trans[0][2] = 0;

database[399].M.first_position_rot[0][0][0] = 1;
database[399].M.first_position_rot[0][0][1] = 0;
database[399].M.first_position_rot[0][0][2] = 0;
database[399].M.first_position_rot[0][1][0] = 0;
database[399].M.first_position_rot[0][1][1] = 0;
database[399].M.first_position_rot[0][1][2] = 0;
database[399].M.first_position_rot[0][2][0] = 0;
database[399].M.first_position_rot[0][2][1] = 0;
database[399].M.first_position_rot[0][2][2] = 1;

}





//has some serous bug. molecule is deformed
void place_molecule_at_general_position(crystal *Xtal, molecule *mol,
	int hall_number)
{
	//declare variables
	//num of atoms in the molecule is N
	int N = mol->num_of_atoms;
	int Z = Xtal->Z;
	int num_atoms_in_cell = N*Z;
	
	//generate random angles for rotation
	//phi, the azimuthal angle
	float phi = 2*PI*uniform_dist_01();
	//theta, angle with x-axis. this should acos(u), -1<u<1
	float u = 2*uniform_dist_01() - 1;
	float theta = acos(u);
	//choose axis
	float axis[]={cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)};
	//random rotation psi around axis 
	float psi = 2*PI*uniform_dist_01();
	
	//create rotation matix around the specified axis with angle psi
	float random_rotation_matrix[3][3];
	rotation_mat_around_axis(random_rotation_matrix, axis, psi);
	//randomly rotate
	molecule_rotate(mol, random_rotation_matrix);
	
	//get symmetry operations from spglib database
	double translations[192][3];
	int rotations[192][3][3];	
	int num_of_operations = spg_get_symmetry_from_database(rotations,
	translations, hall_number);
	Z = num_of_operations;
		
	//for debugging with spglib
	double positions[num_atoms_in_cell][3];
	int types[num_atoms_in_cell];

	//calculate fractional coordinates
	float inverse_lattice_vectors[3][3];
	float *mol_Xfrac = (float *)malloc(N*sizeof(float));
	float *mol_Yfrac = (float *)malloc(N*sizeof(float));
	float *mol_Zfrac = (float *)malloc(N*sizeof(float));
	inverse_mat3b3(inverse_lattice_vectors, Xtal->lattice_vectors);
	for(int i = 0; i < N; i++)
	{
		float temp[3] = {(*mol).X[i], (*mol).Y[i], (*mol).Z[i]};
		vector3_mat3b3_multiply(inverse_lattice_vectors, temp, temp);
		mol_Xfrac[i] = temp[0];
		mol_Yfrac[i] = temp[1];
		mol_Zfrac[i] = temp[2];		
	}
	
	//add the random shift to place the molecule in asym cell
	float rand_fracx = uniform_dist_01();
	float rand_fracy = uniform_dist_01();
	float rand_fracz = uniform_dist_01();
	for(int i = 0; i < N; i++)
	{
		mol_Xfrac[i] += rand_fracx;
		mol_Yfrac[i] += rand_fracy;
		mol_Zfrac[i] += rand_fracz; 
	}
	
	//debug
	//float P[3][3]={{ 2./3,-1./3,-1./3 }, { 1./3, 1./3,-2./3 }, { 1./3, 1./3, 1./3 }};
	//float P[3][3] = {1,0,1,-1,1,1,0,-1,1};
	float P[3][3] = {1,0,0,0,1,0,0,0,1};
	float Pinv[3][3];
	inverse_mat3b3(Pinv, P);
	//print_mat3b3(P);
	//print_mat3b3(Pinv);
	
	float temp_mat[3][3];
	float lattice_vectors_transpose[3][3];
	copy_mat3b3_mat3b3(lattice_vectors_transpose, Xtal->lattice_vectors);
	mat3b3_transpose(lattice_vectors_transpose, lattice_vectors_transpose);
	
	//i loops over all the operations
	for(int i = 0; i < num_of_operations; i++)
	{
		int molecule_index = i*N; 
		//deosnt work float rot[3][3] = rotations[i][3][3];
		int rot[3][3];
		copy_intmat3b3_intmat3b3bN(rot, rotations, i);
	
		//printf("rotations = %d \n",rot[0][0]);
		
		//copy ith translations
		float trans[3];
		//copy_vector3_vector3bN(trans, translations, i);
		//printf("trans= %f %f %f\n", translations[i][0], translations[i][1], translations[i][2]);
		trans[0] = translations[i][0];
		trans[1] = translations[i][1];
		trans[2] = translations[i][2];
		//print_vec3(trans);
		
		//debug
		//rot=pinv rot P
		copy_floatmat3b3_intmat3b3(temp_mat, rot);
		mat3b3_mat3b3_multiply(temp_mat, P, temp_mat);
		mat3b3_mat3b3_multiply(Pinv, temp_mat, temp_mat);
		//print_mat3b3(temp_mat);
		//trans Pinv trans
		vector3_mat3b3_multiply(Pinv, trans, trans);
		
		
		//printf("i index = %d \n", i);
		//j loops over symmetry operations
		for(int j = 0; j < N; j++)
		{
			float atomj_array[3];
			atomj_array[0] = mol_Xfrac[j];
			atomj_array[1] = mol_Yfrac[j];
			atomj_array[2] = mol_Zfrac[j]; 
			//print_vec3(atomj_array);
			//print_mat3b3(rot);
			//printf("ist element = %d ", rot[0][1]);
			//vector3_intmat3b3_multiply(rot, atomj_array, atomj_array);
			vector3_mat3b3_multiply(temp_mat, atomj_array, atomj_array);
			vector3_add(trans, atomj_array, atomj_array);
			//print_vec3(atomj_array);
			
			//temp_array has the fractional cord of jth atom in ith mol
			//convert to cartesian
			vector3_mat3b3_multiply(lattice_vectors_transpose,
				atomj_array, atomj_array);
			//print_vec3(atomj_array);	
			//copy to the structure
			Xtal->Xcord[molecule_index+j] = atomj_array[0];
			Xtal->Ycord[molecule_index+j] = atomj_array[1];
			Xtal->Zcord[molecule_index+j] = atomj_array[2];
			Xtal->atoms[molecule_index+j] = (*mol).atoms[j];
			
			//types[molecule_index+j] = 1;
			
			//printf("%f %f %f  \n", random_crystal.Xcord[molecule_index+j], random_crystal.Ycord[molecule_index+j], random_crystal.Zcord[molecule_index+j]);
			
			//printf("j index = %d \n", j);
			//printf("%d molecule %d atom = %f %f %f \n", i+1, j+1, mol->X[j], mol->Y[j], mol->Z[j] );
			//print_vec3(atomj_array);
		}
	}
	
	Xtal->Z = Z;

	free(mol_Xfrac);
	free(mol_Yfrac);
	free(mol_Zfrac);
		
}
*/
/*
void place_molecule_at_inversion_center(crystal *Xtal, molecule *mol,
	int hall_number)
{
	load_database();
	
	//check if spg supports inversion center
	if(database[hall_number-1].I.number_of_positions == 0)
	{
		printf("**inversion center not supported by hall number %d**\n",
			hall_number);
		return ;
	}

	//declare variables
	//num of atoms in the molecule is N
	int N = mol->num_of_atoms;
	int Z = Xtal->Z;
	int num_atoms_in_cell = N*Z;
	
	//generate random angles for rotation
	//phi, the azimuthal angle
	float phi = 2*PI*((float)rand()/(float)(RAND_MAX));
	//theta, angle with x-axis. this should acos(u), -1<u<1
	float u = 2*((float)rand()/(float)(RAND_MAX)) - 1;
	float theta = acos(u);
	//choose axis
	float axis[]={cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)};
	//random rotation psi around axis 
	float psi = 2*PI*((float)rand()/(float)(RAND_MAX));
	
	//create rotation matix around the specified axis with angle psi
	float random_rotation_matrix[3][3];
	rotation_mat_around_axis(random_rotation_matrix, axis, psi);
	//randomly rotate
	molecule_rotate(mol, random_rotation_matrix);
	
	//get symmetry operations from spglib database
	double translations[192][3];
	int rotations[192][3][3];	
	int num_of_operations = spg_get_symmetry_from_database(rotations,
	translations, hall_number);
	Z = num_of_operations;

	//calculate fractional coordinates
	float inverse_lattice_vectors[3][3];
	float *mol_Xfrac = (float *)malloc(N*sizeof(float));
	float *mol_Yfrac = (float *)malloc(N*sizeof(float));
	float *mol_Zfrac = (float *)malloc(N*sizeof(float));
	inverse_mat3b3(inverse_lattice_vectors, Xtal->lattice_vectors);
	//mat3b3_transpose(inverse_lattice_vectors, inverse_lattice_vectors);
	//print_mat3b3(inverse_lattice_vectors);
	for(int i = 0; i < N; i++)
	{
		float temp[3] = {(*mol).X[i], (*mol).Y[i], (*mol).Z[i]};
		vector3_mat3b3_multiply(inverse_lattice_vectors, temp, temp);
		mol_Xfrac[i] = temp[0];
		mol_Yfrac[i] = temp[1];
		mol_Zfrac[i] = temp[2];		
	}
	
	//debug
	//float P[3][3]={{ 2./3,-1./3,-1./3 }, { 1./3, 1./3,-2./3 }, { 1./3, 1./3, 1./3 }};
	//float P[3][3] = {1,0,1,-1,1,1,0,-1,1};
	float P[3][3] = {1,0,0,0,1,0,0,0,1};
	float Pinv[3][3];
	inverse_mat3b3(Pinv, P);
	//print_mat3b3(P);
	//print_mat3b3(Pinv);
	
	//select random inversion center and place molecule
	//int random_select = rand() % database[hall_number-1].I.number_of_positions;
	int random_select = 3;
	//printf("mol_frac = %f %f %f \n", mol_Xfrac[0], mol_Yfrac[0], mol_Zfrac[0]);
	for(int i = 0; i < N; i++)
	{
		mol_Xfrac[i] += database[hall_number-1].I.positions[random_select][0];
		mol_Yfrac[i] += database[hall_number-1].I.positions[random_select][1];
		mol_Zfrac[i] += database[hall_number-1].I.positions[random_select][2];
		
		//printf("mol_frac = %f %f %f \n", mol_Xfrac[i], mol_Yfrac[i], mol_Zfrac[i]);
	}
	

	float temp_mat[3][3];
	float lattice_vectors_transpose[3][3];
	copy_mat3b3_mat3b3(lattice_vectors_transpose, Xtal->lattice_vectors);
	mat3b3_transpose(lattice_vectors_transpose, lattice_vectors_transpose);
	//i loops over all the operations
	for(int i = 0; i < num_of_operations; i++)
	{
		int molecule_index = i*N; 
		//deosnt work float rot[3][3] = rotations[i][3][3];
		int rot[3][3];
		copy_intmat3b3_intmat3b3bN(rot, rotations, i);
	
		//printf("rotations = %d \n",rot[0][0]);
		
		//copy ith translations
		float trans[3];
		//copy_vector3_vector3bN(trans, translations, i);
		//printf("trans= %f %f %f\n", translations[i][0], translations[i][1], translations[i][2]);
		trans[0] = translations[i][0];
		trans[1] = translations[i][1];
		trans[2] = translations[i][2];
		//print_vec3(trans);
		
		//debug
		//rot=pinv rot P
		copy_floatmat3b3_intmat3b3(temp_mat, rot);
		mat3b3_mat3b3_multiply(temp_mat, P, temp_mat);
		mat3b3_mat3b3_multiply(Pinv, temp_mat, temp_mat);
		//print_mat3b3(temp_mat);
		//trans Pinv trans
		vector3_mat3b3_multiply(Pinv, trans, trans);
		
		
		//printf("i index = %d \n", i);
		//j loops over all atoms in a molecule
		for(int j = 0; j < N; j++)
		{
			float atomj_array[3];
			atomj_array[0] = mol_Xfrac[j];
			atomj_array[1] = mol_Yfrac[j];
			atomj_array[2] = mol_Zfrac[j]; 
			//print_vec3(atomj_array);
			//print_mat3b3(rot);
			//printf("ist element = %d ", rot[0][1]);
			//vector3_intmat3b3_multiply(rot, atomj_array, atomj_array);
			vector3_mat3b3_multiply(temp_mat, atomj_array, atomj_array);
			vector3_add(trans, atomj_array, atomj_array);
			//print_vec3(atomj_array);
			
			//temp_array has the fractional cord of jth atom in ith mol
			//convert to cartesian
			vector3_mat3b3_multiply(lattice_vectors_transpose,
				atomj_array, atomj_array);
			//print_vec3(atomj_array);	
			//copy to the structure
			Xtal->Xcord[molecule_index+j] = atomj_array[0];
			Xtal->Ycord[molecule_index+j] = atomj_array[1];
			Xtal->Zcord[molecule_index+j] = atomj_array[2];
			Xtal->atoms[molecule_index+j] = (*mol).atoms[j];
			
			//types[molecule_index+j] = 1;
			
			//printf("%f %f %f  \n", Xtal->Xcord[molecule_index+j], Xtal->Ycord[molecule_index+j], Xtal->Zcord[molecule_index+j]);
			
			//printf("j index = %d \n", j);
			//printf("%d molecule %d atom = %f %f %f \n", i+1, j+1, mol->X[j], mol->Y[j], mol->Z[j] );
			//print_vec3(atomj_array);
		}
	}
	
	Xtal->Z = Z;
	printf("nmpc  = %d, napm = %d \n", Xtal->Z, Xtal->num_atoms_in_molecule);
	
	//remove_close_molecules(Xtal);
	


	free(mol_Xfrac);
	free(mol_Yfrac);
	free(mol_Zfrac);
	
}
*/
/*
void place_molecule_at_mirror_plane(crystal *Xtal, molecule *mol,
	int hall_number)
{
	load_database();
	
	//debug : check if spg supports inversion center
	if(database[hall_number-1].M.number_of_positions == 0)
	{
		printf("**mirror not supported by hall number %d**\n",
			hall_number);
		return ;
	}

	//declare variables
	//num of atoms in the molecule is N
	int N = mol->num_of_atoms;
	int Z = Xtal->Z;
	int num_atoms_in_cell = N*Z;
	
	//generate random angles for rotation
	//phi, the azimuthal angle
	float phi = 2*PI*((float)rand()/(float)(RAND_MAX));
	//theta, angle with x-axis. this should acos(u), -1<u<1
	float u = 2*((float)rand()/(float)(RAND_MAX)) - 1;
	float theta = acos(u);
	//choose axis
	float axis[]={cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)};
	//random rotation psi around axis 
	float psi = 2*PI*((float)rand()/(float)(RAND_MAX));
	
	//create rotation matix around the specified axis with angle psi
	float random_rotation_matrix[3][3] = {1,0,0,0,1,0,0,0,1};
	rotation_mat_around_axis(random_rotation_matrix, axis, psi);
	//randomly rotate
	molecule_rotate(mol, random_rotation_matrix);

	
	//get symmetry operations from spglib database
	double translations[192][3];
	int rotations[192][3][3];	
	int num_of_operations = spg_get_symmetry_from_database(rotations,
	translations, hall_number);
	Z = num_of_operations;

	//calculate fractional coordinates
	float inverse_lattice_vectors[3][3];
	float *mol_Xfrac = (float *)malloc(N*sizeof(float));
	float *mol_Yfrac = (float *)malloc(N*sizeof(float));
	float *mol_Zfrac = (float *)malloc(N*sizeof(float));
	inverse_mat3b3(inverse_lattice_vectors, Xtal->lattice_vectors);
	//mat3b3_transpose(inverse_lattice_vectors, inverse_lattice_vectors);
	//print_mat3b3(inverse_lattice_vectors);
	for(int i = 0; i < N; i++)
	{
		float temp[3] = {(*mol).X[i], (*mol).Y[i], (*mol).Z[i]};
		vector3_mat3b3_multiply(inverse_lattice_vectors, temp, temp);
		mol_Xfrac[i] = temp[0];
		mol_Yfrac[i] = temp[1];
		mol_Zfrac[i] = temp[2];		
	}
	
	//select random mirror plane and place first molecule
	//rot_mat and trans_vec specifies the location of plane 
	int random_select = 0;
	float rand_fracx = ((float)rand()/(float)(RAND_MAX));
	float rand_fracy = ((float)rand()/(float)(RAND_MAX));
	float rand_fracz = ((float)rand()/(float)(RAND_MAX));
	float rand_frac_array[3] = {rand_fracx,rand_fracy,rand_fracz};
	float rot_mat[3][3], trans_vec[3];
	copy_mat3b3_mat3b3bN(rot_mat, database[hall_number-1].M.first_position_rot,random_select);
	copy_vector3_vector3bN(trans_vec, database[hall_number-1].M.first_position_trans,random_select);
	vector3_mat3b3_multiply(rot_mat,rand_frac_array,rand_frac_array);
	vector3_add(trans_vec,rand_frac_array,rand_frac_array);	
	//debug
	print_vec3(rand_frac_array);
	
	//fractional coordinates of the first molecule placed on a mirror
	for(int i = 0; i < N; i++)
	{	
		mol_Xfrac[i] += rand_frac_array[0];
		mol_Yfrac[i] += rand_frac_array[1];
		mol_Zfrac[i] += rand_frac_array[2];	
		//printf("mol_frac = %f %f %f \n", mol_Xfrac[i], mol_Yfrac[i], mol_Zfrac[i]);
	}
	
	float lattice_vectors_transpose[3][3];
	copy_mat3b3_mat3b3(lattice_vectors_transpose, Xtal->lattice_vectors);
	mat3b3_transpose(lattice_vectors_transpose, lattice_vectors_transpose);
	
	//i loops over all the operations
	for(int i = 0; i < num_of_operations; i++)
	{
		int molecule_index = i*N; 
		int rot[3][3];
		copy_intmat3b3_intmat3b3bN(rot, rotations, i);

		
		//copy ith translations
		float trans[3];
		float temp_rot_mat[3][3];
		//copy_vector3_vector3bN(trans, translations, i);
		//printf("trans= %f %f %f\n", translations[i][0], translations[i][1], translations[i][2]);
		trans[0] = translations[i][0];
		trans[1] = translations[i][1];
		trans[2] = translations[i][2];
		copy_floatmat3b3_intmat3b3(temp_rot_mat, rot);
	
		//printf("i index = %d \n", i);
		//j loops over atoms in the molecule
		for(int j = 0; j < N; j++)
		{
			float atomj_array[3];
			atomj_array[0] = mol_Xfrac[j];
			atomj_array[1] = mol_Yfrac[j];
			atomj_array[2] = mol_Zfrac[j]; 
			//print_vec3(atomj_array);
			//print_mat3b3(rot);
			//printf("ist element = %d ", rot[0][1]);
			//vector3_intmat3b3_multiply(rot, atomj_array, atomj_array);
			vector3_mat3b3_multiply(temp_rot_mat, atomj_array, atomj_array);
			vector3_add(trans, atomj_array, atomj_array);
			//print_vec3(atomj_array);
			
			//temp_array has the fractional cord of jth atom in ith mol
			//convert to cartesian
			vector3_mat3b3_multiply(lattice_vectors_transpose,
				atomj_array, atomj_array);
			
			//copy to the structure
			Xtal->Xcord[molecule_index+j] = atomj_array[0];
			Xtal->Ycord[molecule_index+j] = atomj_array[1];
			Xtal->Zcord[molecule_index+j] = atomj_array[2];
			Xtal->atoms[molecule_index+j] = (*mol).atoms[j];
			
			//types[molecule_index+j] = 1;
			
			//printf("%f %f %f  \n", random_crystal.Xcord[molecule_index+j], random_crystal.Ycord[molecule_index+j], random_crystal.Zcord[molecule_index+j]);
			
			//printf("j index = %d \n", j);
			//printf("%d molecule %d atom = %f %f %f \n", i+1, j+1, mol->X[j], mol->Y[j], mol->Z[j] );
			//print_vec3(atomj_array);
		}
	}
	
	Xtal->Z = Z;
	//print_mat3b3bN(molecule_orientations,Z);
	
	//rot_mat has the information of orientation of the wyckoff position
	//use it to align the first molecule
	remove_close_molecules_mirror_plane(Xtal, rot_mat);

	free(mol_Xfrac);
	free(mol_Yfrac);
	free(mol_Zfrac);
	
}
*/
/*
void align_molecule_mirror_plane(crystal *Xtal, molecule *mol,
	int hall_number)
{
	
	load_database();
	
	//debug : check if spg supports  a mirror
	if(database[hall_number-1].M.number_of_positions == 0)
	{
		printf("**mirror not supported by hall number %d**\n",
			hall_number);
		return ;
	}
	
	//declare variables
	//num of atoms in the molecule is N
	int N = mol->num_of_atoms;
	int Z = Xtal->Z;
	int num_atoms_in_cell = N*Z;
	
	//select random mirror plane and place first molecule
	//rot_mat and trans_vec specifies the location of plane 
	int random_select = 0;
	float rand_fracx = ((float)rand()/(float)(RAND_MAX));
	float rand_fracy = ((float)rand()/(float)(RAND_MAX));
	float rand_fracz = ((float)rand()/(float)(RAND_MAX));
	float rand_frac_array[3] = {rand_fracx,rand_fracy,rand_fracz};
	float rot_mat[3][3], trans_vec[3];
	copy_mat3b3_mat3b3bN(rot_mat, database[hall_number-1].M.first_position_rot,random_select);
	copy_vector3_vector3bN(trans_vec, database[hall_number-1].M.first_position_trans,random_select);
	vector3_mat3b3_multiply(rot_mat,rand_frac_array,rand_frac_array);
	vector3_add(trans_vec,rand_frac_array,rand_frac_array);	
	//debug
	//print_vec3(rand_frac_array);
	
	//find normal from rot_mat
	float rot_mat_transpose[3][3], normal_frac[3];
	mat3b3_transpose(rot_mat_transpose,rot_mat);
	for(int i = 0; i < 3; i++)
	{
		float vec1[3],vec2[3];
		copy_vector3_mat3b3(vec1, rot_mat_transpose, i % 3);
		copy_vector3_mat3b3(vec2, rot_mat_transpose, (i+1) % 3);
		cross_vector3_vector3(normal_frac, vec1, vec2);
		if (vector3_norm(normal_frac) > 0.1 )
			break;
		if (i == 2)
			printf("***failed to find normal!!");
	}
	
	print_vec3(normal_frac);
	//find normal in cartesian
	float normal_cartesian[3];
	vector3_mat3b3_multiply(Xtal->lattice_vectors, normal_frac, normal_cartesian);
	print_vec3(normal_cartesian);
	
	//select mirror plane
	float molecular_axis[3][3]={{1,0,0},{0,1,0},{0,0,1}};
	int random_select2 = 0;
	float molecule_plane[3];
	copy_vector3_mat3b3(molecule_plane, molecular_axis, random_select2);
	
	//create rotation matrix to align molecule
	//1)find axis of rotation
	float rotation_axis[3];
	cross_vector3_vector3(rotation_axis, molecule_plane,
		normal_cartesian);
	normalise_vector3(rotation_axis);
	print_vec3(rotation_axis);
	//2)find angle of rotation
	double rotation_angle;
	normalise_vector3(molecule_plane);
	normalise_vector3(normal_cartesian);
	rotation_angle = acos(dot_vector3_vector3(molecule_plane,
		normal_cartesian));
	//3)find rotation_matrix
	float rotation_matrix[3][3];
	rotation_mat_around_axis(rotation_matrix, rotation_axis,
		rotation_angle);	
	
	//debug
	printf("rotation matrix = \n");
	print_mat3b3(rotation_matrix);
	
	//randomly rotate molecule along the molecule plane normal
	//generate random angles for rotation
	//random rotation psi around axis 
	float psi = 2*PI*((float)rand()/(float)(RAND_MAX));
	float random_rotation_matrix[3][3]; 
	rotation_mat_around_axis(random_rotation_matrix,
		molecule_plane, psi);
	//molecule_rotate(mol, random_rotation_matrix);
	
	//rotate molecule and place
	//calculate fractional coordinates of molecule at origin
	//ans store it in mol_Xfrac etc
	float inverse_lattice_vectors[3][3];
	float *mol_Xfrac = (float *)malloc(N*sizeof(float));
	float *mol_Yfrac = (float *)malloc(N*sizeof(float));
	float *mol_Zfrac = (float *)malloc(N*sizeof(float));
	inverse_mat3b3(inverse_lattice_vectors, Xtal->lattice_vectors);
	for(int i = 0; i < N; i++)
	{
		float atom_i[3] = {(*mol).X[i], (*mol).Y[i], (*mol).Z[i]};
		vector3_mat3b3_multiply(inverse_lattice_vectors,
			atom_i, atom_i);
		//vector3_mat3b3_multiply(rotation_matrix, atom_i, atom_i);
		vector3_add(atom_i, rand_frac_array, atom_i);
		mol_Xfrac[i] = atom_i[0];
		mol_Yfrac[i] = atom_i[1];
		mol_Zfrac[i] = atom_i[2];
		//printf("atom_frac %f %f %f C \n", atom_i[0],atom_i[1], atom_i[2]);
	}

	apply_all_symmetry_ops(Xtal, mol, mol_Xfrac, mol_Yfrac, mol_Zfrac, 
		N, hall_number);
	remove_close_molecules(Xtal);
	print_crystal(Xtal);
	
	free(mol_Xfrac);
	free(mol_Yfrac);
	free(mol_Zfrac);
	//exit(0);
}
*/
void apply_all_symmetry_ops(crystal *xtal,
							molecule *mol,
							float* mol_Xfrac,
							float *mol_Yfrac,
							float *mol_Zfrac,
							int N,
							int hall_number)
{
	//get symmetry operations from spglib database
	double translations[192][3];
	int rotations[192][3][3];	
	int num_of_operations = spg_get_symmetry_from_database(rotations,
	translations, hall_number);
	int Z = num_of_operations;
	xtal->Z = Z;
	
	//get lattice_vector inverse and transpose
	float inverse_lattice_vectors[3][3];
	inverse_mat3b3(inverse_lattice_vectors, xtal->lattice_vectors);
	float lattice_vectors_transpose[3][3];
	copy_mat3b3_mat3b3(lattice_vectors_transpose, xtal->lattice_vectors);
	mat3b3_transpose(lattice_vectors_transpose, lattice_vectors_transpose);
	
	//i loops over all the operations
	for(int i = 0; i < num_of_operations; i++)
	{
		int molecule_index = i*N; 
		
		//rot and trans are ith symmetry operation
		int rot[3][3];
		copy_intmat3b3_intmat3b3bN(rot, rotations, i);
		float trans[3];
		trans[0] = translations[i][0];
		trans[1] = translations[i][1];
		trans[2] = translations[i][2];
		
		//j loops over all atoms in a molecule
		for(int j = 0; j < N; j++)
		{
			float atomj_array[3];
			atomj_array[0] = mol_Xfrac[j];
			atomj_array[1] = mol_Yfrac[j];
			atomj_array[2] = mol_Zfrac[j]; 
			
			//apply operation
			vector3_intmat3b3_multiply(rot, atomj_array, atomj_array);
			vector3_add(trans, atomj_array, atomj_array);
			
			//temp_array has the fractional cord of jth atom in ith mol
			//convert to cartesian
			vector3_mat3b3_multiply(lattice_vectors_transpose,
				atomj_array, atomj_array);
			
			//copy to the structure
			xtal->Xcord[molecule_index+j] = atomj_array[0];
			xtal->Ycord[molecule_index+j] = atomj_array[1];
			xtal->Zcord[molecule_index+j] = atomj_array[2];
			xtal->atoms[2*(molecule_index+j)] = (*mol).atoms[2*j];
            xtal->atoms[2*(molecule_index+j)+1] = (*mol).atoms[2*j+1];	
		}	
	}
	bring_all_molecules_to_first_cell(xtal);
	xtal->num_atoms_in_molecule = mol->num_of_atoms;
}


/*function for placing molecule at any position for any spacegroup
 * molecule is automatically aligned according to the requirement 
 * of the wyckoff position. for general positions, molecule is rotated
 * randomly. generalizes all the functions written before.
 * assumes symmetry directions to be along cartesian axes of molecule
 * #TODO: create a molecule symm structure to store molecule symmetry
 */
void place_molecule_at_position(crystal *Xtal,
								molecule *mol,
								int hall_number,
								int spg,
								int wyckoff_pos)
{
	//load_database();
	
	//declare variables
	//num of atoms in the molecule is N
	int N = mol->num_of_atoms;
	int Z = Xtal->Z;
	//int num_atoms_in_cell = N*Z;
	//int general_position = 0;
	
	
	/*from a database get the primary and secondary axes of the
	 * wyckoff position; should be calculated from wyckoff alphabet 
	 * number. After getting the orientation from the database,
	 * multiply with lattice vectors to get directions in cartesian 
	 * axes. #TODO
	 */
	//for now assume general position  
	float primary_axis[3] = {spg_positions[spg].primary_orientation[wyckoff_pos][0],
							spg_positions[spg].primary_orientation[wyckoff_pos][1],
							spg_positions[spg].primary_orientation[wyckoff_pos][2]
							};
		 
	float secondary_axis[3] = {spg_positions[spg].secondary_orientation[wyckoff_pos][0],
							   spg_positions[spg].secondary_orientation[wyckoff_pos][1],
							   spg_positions[spg].secondary_orientation[wyckoff_pos][2]
							  };
	//From a database, get first coordinates where mol is to be placed
	//#TODO: get rot_mat and trans_vec from a databse of first position.
	//for now assume general position
	float rot_mat[3][3] = {{0,0,1},{0,1,0},{0,0,1}};
	copy_mat3b3_intmat3b3bN(rot_mat, spg_positions[spg].first_position_rot, wyckoff_pos);
	float trans_vec[3] = {0, 0 ,0};
	copy_vector3_vector3bN(trans_vec, spg_positions[spg].first_position_trans, wyckoff_pos);
	
	////check if general position
	//if( 	check_vec3_isNull(primary_axis, 1e-5) &&
			//check_vec3_isNull(secondary_axis, 1e-5) )
		//general_position = 1;
	
	/*This should be obtained from molecule symm structure. this is the
	 * primary and secondary axes of the molecule for the corresponding
	 * point group. for point group 1 or i, it doesnt matter.
	 * #TODO
	 */
	float mol_primary_axis[3] = {0, 0, 0};
	float mol_secondary_axis[3] = {0, 0, 0}; 
	
	//create random axes if primary is null-> only for 1 and i groups
	if( check_vec3_isNull(primary_axis, 1e-5) )
	{
		//phi, the azimuthal angle
		float phi = 2*PI*uniform_dist_01();
		//theta, angle with x-axis. this should acos(u), -1<u<1
		float u = 2*uniform_dist_01() - 1;
		float theta = acos(u);
		//compute random axis with theta and phi
		float axis[] = {cos(phi)*sin(theta),
						sin(phi)*sin(theta),
						cos(theta)};
		copy_vector3_vector3(primary_axis, axis);
		mol_primary_axis[0] = 0;
		mol_primary_axis[1] = 0;
		mol_primary_axis[2] = 1;
	}

	
	//rotate molecule s.t. primary axis coincides with the primary axis 
	//of the molecule
	float rotation_matrix[3][3];
	rotation_matrix_from_vectors(rotation_matrix,
								 primary_axis,
								 mol_primary_axis);
	molecule_rotate(mol, rotation_matrix);
	//debug 
	//printf("det of R1 = %f \n",det_mat3b3(rotation_matrix) );

	//rotate molecule along the the primary axis s.t. the secondary axis
	//of the molecule coincides with secondary of wyckoff position
	//if secondary is null, choose random angle.
	float angle;
	if (check_vec3_isNull(secondary_axis, 1e-5))
		angle = 2*PI*uniform_dist_01();
	else
	{
		normalise_vector3(mol_secondary_axis);
		normalise_vector3(secondary_axis);
		angle = acos(dot_vector3_vector3(secondary_axis,
										 mol_secondary_axis));
	}
	rotation_mat_around_axis(rotation_matrix, primary_axis, angle);
	molecule_rotate(mol, rotation_matrix);
	
	//debug 
	//printf("det of R2 = %f \n",det_mat3b3(rotation_matrix) );
	
	//for converting to fractional coordinates and
	// place the first molecule
	float inverse_lattice_vectors[3][3];
	inverse_mat3b3(inverse_lattice_vectors, Xtal->lattice_vectors);
	mat3b3_transpose(inverse_lattice_vectors,inverse_lattice_vectors);
	float mol_Xfrac[N]; //stores fractional coordinates of first mol
	float mol_Yfrac[N];
	float mol_Zfrac[N];
	

	float rand_fracx = uniform_dist_01();
	float rand_fracy = uniform_dist_01();
	float rand_fracz = uniform_dist_01();;
	float rand_frac_array[3] = {rand_fracx,rand_fracy,rand_fracz};
	vector3_mat3b3_multiply(rot_mat,rand_frac_array,rand_frac_array);
	vector3_add(trans_vec,rand_frac_array,rand_frac_array);	
	//rand_frac_array is the postion of the first mol
	
	//compute fractional coordiates of the first position and
	//move molecule to the position
	for(int i = 0; i < N; i++)
	{
		float atom_i[3] = {(*mol).X[i], (*mol).Y[i], (*mol).Z[i]};
		vector3_mat3b3_multiply(inverse_lattice_vectors,
								atom_i,
								atom_i);
		vector3_add(atom_i, rand_frac_array, atom_i); //translate
		mol_Xfrac[i] = atom_i[0];
		mol_Yfrac[i] = atom_i[1];
		mol_Zfrac[i] = atom_i[2];

	} 
	//now mol_frac has the first mol in fractional coordinates
	// at the given wyckoff position.
	

	//now apply all the symmetry operations of spg
	apply_all_symmetry_ops(Xtal,
						   mol,
						   mol_Xfrac,
						   mol_Yfrac,
						   mol_Zfrac,
						   N,
						   hall_number);
						 

	//exit(0);
	//remove close molecules sitting on top of each other
	//remove_close_molecules(Xtal);
	//print_crystal(Xtal);

	//combine_close_molecules(Xtal);
	//debug
	//print_crystal(Xtal);
	
	int order = Z / Xtal->Z;
	//printf("placed molecule in Wyckoff position of order: %d\n",order);
}


/*function for placing molecule at any position for any spacegroup
 * molecule is automatically aligned according to the requirement 
 * of the wyckoff position. for general positions, molecule is rotated
 * randomly. generalizes all the functions written before.
 * assumes symmetry directions to be along cartesian axes of molecule
 * #TODO: create a molecule symm structure to store molecule symmetry
 */
int auto_align_and_generate_at_position(crystal *Xtal,
								molecule *mol,
								int hall_number,
								int spg,
								int wyckoff_pos,
								float *mol_axes,
								int num_axes, 
								int* overlap_list)
{

	//declare variables
	//num of atoms in the molecule is N
	int N = mol->num_of_atoms;
	int Z = Xtal->Z;
	int order = spg_positions[spg-1].multiplicity[0]/\
		spg_positions[spg-1].multiplicity[wyckoff_pos];
	//int num_atoms_in_cell = N*Z;
	//int general_position = 0;
	//From a database, get first coordinates where mol is to be placed
	//#TODO: get rot_mat and trans_vec from a databse of first position.
	float rot_mat[3][3] = {{0,0,1},{0,1,0},{0,0,1}};
	copy_mat3b3_intmat3b3bN(rot_mat, 
							spg_positions[spg-1].first_position_rot,
							wyckoff_pos);
	float trans_vec[3] = {0, 0 ,0};
	copy_vector3_vector3bN(trans_vec,
						  spg_positions[spg-1].first_position_trans,
						  wyckoff_pos);

	//if general position, rotate molecule randomly
	if ( get_degrees_of_freedom(spg, wyckoff_pos) == 2)
	{
		//create rotation matix around the specified axis with angle psi
		float random_rotation_matrix[3][3];
		generate_random_rotation_matrix(random_rotation_matrix);
		//randomly rotate
		molecule_rotate(mol, random_rotation_matrix);
	}
	
	/*
	
	//create random axes if primary is null-> only for 1 and i groups
	if( check_vec3_isNull(primary_axis, 1e-5) )
	{
		//phi, the azimuthal angle
		float phi = 2*PI*uniform_dist_01();
		//theta, angle with x-axis. this should acos(u), -1<u<1
		float u = 2*uniform_dist_01() - 1;
		float theta = acos(u);
		//compute random axis with theta and phi
		float axis[] = {cos(phi)*sin(theta),
						sin(phi)*sin(theta),
						cos(theta)};
		copy_vector3_vector3(primary_axis, axis);
		mol_primary_axis[0] = 0;
		mol_primary_axis[1] = 0;
		mol_primary_axis[2] = 1;
	}

	
	//rotate molecule s.t. primary axis coincides with the primary axis 
	//of the molecule
	float rotation_matrix[3][3];
	rotation_matrix_from_vectors(rotation_matrix,
								 primary_axis,
								 mol_primary_axis);
	molecule_rotate(mol, rotation_matrix);
	//debug 
	//printf("det of R1 = %f \n",det_mat3b3(rotation_matrix) );
	
	//rotate molecule along the the primary axis s.t. the secondary axis
	//of the molecule coincides with secondary of wyckoff position
	//if secondary is null, choose random angle.
	float angle;
	if (check_vec3_isNull(secondary_axis, 1e-5))
		angle = 2*PI*uniform_dist_01();
	else
	{
		normalise_vector3(mol_secondary_axis);
		normalise_vector3(secondary_axis);
		angle = acos(dot_vector3_vector3(secondary_axis,
										 mol_secondary_axis));
	}
	rotation_mat_around_axis(rotation_matrix, primary_axis, angle);
	molecule_rotate(mol, rotation_matrix);
	*/
	//debug 
	//printf("det of R2 = %f \n",det_mat3b3(rotation_matrix) );
	
	//for converting to fractional coordinates and
	// place the first molecule
	float inverse_lattice_vectors[3][3];
	inverse_mat3b3(inverse_lattice_vectors, Xtal->lattice_vectors);
	mat3b3_transpose(inverse_lattice_vectors,inverse_lattice_vectors);
	float mol_Xfrac[N]; //stores fractional coordinates of first mol
	float mol_Yfrac[N];
	float mol_Zfrac[N];
	
	float rand_frac_array[3] = {uniform_dist_01(),\
								uniform_dist_01(),\
								uniform_dist_01()};
	vector3_mat3b3_multiply(rot_mat,rand_frac_array,rand_frac_array);
	vector3_add(trans_vec,rand_frac_array,rand_frac_array);	
	//rand_frac_array is the postion of the first mol
	
	//compute fractional coordiates of the first position and
	//move molecule to the position
	for(int i = 0; i < N; i++)
	{
		float atom_i[3] = {(*mol).X[i], (*mol).Y[i], (*mol).Z[i]};
		vector3_mat3b3_multiply(inverse_lattice_vectors,
								atom_i,
								atom_i);
		vector3_add(atom_i, rand_frac_array, atom_i); //translate
		mol_Xfrac[i] = atom_i[0];
		mol_Yfrac[i] = atom_i[1];
		mol_Zfrac[i] = atom_i[2];

	} 
	//now mol_frac has the first mol in fractional coordinates
	// at the given wyckoff position.
	

	//now apply all the symmetry operations of spg
	apply_all_symmetry_ops(Xtal,
						   mol,
						   mol_Xfrac,
						   mol_Yfrac,
						   mol_Zfrac,
						   N,
						   hall_number);
	
	//check if inv center or genral position
	int dof = get_degrees_of_freedom(spg, wyckoff_pos);
	if ( dof == 2)
	{	
		//if inv centre remove overlap molecules
		if( strcmp(spg_positions[spg-1].site_symmetry[wyckoff_pos], "-1") == 0 ) 
			combine_close_molecules(Xtal);
		return 1;
	}
						   
	Xtal->spg = spg;		 
	bring_molecules_to_origin(Xtal);

	int result = align_using_std_orientations(Xtal, mol, hall_number, 
		mol_axes, num_axes, rand_frac_array, overlap_list, order, dof);
		
		
	if (result == 0)
	{
		//printf("spg: %d , pos: %d unable to align\n", spg, wyckoff_pos);
		return 0;
	}
	
    //printf("spg: %d , pos: %d aligned\n", spg, wyckoff_pos);
	return 1;
	//exit(0);
	//remove close molecules sitting on top of each other
	//remove_close_molecules(Xtal);
	//print_crystal(Xtal);

	//combine_close_molecules(Xtal);
	//debug
	//print_crystal(Xtal);
	
	//int order_detected = Z / Xtal->Z;
	//printf("placed molecule in Wyckoff position of order: %d\n",order);
}


int align_using_std_orientations(crystal* xtal_1,\
								molecule* mol,\
							  int hall_number,\
							  float *mol_axes,\
							  int num_axes,\
							  float first_com[3],\
							  int overlap_list[],
							  int len_overlap_list,
							  int dof)
{
	//take two pairs
	//rotate molecule to average position
	int N = mol->num_of_atoms;
	int Z = xtal_1->Z;
	float rotation_matrix[3][3];
	float mol_Xfrac[N]; //stores fractional ...
	float mol_Yfrac[N]; //coordinates of first mol
	float mol_Zfrac[N];
	float lattice_vectors[3][3];
	float inv_lattice_vectors[3][3];
	copy_mat3b3_mat3b3(lattice_vectors, xtal_1->lattice_vectors);
	inverse_mat3b3(inv_lattice_vectors, lattice_vectors);
	mat3b3_transpose(inv_lattice_vectors, inv_lattice_vectors);
	
	//create a temp xtal
	crystal *xtal = (crystal *)malloc(sizeof(crystal) );
	allocate_xtal(xtal, Z, N);
	copy_xtal(xtal, xtal_1);
	//debug
	int spg = xtal_1->spg;
	//generate_lattice(xtal->lattice_vectors, spg, 60, 120, 1000000);
	//generate_fake_lattice(xtal->lattice_vectors, spg);
	
	//rotate about axis if allowed
	float random_rotation_about_axis[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
	
	//pick a mol axis
	array_shuffler_1(mol_axes, num_axes );
	for (int i = 0; i < num_axes; i++)
	{
		float mol_axis[3] = {*(mol_axes+3*i + 0),
							 *(mol_axes+3*i + 1),
							 *(mol_axes+3*i + 2) };
		//pick a viewing direction 
		array_shuffler_2(viewing_directions, num_viewing_direction);
		for (int j = 0; j < num_viewing_direction; j++)
		{
			float view_dir[3] = {viewing_directions[j][0],
								 viewing_directions[j][1],
								 viewing_directions[j][2] };
			normalise_vector3(view_dir);
			
			//create rotation matrix
			rotation_matrix_from_vectors(rotation_matrix,mol_axis, view_dir);
			//rotate and store first molecule in fractional coord 
			//in mol_frac
			
			if(dof == 1)
			{
				float psi = 2*PI*uniform_dist_01();
				rotation_mat_around_axis(random_rotation_about_axis, view_dir, psi);
			}
			
			
			for(int z = 0; z < N; z++)
			{
				float temp[3] = {xtal_1->Xcord[z] ,
								 xtal_1->Ycord[z] ,
								 xtal_1->Zcord[z] };
				vector3_mat3b3_multiply(rotation_matrix, temp, temp);
				
				vector3_mat3b3_multiply(random_rotation_about_axis, temp, temp);
				//vector3_add(temp,com,temp);
				//convert to frac
				vector3_mat3b3_multiply(inv_lattice_vectors,
										temp,
										temp );
				vector3_add(temp, first_com, temp);
				mol_Xfrac[z] = temp[0];
				mol_Yfrac[z] = temp[1];
				mol_Zfrac[z] = temp[2];
			}
			
			apply_all_symmetry_ops(	xtal,
									mol,
									mol_Xfrac,
									mol_Yfrac,
									mol_Zfrac,
									N,
									hall_number);
			
			//bring_molecules_to_origin(xtal);
			
			int result = check_overlap_xtal(xtal,
											overlap_list,
											len_overlap_list,
											N);
					
			if (result)
			{
				int Z_gen = xtal->Z;
				//printf("i=%d, j1=%d, j2=%d, k1=%d, k2=%d \n", i,j1,j2,k1,k2);
				//remove_close_molecules(xtal);
				combine_close_molecules(xtal);
				int Z_return = xtal->Z;
				//len_overlap_list == order, Z/len = multiplicity
				if(Z_return != Z_gen/len_overlap_list)
					{continue;}

				//detect_spg_using_spglib(xtal);
				//print_crystal(xtal);
				//printf("result 1\n");
				copy_xtal(xtal_1, xtal);
				free_xtal(xtal);
				return Z_return;
			}
			  
		}//viewing dir
	}//end of loop mol axis
	
	free_xtal(xtal);
	return 0;
}



int get_degrees_of_freedom(int spg, int pos)
{
	char pg[6];
	strcpy(pg, spg_positions[spg-1].site_symmetry[pos]);
	//general_position or inversion_centre: complete freedom
	if(!strcmp(pg,"-1") || !strcmp(pg,"1") )
		return 2;
	else if (!strcmp(pg,"2") || !strcmp(pg,"3") || !strcmp(pg,"m") ||
			 !strcmp(pg,"4") || !strcmp(pg,"6") || !strcmp(pg, "2/m") ||
			 !strcmp(pg,"-4")|| !strcmp(pg,"4/m") || !strcmp(pg,"6/m") ||
			 !strcmp(pg,"-3") || !strcmp(pg,"-6")
			 )
		{return 1;
			}
	else
		{return 0;}
	
}
