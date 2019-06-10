#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include "read_input.h"
#include "spg_generation.h"
#include "combinatorics.h"
#include "lattice_generator.h"
#include "algebra.h"
#include "crystal_utils.h"
#include "molecule_utils.h"
#include "molecule_placement.h"
#include "spglib.h"

#define PI 3.141592653

extern unsigned int *seed2;
#pragma omp threadprivate(seed2)


int generate_crystal(crystal* random_crystal, molecule* mol,float volume,
	float Z, float Zp_max, int spg, COMPATIBLE_SPG compatible_spg[],
	int len_compatible_spg, int compatible_spg_index, float* mol_axes,
	int num_axes)
{
	
	//crystal random_crystal;
	float max_angle = 30 * PI/180;
	float min_angle = 150 * PI/180;

	//recenter molecule to origin
	recenter_molecule(mol);
	random_crystal->Z = Z;
	int N = mol->num_of_atoms;

	//copy molecules to an array to save it. molecule might deform
	//upon many rotations.	
	float Xm[N]; //molecule X coordinate
	float Ym[N];
	float Zm[N];
	copy_positions_to_array(mol, Xm, Ym, Zm);

	int hall_number;
	hall_number = hall_number_from_spg(spg);
	//printf("hall = %d \n", hall_number);
	//printf("attempted spg = %d \n", spg);
	//Z=4;
	//returns lattice vector for given spg and volume 
	generate_lattice(random_crystal->lattice_vectors, spg, max_angle, min_angle, volume);
	
	//find a random pos
	int pos_index = rand_r(seed2) % compatible_spg[compatible_spg_index].num_allowed_pos;
	
	int pos = compatible_spg[compatible_spg_index].allowed_pos[pos_index];
	
	//place, align and attempt to generate crystal at position pos
	int result = auto_align_and_generate_at_position(random_crystal, mol, hall_number, spg, 
		pos,mol_axes, num_axes, compatible_spg[compatible_spg_index].pos_overlap_list[pos_index]);
	

	
	
	
	
	//random_crystal->Z = Z;
	/*
	random_crystal->lattice_vectors[0][0] =45;
	random_crystal->lattice_vectors[0][1] =0;
	random_crystal->lattice_vectors[0][2] =0;
	
	random_crystal->lattice_vectors[1][0] =0;
	random_crystal->lattice_vectors[1][1] =30;
	random_crystal->lattice_vectors[1][2] =0;
	
	random_crystal->lattice_vectors[2][0] =20;
	random_crystal->lattice_vectors[2][1] =0;
	random_crystal->lattice_vectors[2][2] = 20;
	*/
		//align_molecule_mirror_plane(random_crystal, mol,hall_number);
	//print_molecule(mol);
	//place_molecule_at_general_position(random_crystal, mol, hall_number);
	//print_crystal(random_crystal);
	//print_molecule(mol);
	//exit(0);
	//place_molecule_at_inversion_center(random_crystal, mol, hall_number);
	//place_molecule_at_mirror_plane(random_crystal, mol, hall_number);
	//place_molecule_at_position(random_crystal, mol, hall_number,5,0);
	copy_positions_to_mol(mol, Xm, Ym, Zm);
	
	if(result == 0)
		return 0;
	else
		return 1;
	//printf("hall =%d ,nmpc  = %d, napm = %d \n", hall_number, random_crystal->Z, random_crystal->num_atoms_in_molecule);

	
	                   /*spglib check*/
	/*                
	convert_xtal_to_fractional(random_crystal);				
	int num_atoms_in_cell = random_crystal->Z * random_crystal-> num_atoms_in_molecule ;
	printf("num in cell = %d \n", num_atoms_in_cell);
	int types[num_atoms_in_cell];
	double positions[num_atoms_in_cell][3];
	for(int i = 0; i < num_atoms_in_cell; i++)
	{
		positions[i][0] = random_crystal->Xcord[i]; 
		positions[i][1] = random_crystal->Ycord[i];
		positions[i][2] = random_crystal->Zcord[i];
		types[i] = 1;
	}
	
	char symbol[21];
	int num_spg;
	double lattice_vector[3][3];
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++ )
		{
				lattice_vector[i][j] = random_crystal->lattice_vectors[j][i];
				//printf("%d, %d = %f /n ", i, j,	lattice_vector[i][j]);
		}
	
	
	
	for(int i = 0; i < num_atoms_in_cell; i++)
		printf("atom_frac %f %f %f  %d\n", positions[i][0], positions[i][1], positions[i][2], types[i]);
	
	
	
	num_spg = spg_get_international(symbol, lattice_vector, positions, types, num_atoms_in_cell, 1e-2);
	
	SpglibError error;
	error = spg_get_error_code();
	printf("%s\n", spg_get_error_message(error));
	
	printf("detected space group = %d \n", num_spg);
	*/
	//add info to rand crystal
	//print_crystal_fractional(random_crystal);
	
	//print_mat3b3(random_rotation_matrix);
	/*
	print_mat3b3(random_crystal.lattice_vectors);
	print_mat3b3(inverse_lattice_vectors);
	//printf("det = %f \n", det_mat3b3(random_rotation_matrix) );
	*/
	//exit(0);
	
}

