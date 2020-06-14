#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
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



int generate_crystal(crystal* random_crystal, molecule* mol,float volume,
	float Z, float Zp_max, int spg, COMPATIBLE_SPG compatible_spg[],
	int len_compatible_spg, int compatible_spg_index)
{
	Zp_max = 192; //stupid argument
    len_compatible_spg += 1; // not needed now; stupid argument 
	
    //crystal random_crystal;
	float max_angle = 30 * PI/180;
	float min_angle = 150 * PI/180;

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
	
	generate_lattice(random_crystal->lattice_vectors, spg, max_angle, min_angle, volume);
	
	//find a random pos
	int pos_index = rand_r(seed2) % compatible_spg[compatible_spg_index].num_allowed_pos;
	int pos = compatible_spg[compatible_spg_index].allowed_pos[pos_index];
	random_crystal->wyckoff_position = pos;
	
	//place, align and attempt to generate crystal at position pos
	int result = auto_align_and_generate_at_position(random_crystal,
							mol,
							hall_number,
							spg, 
							pos_index,
							compatible_spg[compatible_spg_index]);
	//copy back to mol
	copy_positions_to_mol(mol, Xm, Ym, Zm);
	
	if(!result)
	{
		
		return 0;
	}
	else
		return 1;
}

int find_num_structure_for_spg(int num_structures, char spg_dist_type[10], int spg, int Z)
{
	if ( strcmp(spg_dist_type, "uniform") == 0)
		return num_structures;

	else if ( strcmp(spg_dist_type, "standard") == 0)
	{
		int order  = spg_positions[spg-1].multiplicity[0]/ Z;
		return (num_structures/order);
	}

	// list from Genarris 1.0 source code
	else if ( strcmp(spg_dist_type, "chiral") == 0 )
	{
		if (  spg == 1 || (spg > 2 && spg < 6) || (spg > 15 && spg < 25) ||
			 (spg > 74 && spg < 81) || (spg > 88 && spg < 99) || (spg > 142 && spg < 147) ||
			 (spg > 148 && spg < 156) || (spg > 167 && spg < 174) ||(spg > 176 && spg < 183) ||
			 (spg > 194 && spg < 200) || (spg > 207 && spg < 215)
			)
		{
			return num_structures;
		}

	 	else
	 	{
	 		return 0;
	 	}
	}

	else if ( strcmp(spg_dist_type, "csd") == 0 )
	{
		const int len = 10;
		// list of high frequency spgs in csd
		// http://pd.chem.ucl.ac.uk/pdnn/symm3/sgpfreq.htm
		int csd_list[ ] = { 14, 2, 19, 15, 3, 61, 62, 33, 9, 1};
		for (int i = 0; i < len; i++)
		{
			if (csd_list[i] == spg)
				return num_structures;
		}

		return 0;
	}

	else
	{
		printf("***ERROR: spg_generation: spg_dist_type not found\n");
		exit(0);
	}
}