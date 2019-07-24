#ifndef SPG_GENERATOR_H
#define SPG_GENERATOR_H



//for storing allowed spacegroup and position
typedef struct 
{
	unsigned int spg;
	unsigned int num_allowed_pos;
	unsigned int *allowed_pos;
	int *pos_overlap_list[16];
	
}COMPATIBLE_SPG;

typedef struct
{
	float lattice_vectors[3][3];
	int num_atoms_in_molecule;
	int Z;
	int Zp;
	float *Xcord;
	float *Ycord;
	float *Zcord;
	char *atoms;
	int spg;
	int wyckoff_position;
	
}crystal;



int generate_crystal(crystal* random_crystal, molecule* mol,float volume,
	float Z, float Zp_max, int spg, COMPATIBLE_SPG compatible_spg[],
	int len_compatible_spg, int compatible_spg_index, float* mol_axes,
	int num_axes);




#endif
