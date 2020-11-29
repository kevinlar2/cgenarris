#ifndef _CHECK_STRUCTURE_H
#define _CHECK_STRUCTURE_H

typedef struct 
{
	float (*L)[3];
	float *com1;
	float *com2;
	float *X;
	float *Y;
	float *Z;
	int index1;
	int num_atoms1;
	int index2;
	int num_atoms2;

}xtal_molecule_pair;


int check_structure(crystal random_crystal, float sr);
int check_structure_with_vdw_matrix(crystal random_crystal,
	float *vdw_matrix,
	int dim1,
	int dim2);

//void vector_cpy(float A[], float B[][3], int index);
float pdist(float T[3][3],
			float T_inv[3][3],
			float x1, 
			float x2,
			float x3,
			float y1,
			float y2,
			float y3  );
			
int structure_checker(crystal *xtal,
	float *vdw_cutoff,
	int total_atoms,
	int *mol_id,
	int num_mols
	);

int fast_screener_vdw(crystal xtal, float *vdw_matrix);

int check_pairwise_if_mol_close( float *vdw_matrix, int total_atoms,
	xtal_molecule_pair *mol_pair, float max_dist);

int check_pairwise(float *vdw_matrix , int total_atoms,
	xtal_molecule_pair *mol_pair, float l_disp[3]);

int check_lower_triangular(float L[3][3]);

float find_mol_len(float *X, float *Y, float *Z, int len);

void find_mol_com(float *X, float *Y, float *Z, int len, float com[3]);

float pdist_appx(float T[3][3],
			float T_inv[3][3],
			float x1, 
			float x2,
			float x3,
			float y1,
			float y2,
			float y3  );

void create_vdw_matrix_from_sr( molecule *mol,
								float *vdw_matrix,
								float sr,
								int Z);
/*
static int xseq_generator(int reset);
static int yseq_generator(int reset);
static int zseq_generator(int reset);
*/



#endif
