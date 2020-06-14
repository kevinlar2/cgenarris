#ifndef _CHECK_STRUCTURE_H
#define _CHECK_STRUCTURE_H

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
			
int structure_checker(float L[3][3],
	float *X ,
	int dim1,
	float *Y,
	int dim2,
	float *Z,
	int dim3,
	float *vdw_cutoff,
	int dim4, 
	int dim5,
	int *mol_id,
	int num_mols
	);

int fast_screener_vdw(crystal xtal, float *vdw_matrix);

int check_pairwise_if_mol_close(float L[3][3],
 float *vdw_matrix, int dim, float *X, float *Y,
 float *Z, int id1, int N1, int id2, int N2,
 float com1[3], float com2[3], float sum_len);

int check_pairwise(float L[3][3], float *vdw_matrix,
int dim, float *X, float *Y, float *Z, int id1,
int id2, int N1, int N2, int i, int j, int k);

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



#endif
