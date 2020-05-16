#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef struct
{
	float lattice_vectors[3][3];
	float *Xcord;
	float *Ycord;
	float *Zcord;
	char *atoms;
	int spg;
	int wyckoff_position;
    int num_atoms_in_molecule;
	int Z;
	int Zp;
	
}crystal;

void print_crystal(crystal* xtal)
{
	
	
	 int N = xtal->num_atoms_in_molecule;
	 int m = xtal->Z;
	 
	 printf("#Z = %d \n", xtal->Z);
	 printf("#napm = %d \n", xtal->num_atoms_in_molecule);
	
	for(int i = 0; i < 3; i++)
	{
		printf("lattice_vector %12f %12f %12f \n",
			xtal->lattice_vectors[i][0], xtal->lattice_vectors[i][1],
			xtal->lattice_vectors[i][2]);
	}
	
	for(int i = 0; i < N*m; i++)
	{
		printf("atom %12f %12f %12f %4c \n", xtal->Xcord[i],
			xtal->Ycord[i],  xtal->Zcord[i],  xtal->atoms[2*i]);
	}
}


void read_xtal_from_file(char *fname, crystal* xtal)
{
	FILE *fileptr;
	size_t len = 0;
	char *line = NULL;
	char *sub_line = NULL;
	int read;
	fileptr = fopen(fname,"r");

	if(!fileptr)
	{
		printf("***ERROR: xtal file not found \n");
		exit(0);
	}

	int total_atoms = 0;
	while ((read = getline(&line, &len, fileptr)) != -1)
	{
        if (strstr(line, "#") != NULL)
            continue;

		sub_line=strtok(line," ");
        //printf("%s \n" , sub_line);
        if(strcmp(sub_line, "atom") == 0)
            total_atoms++;
        else
            continue;
    }
    fclose(fileptr);

    int N = total_atoms;
    printf("total atoms = %d\n", total_atoms);
    xtal->Xcord = (float *)malloc(N*sizeof(float));
    
    xtal->Ycord = (float *)malloc(N*sizeof(float));
    xtal->Zcord = (float *)malloc(N*sizeof(float));
    xtal->atoms = (char *)malloc(2*N*sizeof(char));

    fileptr = fopen(fname,"r");
    int i = 0; 
    total_atoms = 0;

    while ((read = getline(&line, &len, fileptr)) != -1)
	{
        if (strstr(line, "#") != NULL)
            continue;

		sub_line=strtok(line," ");
        
        if(strcmp(sub_line, "atom") == 0)
         {
        	sub_line=strtok(NULL," ");
			(*xtal).Xcord[total_atoms]=atof(sub_line);
			sub_line=strtok(NULL," ");
			(*xtal).Ycord[total_atoms]=atof(sub_line);
			sub_line=strtok(NULL," ");
			(*xtal).Zcord[total_atoms]=atof(sub_line);
			sub_line=strtok(NULL," ");
			(*xtal).atoms[2*total_atoms]=*sub_line;
        	if(*(sub_line+1) == '\n' || *(sub_line+1) == ' ' || *(sub_line+1) == '\0' )
            	(*xtal).atoms[2*total_atoms+1]=' ';            
        	else
            	(*xtal).atoms[2*total_atoms+1]=*(sub_line+1);

            total_atoms++;
        }
        
        if (strcmp(sub_line, "lattice_vector") == 0)
        {
        	sub_line=strtok(NULL," ");
        	xtal->lattice_vectors[i][0] = atof(sub_line);
        	sub_line=strtok(NULL," ");
        	xtal->lattice_vectors[i][1] = atof(sub_line);
        	sub_line=strtok(NULL," ");
        	xtal->lattice_vectors[i][2] = atof(sub_line);

        	i++;
        } 

        else
            continue;
   	}
    fclose(fileptr);

    xtal->Z = 0;
    xtal->spg = 0;

}

int main()
{
    crystal xtal;
    float temp = 0;
    read_xtal_from_file("test.in", &xtal);
    xtal.num_atoms_in_molecule = 30;
    xtal.Z = 2;
    print_crystal(  &xtal);
    int mol_id[] = {0, 30};
    int num_mols = 2;
    structure_checker(xtal.lattice_vectors, xtal.Xcord, 60, xtal.Ycord, 60, xtal.Zcord, 60, &temp, 60, 60, mol_id, num_mols);

}


#define CONSTANT_TOL 0.001
#define MAXVDW 2.1
#define NORM2(a,b) ( (a[0]-b[0])*(a[0]-b[0]) +\
 					 (a[1]-b[1])*(a[1]-b[1]) +\
  					 (a[2]-b[2])*(a[2]-b[2]) )

#define NORMSQR(a) ( a[0]*a[0] + a[1]*a[1] + a[2]*a[2] )

float find_mol_len(float *X, float *Y, float *Z, int len);
int check_pairwise_if_mol_close(float L[3][3], float *vdw_matrix,\
 float *X, float *Y, float *Z, int id1, int N1, int id2, int N2,\
 float com1[3], float com2[3], float sum_len);


/* for pygenarris API only
* Assume lattice to be lower triangular in row major form
*/
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
	)
{
	//preliminary checks
	if(dim1 != dim2 || dim2 != dim3 || dim3 != dim1)
	{
		printf("***ERROR: length of coordinate arrays don't match\n");
		return -1;
	}

	if(dim4 != dim5)
	{
		printf("***ERROR: vdw cutoff matrix should be a square matrix\n");
		return -1;
	}

	if(!check_lower_triangular(L))
	{
		printf("***ERROR: lattice_vectors should be a lower triangular matrix in row major form\n");
		return -1;
	}	

	//find number of atoms in each molecule
	int num_atoms_in_molecule[num_mols];
	for(int i = 0; i < num_mols -1; i++)
	{
		num_atoms_in_molecule[i] = mol_id[i+1] - mol_id[i];
		printf("%s = %d\n", "num_atoms_molecule = ", num_atoms_in_molecule[i]);
	}
	//for last mol:
	num_atoms_in_molecule[num_mols-1] = dim1 - mol_id[num_mols-1];

	// estimate molecule lengths
	float mol_len[num_mols];
	for(int i = 0 ; i < num_mols; i++)
	{
		mol_len[i] = find_mol_len(X + mol_id[i],
								 Y + mol_id[i],
								 Z + mol_id[i],
								 num_atoms_in_molecule[i]);
		printf("molLen = %f\n", mol_len[i]);
	}

	// find com
	float com[num_mols * 3];
	for(int i = 0; i < num_mols; i++)
	{
		printf("num_atoms_in_molecule = %d\n", num_atoms_in_molecule[i]);
		find_mol_com(X + mol_id[i], Y +  mol_id[i], Z +  mol_id[i], num_atoms_in_molecule[i], com + 3*i);
		printf("i = %d , com  === %f, %f, %f\n",i, com[3*i], com[3*i +1], com[3*i+2] );
	}

	for(int i = 0; i < num_mols; i++)
	{
		for(int j = i; j < num_mols; j++)
		{
			//check pair of selcted molecules.
			printf("i, j = %d, %d \n", i, j );
			float sum_len = mol_len[i] + mol_len[j];
			check_pairwise_if_mol_close(L, vdw_cutoff, X, Y, Z, mol_id[i], num_atoms_in_molecule[i], mol_id[j], num_atoms_in_molecule[j], com+3*i, com+3*j, sum_len);
			exit(0);
		}
	}

}


int check_pairwise_if_mol_close(float L[3][3], float *vdw_matrix, float *X, float *Y, float *Z, int id1, int N1, int id2, int N2, float com1[3], float com2[3], float sum_len)
{
	int xmax = sum_len/L[0][0] + 1;
	int ymax = (sum_len + xmax*L[1][0])/L[1][1] + 1;
	int zmax = (sum_len + xmax*L[2][0] + ymax*L[2][0])/L[2][2] + 1;
	printf("xmax, ymax, zmax = %d %d %d\n", xmax, ymax, zmax);

	for(int i = -xmax; i < xmax; i++)
	for(int j = -ymax; j < ymax; j++)
	for(int k = -zmax; k < zmax; k++)
	{
		float dist[3] = {com1[0] - com2[0] + i*L[0][0] + j*L[1][0] + k*L[2][0], 
						 com1[1] - com2[1]             + j*L[1][1] + k*L[2][1],
						 com1[2] - com2[2]                         + k*L[2][2]};

		
		if (NORMSQR(dist) < (sum_len/2 + MAXVDW)*(sum_len/2 + MAXVDW))
		{
			printf("i, j, k = %d, %d, %d dist = %f %f %f\n",i, j, k, dist[0], dist[1], dist[2] );
			printf("sumlen = %f\n", sum_len);
			printf("dist = %f\n", NORMSQR(dist));
			printf("close; ijk = %d %d %d\n", i,j,k );
			//int ret_val = check_pairwise(L, vdw_matrix, X, Y, Z, id1, id2, N1, N2, i, j, k);
			//if(!ret_val)
			//	return 0;
		}
		
	}
	return 1;
}

int check_pairwise(float L[3][3], float *vdw_matrix, float *X, float *Y, float *Z, int id1, int N1, int id2, int N2, int i, int j, int k)
{
	for(int i1 = 0; i1 < N1; i1++)
	for(int i2 = 0; i2 < N2; i2++)
	{
		int atom1 = i1 + id1;
		int atom2 = i2 + id2;

		float dist[3] = {X[atom1] - X[atom2] + i*L[0][0],
						 Y[atom1] - Y[atom2] + i*L[1][0] + j*L[1][1],
						 Z[atom1] - Z[atom2] + i*L[1][0] + j*L[1][1]+ k*L[2][2]
						};

		if (NORMSQR(dist) < *(vdw_matrix) )
		{
			int i = 0 ;
		}

	}

	return 0;
}

check_lower_triangular(float L[3][3])
{
	if (L[0][1] < CONSTANT_TOL &&
	    L[0][2] < CONSTANT_TOL &&
	    L[1][2] < CONSTANT_TOL
	    )
	    return 1;

	return 0;
}





float find_mol_len(float *X, float *Y, float *Z, int len)
{
	float com[3] = {0, 0, 0};
	find_mol_com(X, Y, Z, len, com);
	printf("com = %f %f %f\n", com[0], com[1], com[2]);
	float first_atom[3] = {X[0], Y[0], Z[0]};
	float max = NORM2(com, first_atom);
	float dist_sqr = 0;

	for(int i = 0; i < len; i++)
	{
		float atom_vec[3] = {X[i], Y[i], Z[i]};
		dist_sqr = NORM2(atom_vec, com);

		if (max < dist_sqr)
			max = dist_sqr;
	}

	return 2*sqrt(max);
}



void find_mol_com(float *X, float *Y, float *Z, int len, float com[3])
{
	com[0] = 0;
	com[1] = 0;
	com[2] = 0;

	for(int i = 0; i < len; i++)
	{
		com[0] += X[i];
		com[1] += Y[i];
		com[2] += Z[i];
	}

	com[0] /= len;
	com[1] /= len;
	com[2] /= len;
}
