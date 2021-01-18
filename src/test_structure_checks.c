#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "crystal.h"
#include "molecule.h"
#include "molecule_utils.h"
#include "crystal_utils.h"
#include "check_structure.h"
#include "read_input.h"


void read_xtal_from_file(char *fname, crystal* xtal);
void print_vdw_matrix(float *vdw_matrix, int dim_vdw_matrix);

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

void print_vdw_matrix(float *vdw_matrix, int dim_vdw_matrix)
{
    for(int i = 0; i < dim_vdw_matrix; i++)
    {
        for(int j = 0; j < dim_vdw_matrix; j++)
        {
            printf("%f ", *(vdw_matrix + i*dim_vdw_matrix +j));
        }
        printf("\n");
    }

}

int main(int argc, char* argv[])
{
    if(argc != 3)
    {
        printf("sr and filename are required arguments\n");
        exit(1);
    }

    float sr = atof(argv[1]);
    printf("Using sr = %f\n", sr);

    printf("Reading from %s\n", argv[2]);
    crystal xtal;
    read_xtal_from_file(argv[2], &xtal);
    xtal.Z = 2;
    xtal.num_atoms_in_molecule = 100;
    //printf("Z = %d\n", xtal.Z );
    //print_crystal(&xtal);

    int dim_vdw_matrix = xtal.num_atoms_in_molecule * xtal.Z ;
    float *vdw_matrix = (float *) malloc( dim_vdw_matrix *
                                dim_vdw_matrix *
                                sizeof(float) ); //square matrix

    molecule mol;
    read_geometry(&mol);
    //print_molecule(&mol);
    create_vdw_matrix_from_sr(&mol, vdw_matrix, sr, xtal.Z);

    int verdict = check_structure_with_vdw_matrix(xtal, vdw_matrix, dim_vdw_matrix, dim_vdw_matrix);
    printf("Result = %d\n", verdict);
}
