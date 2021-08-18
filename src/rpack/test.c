// test example for the rigid-body optimizer

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../crystal.h"
#include "../crystal_utils.h"
#include "../cocrystal.h"
#include "../spglib.h"
#include "../algebra.h"
#include "rigid_press.h"

unsigned int seed2;

void run_example(char *xdir, int natoms_per_mol,int  Z, int cell_type, int spg);

void read_vector(int size, char *path, float *vec)
{
    FILE *file_ptr = fopen(path, "r");

    for(int i=0 ; i<size ; i++)
    { fscanf(file_ptr, "%f", vec+i); }

    fclose(file_ptr);
}

void read_crystal(char *path, crystal *xtl)
{
    char line[1000], dummy[1000];
    FILE *file_ptr = fopen(path, "r");

    int i=0, j=0;
    while(!feof(file_ptr))
    {
        char *out = fgets(line, 1000, file_ptr);
        if(out == NULL)
        { break; }

        switch(line[0])
        {
            case 'l':
            sscanf(line, "%s %f %f %f", dummy, xtl->lattice_vectors[j], xtl->lattice_vectors[j]+1, xtl->lattice_vectors[j]+2);
            j++;
            break;

            case 'a':
            sscanf(line, "%s %f %f %f %s", dummy, xtl->Xcord+i, xtl->Ycord+i, xtl->Zcord+i, xtl->atoms+2*i);
            i++;
            break;

            case '#':
            break;

            default:
            printf("ERROR: parsing line (%s)\n", line);
            exit(1);
        }
    }
    fclose(file_ptr);
}

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
        printf("atom %12f %12f %12f  %c \n", xtal->Xcord[i],
            xtal->Ycord[i],  xtal->Zcord[i],  xtal->atoms[2*i]);
    }
}


void run_example(char *xdir, int natoms_per_mol, int  Z, int cell_type, int spg)
{
    crystal xtl;
    float *cutmat;

    printf("Running example: %s\n", xdir);

    xtl.num_atoms_in_molecule = natoms_per_mol;
    xtl.Z = Z;
    int num_atoms = xtl.Z*xtl.num_atoms_in_molecule;
    xtl.Xcord = (float*)malloc(sizeof(float)*num_atoms);
    xtl.Ycord = (float*)malloc(sizeof(float)*num_atoms);
    xtl.Zcord = (float*)malloc(sizeof(float)*num_atoms);
    xtl.atoms = (char*)malloc(sizeof(char)*(2*num_atoms+1));
    cutmat = (float*)malloc(sizeof(float)*num_atoms*num_atoms);
    for(int i=0 ; i<=2*num_atoms ; i++)
    { xtl.atoms[i] = ' '; }

    char geo_path[1000];
    strcpy(geo_path, xdir);
    strcat(geo_path, "/geometry.in");
    read_crystal(geo_path, &xtl);
    char cutmat_path[1000];
    strcpy(cutmat_path, xdir);
    strcat(cutmat_path, "/cutoff_matrix.txt");
    read_vector(num_atoms*num_atoms, cutmat_path, cutmat);
    //print_crystal(&xtl);

    Opt_settings set;
    set.cell_family = cell_type;
    set.max_iteration = 400;
    set.spg = spg;
    Opt_status status = optimize_crystal(&xtl, cutmat, set);

    if(status == SUCCESS)
    { print_crystal(&xtl); }
    else if(status  == ITER_LIMIT)
    { printf("Optimization failed: Max iterations reached\n"); }
    else
    { printf("Optimization failed\n"); }

    int spg_detect = detect_spg_using_spglib(&xtl);
    printf("spg = %d\n", spg_detect);
    
    free(xtl.Xcord);
    free(xtl.Ycord);
    free(xtl.Zcord);
    free(xtl.atoms);
    free(cutmat);

}


int main(void)
{
    // run_example("sample_structures/Example1", 12, 2, TRICLINIC, 2);
  //run_example("sample_structures/fast_opt", 30, 2, MONOCLINIC);
  //run_example("sample_structures/failed_1", 30, 2, MONOCLINIC);
    run_example("sample_structures/Example2", 12, 4, MONOCLINIC, 13);
  //run_example("sample_structures/Example5", 30, 4, MONOCLINIC);
  // run_example("sample_structures/Example4", 30, 2, TRICLINIC);
    return 0;
}


int detect_spg_using_spglib(crystal* xtal)
{
    float tol = 1e-3;
    //print_crystal(xtal);
    convert_xtal_to_fractional(xtal);
    //variable declarations
    int num_atoms_in_cell = xtal->Z * xtal-> num_atoms_in_molecule ;
    int types[num_atoms_in_cell];
    double positions[num_atoms_in_cell][3];
    char atom[num_atoms_in_cell*2];
    char symbol[21];
    int num_spg;
    double lattice_vector[3][3];
    //printf("total atoms = %d \n", num_atoms_in_cell);
    for(int i = 0; i < num_atoms_in_cell; i++)
    {
        positions[i][0] = xtal->Xcord[i];
        positions[i][1] = xtal->Ycord[i];
        positions[i][2] = xtal->Zcord[i];
        atom[2*i] = xtal->atoms[2*i];
        atom[2*i+1] = xtal->atoms[2*i+1];

        if      (atom[2*i] == 'C' )
        types[i] = 6;
        else if (atom[2*i] == 'H' )
        types[i] = 1;
        else if (atom[2*i] == 'N' && atom[2*i+1] == ' ')
        types[i] = 7;
        else if (atom[2*i] == 'O' && atom[2*i+1] == ' ')
        types[i] = 8;
        else if (atom[2*i] == 'F' && atom[2*i+1] == ' ')
        types[i] = 9;
        else if (atom[2*i] == 'P' && atom[2*i+1] == ' ')
        types[i] = 15;
        else if (atom[2*i] == 'S' && atom[2*i+1] == ' ')
        types[i] = 16;
        else if (atom[2*i] == 'B' && atom[2*i+1] == 'r')
        types[i] = 35;
        else if (atom[2*i] == 'I' && atom[2*i+1] == ' ')
        types[i] = 53;
        else if (atom[2*i] == 'B' && atom[2*i+1] == ' ')
        types[i] = 5;
        else if (atom[2*i] == 'H' && atom[2*i+1] == 'e')
        types[i] = 2;
        else if (atom[2*i] == 'N' && atom[2*i+1] == 'e')
        types[i] = 10;
        else if (atom[2*i] == 'K' && atom[2*i+1] == 'r')
        types[i] = 36;
        else if (atom[2*i] == 'S' && atom[2*i+1] == 'i')
        types[i] = 14;
        else
        {printf("***ERROR: spglib detector: atom not found -> %c%c\n", atom[2*i], atom[2*i+1]);exit(0);}

    }
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++ )
        {
            //note transpose
            lattice_vector[i][j] = xtal->lattice_vectors[j][i];
        }

    num_spg = spg_get_international(symbol,
                                    lattice_vector,
                                    positions,
                                    types,
                                    num_atoms_in_cell,
                                    tol);
    //SpglibError error;
    //error = spg_get_error_code();
    //printf("#SPGLIB says %s\n", spg_get_error_message(error));
    //printf("#SPGLIB detected space group = %d\n",
    //                                          num_spg);
    convert_xtal_to_cartesian(xtal);
    return num_spg;

}

void convert_xtal_to_cartesian(crystal *xtal)
{
    float lattice_vectors_transpose[3][3];
    copy_mat3b3_mat3b3(lattice_vectors_transpose,
                       xtal->lattice_vectors);
    mat3b3_transpose(lattice_vectors_transpose,
                     lattice_vectors_transpose);
    int total_atoms = xtal->Z *xtal->num_atoms_in_molecule;

    for(int i = 0; i < total_atoms; i++)
    {
        float temp[3] = {xtal->Xcord[i],
                         xtal->Ycord[i],
                         xtal->Zcord[i] };
        vector3_mat3b3_multiply(lattice_vectors_transpose, temp, temp );
        xtal->Xcord[i] = temp[0];
        xtal->Ycord[i] = temp[1];
        xtal->Zcord[i] = temp[2];
    }
}

void convert_xtal_to_fractional(crystal *xtal)
{
    float lattice_vectors[3][3];
    float inv_lattice_vectors[3][3];
    copy_mat3b3_mat3b3(lattice_vectors, xtal->lattice_vectors);
    inverse_mat3b3(inv_lattice_vectors, lattice_vectors);
    mat3b3_transpose(inv_lattice_vectors, inv_lattice_vectors);

    int total_atoms = xtal->Z *xtal->num_atoms_in_molecule;

    for(int i = 0; i < total_atoms; i++)
    {
        float temp[3] = {xtal->Xcord[i], xtal->Ycord[i], xtal->Zcord[i] };
        vector3_mat3b3_multiply(inv_lattice_vectors, temp, temp );
        xtal->Xcord[i] = temp[0];
        xtal->Ycord[i] = temp[1];
        xtal->Zcord[i] = temp[2];
    }
}

