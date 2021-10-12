// test example for the rigid-body optimizer

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../../crystal.h"
#include "../../crystal_utils.h"
#include "../../cocrystal.h"
#include "../../spglib.h"
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


void run_example(char *xdir, int natoms_per_mol, int  Z, int cell_type, int spg)
{
    crystal xtl;
    float *cutmat;

    printf("---------------------------------\n");
    printf("Running example: %s\n", xdir);
    time_t start = time(NULL);

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
    set.max_iteration = 4000;
    set.spg = spg;
    int placeholder;
    Opt_status status = optimize_crystal(&xtl, cutmat, placeholder, set);

    time_t end = time(NULL);
    if(status == SUCCESS)
    { print_crystal(&xtl); }
    else if(status  == ITER_LIMIT)
    { printf("Optimization failed: Max iterations reached\n"); }
    else
    { printf("Optimization failed\n"); }

    int spg_detect = detect_spg_using_spglib(&xtl);
    printf("spglib detected spacegroup = %d\n", spg_detect);

    printf("Completed optmization in %.2f seconds\n", difftime(end, start));
    printf("---------------------------------\n\n");
    free(xtl.Xcord);
    free(xtl.Ycord);
    free(xtl.Zcord);
    free(xtl.atoms);
    free(cutmat);

}


int main(void)
{
    /*
    run_example("sample_structures/Example1", 12, 2, TRICLINIC, 2);
    run_example("sample_structures/fast_opt", 30, 2, MONOCLINIC, 12);
    run_example("sample_structures/failed_1", 30, 2, MONOCLINIC, 6);
    run_example("sample_structures/Example2", 12, 4, MONOCLINIC, 13);
    run_example("sample_structures/Example5", 30, 4, MONOCLINIC, 9);
    run_example("sample_structures/Example4", 30, 2, TRICLINIC, 2);
    */
    //run_example("../sample_structures/Example5", 30, 4, MONOCLINIC, 9);
    //run_example("../sample_structures/Example1", 12, 2, TRICLINIC, 2);
    //run_example("../sample_structures/Example2", 12, 4, MONOCLINIC, 13);
    //run_example("../sample_structures/Example2", 12, 4, MONOCLINIC, 13);
    run_example("../sample_structures/Example6", 12, 4, TRICLINIC,0);
    //run_example("../sample_structures/Example2", 12, 4, MONOCLINIC, 13);
    //run_example("../sample_structures/failed_2", 30, 2, TRICLINIC, 0);
    return 0;
}


