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

void read_cocrystal(char *path, cocrystal *xtl)
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

void cxtal_print(cocrystal *cxtal, FILE* out, int fractional)
{
    static int struct_num = 1;

    fprintf(out, "####### BEGIN STRUCTURE #######\n");
    fprintf(out, "#structure_number = %d\n", struct_num);
    struct_num++;

    fprintf(out, "#Number of atoms = %d\n", cxtal->n_atoms);
    fprintf(out, "#Number of molecules = %d\n", cxtal->n_mols);
    fprintf(out, "#Number of molecule types = %d\n", cxtal->n_mol_types);
    fprintf(out, "#Z = %d\n", cxtal->Z);
    fprintf(out, "#attempted_spacegroup = %d\n", cxtal->spg);

    // Print mol index
    fprintf(out, "#mol_index = ");
    for(int i = 0; i < cxtal->n_mols; i++)
        fprintf(out, "%d ", cxtal->mol_index[i]);
    fprintf(out, "\n#mol_types = ");
    for(int i = 0; i < cxtal->n_mols; i++)
        fprintf(out, "%d ", cxtal->mol_types[i]);
    fprintf(out, "\n");

    for(int ltt = 0; ltt < 3; ltt++)
    {
        fprintf(out,"lattice_vector %12f %12f %12f \n",
                cxtal->lattice_vectors[ltt][0],
                cxtal->lattice_vectors[ltt][1],
                cxtal->lattice_vectors[ltt][2]);
    }

    for(int at = 0; at < cxtal->n_atoms; at++)
    {
        if(!fractional)
            fprintf(out, "atom %12f %12f %12f  %c%c \n", cxtal->Xcord[at],
                    cxtal->Ycord[at],  cxtal->Zcord[at], cxtal->atoms[2*at],
                    cxtal->atoms[2*at + 1]);
        else
            fprintf(out, "atom_frac %12f %12f %12f  %c%c \n", cxtal->Xcord[at],
                    cxtal->Ycord[at],  cxtal->Zcord[at], cxtal->atoms[2*at],
                    cxtal->atoms[2*at + 1]);

    }

    fprintf(out, "#######  END  STRUCTURE #######\n\n");
    fflush(out);
}


run_cocrystal_example(char *xdir, int natom1, int natom2, int  Z)
{
    int total_atoms = (natom1 + natom2) * Z;
    cocrystal *cxtal = (cocrystal *)malloc(sizeof(cocrystal));
    int fbytes = total_atoms * sizeof(float);
    cxtal->Xcord = malloc(fbytes);
    cxtal->Ycord = malloc(fbytes);
    cxtal->Zcord = malloc(fbytes);
    cxtal->com   = malloc(3 * cxtal->n_mols * sizeof(float));
    cxtal->atoms = malloc(total_atoms * sizeof(char) *2);
    cxtal->n_atoms = total_atoms;
    cxtal->n_mol_types = 2;

    int tbytes = 2 * sizeof(int);  // Total number of mol types
    int mbytes = 2 * Z * sizeof(int);  // Total num of molecules
    cxtal->stoic             = (int *) malloc(tbytes);
    cxtal->n_atoms_in_mol    = (int *) malloc(tbytes);
    cxtal->wyckoff_position  = (int *) malloc(tbytes);
    cxtal->mol_index = (int *) malloc(mbytes);
    cxtal->mol_types = (int *) malloc(mbytes);

    // Assign values
    int stoic[2] = {1, 1};
    memcpy(cxtal->stoic, stoic, 2*sizeof(int));
    cxtal->Z = Z;
    int n_atoms_in_mol[2] = {natom1, natom2};
    memcpy(cxtal->n_atoms_in_mol, n_atoms_in_mol, 2*sizeof(int));
    cxtal->n_mols = 2*Z;
    cxtal->spg = 0;

    // Get molecule index
    int at = 0;
    int mol_id = 0;
    for(int z = 0; z < Z; z++)
    {
        for(int m = 0; m < 2; m++)
        {
            for(int st = 0; st < cxtal->stoic[m]; st++)
            {
		printf("at = %d, mol_id = %d\n", at, mol_id);
                cxtal->mol_types[mol_id] = m;
                cxtal->mol_index[mol_id] = at;
                at += n_atoms_in_mol[m];
                mol_id++;
            }
        }
    }

    // Read cocrystal
    char geo_path[1000];
    strcpy(geo_path, xdir);
    strcat(geo_path, "/geometry.in");
    read_cocrystal(geo_path, cxtal);
    cxtal_print(cxtal, stdout, 0);
    
    char cutmat_path[1000];
    float *cutmat = malloc(sizeof(float)*total_atoms*total_atoms);
    strcpy(cutmat_path, xdir);
    strcat(cutmat_path, "/cutoff_matrix.txt");
    read_vector(total_atoms*total_atoms, cutmat_path, cutmat);
    Opt_settings set;
    set.cell_family = 0;
    set.max_iteration = 400;
    set.spg = 0;
    int placeholder;
    printf("Starting optimization\n");
    Opt_status status = optimize_cocrystal(cxtal, cutmat, set);

    printf("Completed optimization\n");
    time_t end = time(NULL);
    if(status == SUCCESS)
    { printf("\ndone\n"); }
    else if(status  == ITER_LIMIT)
    { printf("Optimization failed: Max iterations reached\n"); }
    else
    { printf("Optimization failed\n"); }

    cxtal_print(cxtal, stdout, 0);

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
    //run_example("../sample_structures/Example5", 30, 4, TRICLINIC, 0);
    run_cocrystal_example("../sample_structures/Example5", 30, 30, 2);
    //run_example("../sample_structures/Example1", 12, 2, TRICLINIC, 0);
    //run_example("../sample_structures/Example2", 12, 4, MONOCLINIC, 13);
    //run_example("../sample_structures/Example2", 12, 4, MONOCLINIC, 13);
    //run_example("../sample_structures/Example6", 12, 4, TRICLINIC,0);
    //run_example("../sample_structures/Example2", 12, 4, MONOCLINIC, 13);
    //run_example("../sample_structures/failed_2", 30, 2, TRICLINIC, 0);
    return 0;
}


