#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "molecule_utils.h"
#include "pygenarris_mpi_utils.h"
#include "randomgen.h"

void print_time(void)
{
    time_t current_time = time(NULL);
    struct tm *loc_time;
    loc_time = localtime (&current_time);
    printf("The local time is now: %s", asctime (loc_time));
}

FILE* open_output_file(int my_rank)
{
    FILE *out_file = NULL;

    if (my_rank == 0)
    {
        out_file = fopen("geometry.out", "w");
        if(!out_file)       //check permissions
        {
            printf("***ERROR: cannot create geometry.out \n");
            exit(EXIT_FAILURE);
        }
    }

    return out_file;
}

void init_random_seed(unsigned int *seed, unsigned int *seed2, int random_seed, int rank)
{
    //random number seeding, different seeds for different threads
    if (random_seed == 0)
    {
        srand((unsigned int)time(NULL));
        random_seed = rand();
    }
    else
    {
        srand((unsigned int) 19023411);
    }

    *seed = (unsigned int)abs(rank*7 + random_seed);  //some random seed private for each threads
    *seed2 = (unsigned int)abs(rank*17 + random_seed);
    init_genrand(*seed);

}

void recenter_molecules(molecule* mol, int mol_types)
{
    for(int i = 0; i < mol_types; i++)
    {
        recenter_molecule(mol + i);
    }
}

float draw_volume(float volume_mean, float volume_std)
{
    float volume;

    do
    {
        volume = normal_dist_ab(volume_mean, volume_std);
    }while(volume < 1);

    return volume;
}

/*
Finds the number of atoms in each molecule from molecule array
and constructs n_atoms_in_mol array of length n_mol_types
*/
void get_n_atoms_in_mol(int *n_atoms_in_mol, molecule *mol, int n_mol_types)
{
    for(int m = 0; m < n_mol_types; m++)
    {
        n_atoms_in_mol[m] = (mol + m)->num_of_atoms;
    }
}


void print_allowed_spg(int *allowed_spg, int num_spg)
{
    printf("Detecting compatible spacegroups using only general Wyckoff positions:\n");
    printf("Number of allowed spg: %d\n", num_spg);
    for(int spg_id = 0; spg_id < num_spg; spg_id++)
    {
        int spg = allowed_spg[spg_id];
        printf("Spacegroup %d is compatible.\n", spg);
    }
    printf("\n");
}

/*
main structure generation loop
*/
int try_crystal_generation(cocrystal cxtal,
                           Settings set,
                           float *volume,
                           float *vdw_matrix,
                           long attempts,
                           long batch_size,
                           int spg)
{
    cxtal.spg = spg;
    // Loop over the batch. bat = batch attempt
    for(long bat = 0; bat < batch_size; bat++)
    {
        // Attempt one generation
        int result = generate_cocrystal(cxtal, set, *volume, vdw_matrix);

        // Alignment failures
        if(!result)
            continue;

        // reset volume
        if( (attempts + bat) % set.vol_attempts == 0)
        {
            *volume = draw_volume(set.vol_mean, set.vol_std);
            fflush(stdout);
        }

        // Structure check
        int verdict = cxtal_check_structure(cxtal, vdw_matrix);
        if(verdict)
            break;
    }

}


/*
Checks if structure generation can be stopped.
Check if enough structures were generated. OR
ran out of attempts
*/
int check_stop_condition()
{
    return 1;
}


/*
Slave processes send structures to master rank
*/
int send_structures()
{
    return 1;
}

/*
Master rank write to output file
*/
int write_structures()
{
    return 1;
}

/*
Print after Genarris exits
*/
void print_exit()
{
    printf("Genarris has completed generation.\n");
}