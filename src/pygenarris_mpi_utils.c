#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "molecule_utils.h"
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

