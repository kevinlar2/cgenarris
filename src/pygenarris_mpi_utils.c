#include <stdio.h>
#include <stdlib.h>
#include <time.h>
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