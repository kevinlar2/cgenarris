#ifndef _PYGENARRIS_MPI_UTILS_H_
#define _PYGENARRIS_MPI_UTILS_H_

void print_time(void);
FILE* open_output_file(int my_rank);
void init_random_seed(unsigned int *seed, unsigned int *seed2,
                      int random_seed, int rank);


#endif //  pygenarris_mpi_utils.h