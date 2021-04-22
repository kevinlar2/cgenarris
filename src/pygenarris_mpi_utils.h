#ifndef _PYGENARRIS_MPI_UTILS_H_
#define _PYGENARRIS_MPI_UTILS_H_

void print_time(void);
FILE* open_output_file(int my_rank);
void init_random_seed(unsigned int *seed, unsigned int *seed2,
                      int random_seed, int rank);
void recenter_molecules(molecule* mol, int mol_types);
float draw_volume(float volume_mean, float volume_std);
int find_total_atoms(molecule* mol, int *stoic, int mol_types);


#endif //  pygenarris_mpi_utils.h