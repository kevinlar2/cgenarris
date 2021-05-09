#ifndef _PYGENARRIS_MPI_UTILS_H_
#define _PYGENARRIS_MPI_UTILS_H_

#include "input_settings.h"
#include "cocrystal.h"
#include "molecule.h"

void print_time(void);
FILE* open_output_file(int my_rank);
void init_random_seed(unsigned int *seed, unsigned int *seed2,
                      int random_seed, int rank);
void recenter_molecules(molecule* mol, int mol_types);
float draw_volume(float volume_mean, float volume_std);
void get_n_atoms_in_mol(int *n_atoms_in_mol, molecule *mol, int n_mol_types);
void print_allowed_spg(int *allowed_spg, int num_spg);
int try_crystal_generation(cocrystal *cxtal,
                           Settings set,
                           molecule *mol,
                           float *volume,
                           long attempts,
                           long batch_size);


#endif //  pygenarris_mpi_utils.h