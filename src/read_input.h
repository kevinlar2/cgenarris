#ifndef READ_INPUT_H
#define READ_INPUT_H

#include <stdio.h>
#include <stdlib.h>

#include "molecule.h"

void read_control(int* num_structures,
                  int* Z,
                  float* Zp_max,
                  float* volume_mean,
                  float* volume_std,
                  float *sr,
                  long *max_attempts,
                  char *spg_dist_type, 
                  int *vol_attempt, 
                  int *random_seed);

void read_geometry(molecule* mol);

void print_input_geometry(molecule* mol);

void print_molecule(molecule *mol);

void print_input_settings(int* num_structures,
                          int* Z,
                          float* Zp_max,
                          float* volume_mean,
                          float* volume_std,
                          float *sr,
                          long *max_attempts,
                          char * spg_dist_type, 
                          int *vol_attempt, 
                          int *random_seed);

#endif 
