#ifndef _COCRYSTAL_UTILS_H_
#define _COCRYSTAL_UTILS_H_

#include "cocrystal.h"

void cxtal_init(cocrystal *cxtal, int *stoic, int *n_atoms_in_mol, int n_mol_types, int Z);
void cxtal_allocate(cocrystal *cxtal, int total_atoms);
void cxtal_print(cocrystal *cxtal, FILE* out);


#endif