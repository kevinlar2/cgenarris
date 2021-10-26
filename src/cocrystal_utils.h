#ifndef _COCRYSTAL_UTILS_H_
#define _COCRYSTAL_UTILS_H_

#include "cocrystal.h"
#include "input_settings.h"

void cxtal_init(cocrystal *cxtal, int *stoic, int *n_atoms_in_mol, int n_mol_types, int Z);
void cxtal_allocate(cocrystal *cxtal, int total_atoms);
void cxtal_print(cocrystal *cxtal, FILE* out, int fractional);
int cxtal_check_structure(cocrystal *cxtal, Settings *set);
float cxtal_get_cell_volume(cocrystal *cxtal);
int cxtal_check_structure(cocrystal *cxtal, Settings *set);
void cxtal_bring_molecules_first_cell(cocrystal *cxtal);
void cxtal_create_from_data(cocrystal *cxtal, double lattice[3][3],
			    double *coords, int total_atom1,
			    char *atoms, int total_atoms2,
			    int Z, int *stoic, int n_types1,
			    int *n_atoms_in_mol, int n_types2, int spg);

void cxtal_get_data(cocrystal *cxtal, double lattice[3][3], double *coords,
		    int total_atom1);

void cxtal_free(cocrystal *cxtal);
#endif
