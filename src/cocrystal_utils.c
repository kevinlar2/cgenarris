#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cocrystal_utils.h"

void cxtal_init(cocrystal *cxtal, int *stoic, int *n_atoms_in_mol, int n_mol_types, int Z)
{
    int n_atoms = 0;
    int n_mol_asym = 0;  // Number of molecules in asymmetric unit

    // Find number of atoms first
    for(int m = 0; m < n_mol_types; m++)
    {
        // Number of atoms in asym unit
        n_atoms += n_atoms_in_mol[m] * stoic[m];
        n_mol_asym += stoic[m];
    }

    // Total atoms in unit cell
    n_atoms *= Z;

    // Allocate memory
    int tbytes = n_mol_types * sizeof(int);  // Total number of mol types
    int mbytes = n_mol_asym * Z * sizeof(int);  // Total num of molecules
    cxtal->stoic             = (int *) malloc(tbytes);
    cxtal->n_atoms_in_mol    = (int *) malloc(tbytes);
    cxtal->wyckoff_position  = (int *) malloc(tbytes);
    cxtal->mol_index = (int *) malloc(mbytes);
    cxtal->mol_types = (int *) malloc(mbytes);
    cxtal_allocate(cxtal, n_atoms);

    // store details to cxtal
    cxtal->Z = Z;
    memcpy(cxtal->stoic, stoic, tbytes);
    memcpy(cxtal->n_atoms_in_mol, n_atoms_in_mol, tbytes);
    cxtal->n_mol_types = n_mol_types;
    cxtal->n_atoms = n_atoms;
    cxtal->n_mols = n_mol_asym * Z;
    cxtal->spg = 0;
    cxtal->Zp = 0;
}


void cxtal_allocate(cocrystal *cxtal, int total_atoms)
{
    cxtal->Xcord = (float *) malloc(total_atoms * sizeof(float));
    cxtal->Ycord = (float *) malloc(total_atoms * sizeof(float));
    cxtal->Zcord = (float *) malloc(total_atoms * sizeof(float));
    cxtal->atoms = (char *)  malloc(total_atoms * sizeof(char) * 2);
}

void cxtal_print(cocrystal *cxtal, FILE* out, int fractional)
{
    printf("#Number of atoms = %d\n", cxtal->n_atoms);
    printf("#Number of molecules = %d\n", cxtal->n_mols);
    printf("#Number of molecule types = %d\n", cxtal->n_mol_types);
    printf("#Z = %d\n", cxtal->Z);
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
}
