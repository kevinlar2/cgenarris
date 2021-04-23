#ifndef _COCRYSTAL_H_
#define _COCRYSTAL_H_

typedef struct
{
    float lattice_vectors[3][3];
    float *Xcord;
    float *Ycord;
    float *Zcord;
    char  *atoms;
    int   *mol_index;
    int   *mol_types;
    int   *n_atoms_in_mol;
    int   *wyckoff_position;
    int   *stoic;
    int   n_mols;
    int   n_mol_types;
    int   n_atoms;
    int   spg;
    int   Z;
    int   Zp;
}cocrystal;


#endif  // Cocrystal.h
