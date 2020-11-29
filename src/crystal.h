#ifndef _CRYSTAL_H
#define _CRYSTAL_H

typedef struct
{
	float lattice_vectors[3][3];
	float *Xcord;
	float *Ycord;
	float *Zcord;
	char *atoms;
	int spg;
	int wyckoff_position;
    int num_atoms_in_molecule;
	int Z;
	int Zp;
	
}crystal;

#endif
