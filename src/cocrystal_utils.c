#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cocrystal_utils.h"
#include "check_structure.h"


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
    cxtal->n_mols = n_mol_asym * Z;
    cxtal_allocate(cxtal, n_atoms);

    // store details to cxtal
    cxtal->Z = Z;
    memcpy(cxtal->stoic, stoic, tbytes);
    memcpy(cxtal->n_atoms_in_mol, n_atoms_in_mol, tbytes);
    cxtal->n_mol_types = n_mol_types;
    cxtal->n_atoms = n_atoms;
    cxtal->spg = 0;
    cxtal->Zp = 0;

    // Get molecule index
    int at = 0;
    int mol_id = 0;
    for(int z = 0; z < Z; z++)
    {
        for(int m = 0; m < n_mol_types; m++)
        {
            for(int st = 0; st < cxtal->stoic[m]; st++)
            {
                cxtal->mol_types[mol_id] = m;
                cxtal->mol_index[mol_id] = at;
                at += n_atoms_in_mol[m];
                mol_id++;
            }
        }
    }
}


void cxtal_allocate(cocrystal *cxtal, int total_atoms)
{
    int fbytes = total_atoms * sizeof(float);
    cxtal->Xcord = malloc(fbytes);
    cxtal->Ycord = malloc(fbytes);
    cxtal->Zcord = malloc(fbytes);
    cxtal->com   = malloc(3 * cxtal->n_mols * sizeof(float));
    cxtal->atoms = malloc(total_atoms * sizeof(char) *2);
}


void cxtal_print(cocrystal *cxtal, FILE* out, int fractional)
{
    fprintf(out, "#Number of atoms = %d\n", cxtal->n_atoms);
    fprintf(out, "#Number of molecules = %d\n", cxtal->n_mols);
    fprintf(out, "#Number of molecule types = %d\n", cxtal->n_mol_types);
    fprintf(out, "#Z = %d\n", cxtal->Z);

    // Print mol index
    fprintf(out, "#mol_index = ");
    for(int i = 0; i < cxtal->n_mols; i++)
        fprintf(out, "%d ", cxtal->mol_index[i]);
    fprintf(out, "\n#mol_types = ");
    for(int i = 0; i < cxtal->n_mols; i++)
        fprintf(out, "%d ", cxtal->mol_types[i]);
    fprintf(out, "\n");

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


void cxtal_compute_COM(cocrystal *cxtal)
{
    for(int mol = 0; mol < cxtal->n_mols; mol++)
    {
        int n_atom_mol = cxtal->n_atoms_in_mol[cxtal->mol_types[mol]];
        float com[3] = {0.0, 0.0, 0.0};

        int start_id = cxtal->mol_index[mol];
        int end_id = cxtal->mol_index[mol] + n_atom_mol;
        for(int mat = start_id; mat < end_id; mat++)
        {
            com[0] += cxtal->Xcord[mat];
            com[1] += cxtal->Ycord[mat];
            com[2] += cxtal->Zcord[mat];
        }
        com[0] /= n_atom_mol;
        com[1] /= n_atom_mol;
        com[2] /= n_atom_mol;

        memcpy(cxtal->com[mol], com, sizeof(float)*3);
    }
}


void cxtal_bring_molecules_first_cell(cocrystal *cxtal)
{
    float inv_lat_vec[3][3];
    float trans_lat_vec[3][3];
    mat3b3_transpose(trans_lat_vec, cxtal->lattice_vectors);
    inverse_mat3b3(inv_lat_vec, trans_lat_vec);

    cxtal_compute_COM(cxtal);

    for(int mol = 0; mol < cxtal->n_mols; mol++)
    {
        // Compute COM shift - integer part of com in fractional cords
        float com_frac[3];
        memcpy(com_frac, cxtal->com[mol], 3*sizeof(float));
        vector3_mat3b3_multiply(inv_lat_vec, com_frac);
        vector3_int(com_frac);
        vector3_mat3b3_multiply(trans_lat_vec, com_frac);

        int n_atom_mol = cxtal->n_atoms_in_mol[cxtal->mol_types[mol]];
        int start_id = cxtal->mol_index[mol];
        int end_id = cxtal->mol_index[mol] + n_atom_mol;
        for(int mat = start_id; mat < end_id; mat++)
        {
            cxtal->Xcord[mat] -= com_frac[0];
            cxtal->Ycord[mat] -= com_frac[1];
            cxtal->Zcord[mat] -= com_frac[2];
        }
    }

    cxtal_compute_COM(cxtal);
}


int cxtal_check_structure(cocrystal *cxtal, Settings *set)
{
    crystal xtal;
    memcpy(xtal.lattice_vectors, cxtal->lattice_vectors, 9*sizeof(float));
    xtal.Xcord = cxtal->Xcord;
    xtal.Ycord = cxtal->Ycord;
    xtal.Zcord = cxtal->Zcord;
    xtal.atoms = cxtal->atoms;

    return structure_checker(&xtal, set->vdw_matrix,
        cxtal->n_atoms, cxtal->mol_index, cxtal->n_mols);

}