#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <stddef.h>
#include "mpi.h"
#include "read_input.h"
#include "spg_generation.h"
#include "combinatorics.h"
#include "check_structure.h"
#include "crystal_utils.h"
#include "molecule_utils.h"
#include "lattice_generator.h"
#include "randomgen.h"
#include "algebra.h"
#include "cgenarris_mpi.h"
#include "pygenarris_mpi.h"

//maximum mulipicity possible
#define ZMAX 192
#define VOL_ATTEMPT  100000
#define GRAIN_SIZE 10000

int *seed;
unsigned int *seed2;
extern float TOL;

void create_vdw_matrix_from_sr(molecule *mol,
								float *vdw_matrix,
								float sr,
								int Z);

int main(int argc, char **argv)
{
	//Initialise MPI
	MPI_Init(&argc, &argv);
	int total_ranks;
    MPI_Comm_size(MPI_COMM_WORLD, &total_ranks);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm world_comm = MPI_COMM_WORLD;

    //variable declarartion
	molecule *mol = (molecule*)malloc(sizeof(molecule));//store molecule

	//added here//////
	int crystal_generation //indicate whether its molecular crystal generation or layer generation
	float interface_area_mean; //for layer
	float interface_area_std; //for layer
	int volume_multiplier; // for layer
	float  lattice_vector_2d[2][3]; //for layer
	///// done/////



	float volume_std;	//standard dev for volumes
	float volume_mean;	//mean volume
	float sr;			//specific radius proportion for structure checks
						//see paper for definition
	float Zp_max;		//Z'' . not implemented
	int Z;				//multiplicity of general position
	int num_structures;	//num of structures per spg
	long max_attempts;	//max attempts per space group
    float tol;
    char spg_dist_type[10];  //spg distribution type
    int vol_attempt;    // no of attempts after which volume is resampled
    int random_seed;    //seed for random gen

    read_geometry(mol);				//read molecule from geometry.in

	if (crystal_generation)     // for molecular crystal generation
	{
		read_control(&num_structures,
				 &Z,
				 &Zp_max,
				 &volume_mean,
				 &volume_std,
				 &sr,
				 &max_attempts,
                 spg_dist_type,
                 &vol_attempt,
                 &random_seed);	//get settings
		tol = TOL;

		int num_atoms_in_molecule = mol->num_of_atoms;
		int dim_vdw_matrix = num_atoms_in_molecule * Z ;
		float *vdw_cutoff_matrix = (float *) malloc( dim_vdw_matrix *
								dim_vdw_matrix *
								sizeof(float) ); //square matrix

		create_vdw_matrix_from_sr(mol, vdw_cutoff_matrix, sr, Z);

		//call the generator from pygenarris_mpi
		mpi_generate_molecular_crystals_with_vdw_cutoff_matrix(
		vdw_cutoff_matrix,
		dim_vdw_matrix,
		dim_vdw_matrix,
		num_structures,
		Z,
		volume_mean,
		volume_std,
		tol,
		max_attempts,
		spg_dist_type,
		vol_attempt,
		random_seed,
		world_comm);

    	MPI_Finalize();

    	return 0;
	}
	else						// for layer generation
	{
		read_control_layer(&num_structures,
				 &Z,
				 &Zp_max,
				 &volume_mean, 
				 &volume_std,
				 &interface_area_mean,
				 &interface_area_std,
				 &volume_multiplier,
				 &sr,
				 &max_attempts,
                 spg_dist_type, lattice_vector_2d);	//get settings
		tol = TOL;
		int num_atoms_in_molecule = mol->num_of_atoms;
		int dim_vdw_matrix = num_atoms_in_molecule * Z ;
		float *vdw_cutoff_matrix = (float *) malloc( dim_vdw_matrix * 
								dim_vdw_matrix *
								sizeof(float) ); //square matrix
	
		create_vdw_matrix_from_sr(mol, vdw_cutoff_matrix, sr, Z);

		mpi_generate_molecular_crystals_with_vdw_cutoff_matrix(
		vdw_cutoff_matrix,
		dim_vdw_matrix,
		dim_vdw_matrix,
		num_structures,
		Z,
		volume_mean,
		volume_std,
		interface_area_mean,
		interface_area_std,
		volume_multiplier,
		tol, 
		max_attempts,
		spg_dist_type,
		lattice_vector_2d,
		world_comm);
	
    	MPI_Finalize();

    	return 0;	
	
	}
    

	
}

