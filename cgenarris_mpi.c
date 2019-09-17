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

//maximum mulipicity possible
#define ZMAX 192
#define VOL_ATTEMPT  100000
#define GRAIN_SIZE 10000

int *seed;
unsigned int *seed2;


//void create_mpi_xtal_type(MPI_Datatype* XTAL_TYPE, int total_atoms);
void send_xtal(MPI_Comm comm, int destination, crystal* xtal, int total_atoms);
void receive_xtal(MPI_Comm comm, int source, crystal* xtal, int total_atoms);

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	int total_ranks;
    MPI_Comm_size(MPI_COMM_WORLD, &total_ranks);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Status *status;
    printf("Rank number %d out of %d reporting \n", my_rank+1, total_ranks);

	//random number seeding
	srand((unsigned int)time(NULL));
	int seed_shift = rand()% 1061 + 7 ;
	//variable declarations
	int stop_flag = 0;	// to stop threads if limit is reached
	int counter = 0;	//counts number of structures
	int spg_index = 0;	//space group to be generated
	FILE *out_file;		//file to output geometries
	if (my_rank == 0)
	{
		out_file = fopen("geometry.out","w");
		if(!out_file)		//check permissions
		{
			printf("***ERROR: cannot create geometry.out \n");
			exit(0);
		}
		//fprintf(out_file, "my_rank=%d\n", my_rank);
	}


	//random number seeding, different seeds for different threads
	seed = (int*)malloc(sizeof(int)); //seed for uniform gen
	seed2 = (unsigned int*)malloc(sizeof(unsigned int)); //seed for random
	*seed += my_rank*7 + seed_shift*13; //some random seed private for each threads
	*seed2 = my_rank*17 + seed_shift*11;
	
	//storing information for compatible space groups
	COMPATIBLE_SPG compatible_spg[230]; 
	int num_compatible_spg = 0;
	int num_axes;	//storing number molecule axes
	float *mol_axes; //storing possible molecular axes
	
	//variable declarartion	
	molecule *mol = (molecule*)malloc(sizeof(molecule));//store molecule
	crystal *random_crystal = (crystal*)malloc(sizeof(crystal));//dummy crystal
	float volume_std;	//standard dev for volumes
	float volume_mean;	//mean volume
	float sr;			//specific radius proportion for structure checks
						//see paper for definition
	float Zp_max;		//Z'' . not implemented
	float volume;		//random volume used of generation
	int Z;				//multiplicity of general position
	int num_structures;	//num of structures
	int spg;			//space group attempted
	int max_attempts;	//max attempts per space group

	//read input from file, read molecular structure from geometry,Z, Zp

	if(my_rank == 0)
	{
		int len = 100;
		char name[len];
		gethostname(name, len);
		printf("PARALLELIZATION INFO:\n");
		printf("---------------------------\n");
		printf("Using MPI for parallelization.\n");
		printf("cgenarris is running on %d ranks \n", total_ranks);
		printf("---------------------------\n");
		printf("\n");
	}



    read_geometry(mol);				//read molecule from geometry.in
    read_control(&num_structures,
				 &Z,
				 &Zp_max,
				 &volume_mean, 
				 &volume_std,
				 &sr,
				 &max_attempts);	//get settings

	
	//recenter molecule to origin
	recenter_molecule(mol);

	if (my_rank == 0)
	{
		print_input_geometry(mol);
		print_input_settings(&num_structures,
							 &Z,
							 &Zp_max,
							 &volume_mean, 
							 &volume_std,
							 &sr,
							 &max_attempts);
	}

	MPI_Barrier(MPI_COMM_WORLD);

    //inititalise volume
    do {volume = normal_dist_ab(volume_mean, volume_std);} while(volume < 0.1);
	int N = mol->num_of_atoms;
	allocate_xtal(random_crystal, ZMAX, N); //allcate memory
	random_crystal->num_atoms_in_molecule = mol->num_of_atoms;
	int total_atoms = mol->num_of_atoms * Z;

	//create xtal type structure for MPI
	//MPI_Datatype XTAL_TYPE;
	//create_mpi_xtal_type(&XTAL_TYPE, mol->num_of_atoms*Z);
	//MPI_Type_commit(&XTAL_TYPE);

	
	if(my_rank == 0)
	{
		printf("COMPATIBLE SPACE GROUP INFO:\n");
		printf("-----------------------------------------------------\n");
		printf("Detecting compatible space groups"
			   " from molecule symmetries...\n");
	}

	//compatible space groups are stored in compatible_spg. see the 
	//object definition in spg_generation.h
	//every thread has its own copy
	find_compatible_spg_positions(mol,
								  Z,
								  compatible_spg,
								  &num_compatible_spg,
								  &mol_axes,
								  &num_axes,
								  my_rank);
			
	if(my_rank == 0)
	{
		printf("\n");
		printf("Total compatible space groups = %d\n", num_compatible_spg);
		printf("Number of molecular axes = %d \n", num_axes);
		printf("-----------------------------------------------------\n\n");
		printf("Starting generation...\n\n");
		sleep(1);
	}
	
	MPI_Barrier(MPI_COMM_WORLD); // wait for other friends to join
	
	/*deprecated
	//find allowed space groups for general position (deprecated)
	//int allowed_spg[230];
	//int num_allowed_spg = 0;
	//find_allowed_spg(allowed_spg, &num_allowed_spg, Z);
	//printf("allowed %d \n", num_compatible_spg);
	*/
	
	
	while( spg_index < num_compatible_spg )
	{
		//counter counts the number of structures generated
		spg = compatible_spg[spg_index].spg; //pick a spg 

		//print attempted space group 
		if (my_rank == 0)
			{printf("Attempting to generate spacegroup number %d....\n", spg);}

		while( counter < num_structures )
		{
			int verdict = 0; //for structure check
			int i = 0; 		 //counts attempts for spg
			//attempts for an spg. Assume all threads run equally fast
			for(; i < max_attempts/total_ranks; i = i + GRAIN_SIZE) 
			{
				for(int j = 0; j < GRAIN_SIZE; j++)
				{
						

				
					// else generate again
					int result = generate_crystal(random_crystal,
									 mol,
									 volume,
									 Z,
									 Zp_max,
									 spg,
									 compatible_spg,
									 num_compatible_spg,
									 spg_index,
									 mol_axes,
									 num_axes);
					
					//alignment failure
					if(!result)
						continue;
					
					//check if molecules are too close with sr	    
					verdict = check_structure(*random_crystal, sr);    
					
					//reset volume after volume attempts
					if( (i+j) % VOL_ATTEMPT == 0 && i+j != 0)
					{
						do {volume = normal_dist_ab(volume_mean, volume_std);} while(volume < 0.1);
						if(my_rank == 0)
							printf("#Rank 1: Completed %d attempts\n", (i+j)*total_ranks);
						*seed = *seed2 + my_rank*rand_r(seed2);
						fflush(stdout);
					}

										//if generation is successful
					if (verdict == 1 )
					{
						random_crystal->spg = spg;
						//print_crystal(random_crystal);
						//print_crystal2file(random_crystal, out_file);
						//fprintf(out_file, "my_rank = %d\n", my_rank);
						
						counter++; //one less structure to generate
						
						printf("#Rank %d:Generation successful after %d attempts.\n",
								my_rank,
								i+j);
						/*			
						printf("#Structure number = %d,\n"
							   "#attempted space group = %d,\n"
							   "#unit cell volume (cubic Angstrom)= %f.\n", 
							   counter, 
							   spg,
							   volume);
						fflush(stdout);
						int spglib_spg = detect_spg_using_spglib(random_crystal);
						printf("#SPGLIB detected space group = %d\n\n",
												                     spglib_spg);
						*/
						do {volume = normal_dist_ab(volume_mean, volume_std);} while(volume < 0.1);
						break;
						//verdict = 1;

					}


				}//end of GRAIN loop

				int found_poll[total_ranks];
				//printf("verdict = %d\n", verdict);
				MPI_Gather(&verdict, 1, MPI_INT, &found_poll, 1, MPI_INT, 0, MPI_COMM_WORLD);

				MPI_Barrier(MPI_COMM_WORLD);
				if (my_rank == 0)
				{
					//print the structure generated by root first to outfile
					if(verdict)
						print_crystal2file(random_crystal, out_file);

					//get and print structures from other ranks
					for(int i = 1; i < total_ranks; i++)
					{
						if (found_poll[i] == 1)
						{
							receive_xtal(MPI_COMM_WORLD, i, random_crystal, total_atoms);
							//print_crystal(random_crystal);
							counter++;
							print_crystal2file(random_crystal, out_file);
						}
					}
				}
				else
				{
					
					if (verdict == 1)
					{
						send_xtal(MPI_COMM_WORLD, 0, random_crystal, total_atoms);
					}
					
					
				}
				MPI_Barrier(MPI_COMM_WORLD);

			}//end of attempt loop	
			
			//if max limit is reached or if some thread hit the limit
			if (i >= max_attempts/total_ranks || stop_flag == 1)
			{	
				//stop other threads
				#pragma omp critical
					{stop_flag = 1;}
				#pragma omp barrier
				{}
				do {volume = normal_dist_ab(volume_mean, volume_std);} while(volume < 0.1);
				if(my_rank== 0)
				{	
					printf("**WARNING: generation failed for space group = %d "
							"after %d attempts. \n",
							spg,
							max_attempts);
					counter = num_structures + 1;
					fflush(stdout);
					//print_crystal(random_crystal);
				}
				MPI_Barrier(MPI_COMM_WORLD);
			}	
			
		}//end of numof structures whileloop
	
		//move to next spacegroup
		MPI_Barrier(MPI_COMM_WORLD);
		if(my_rank == 0 )
		{	
			 counter = 0;
			 spg_index++;
			 stop_flag = 0;
			 printf("#space group counter reset. Moving to next space group...\n\n");
			 fflush(stdout);
		} 
		MPI_Barrier(MPI_COMM_WORLD);
		
	}//end of spg while loop

																	    //omp block ends here.
	fclose(out_file);
	//MPI_Type_free(&XTAL_TYPE);
	MPI_Finalize();
	printf("Generation completed!\n");

}


void send_xtal(MPI_Comm comm, int destination, crystal* xtal, int total_atoms)
{
	MPI_Send(xtal->lattice_vectors, 9, MPI_FLOAT , destination, 1, comm);
	MPI_Send(xtal->Xcord, total_atoms, MPI_FLOAT , destination, 2, comm);
	MPI_Send(xtal->Ycord, total_atoms, MPI_FLOAT , destination, 3, comm);
	MPI_Send(xtal->Zcord, total_atoms, MPI_FLOAT , destination, 4, comm);
	MPI_Send(xtal->atoms, 2*total_atoms, MPI_CHAR , destination, 5, comm);
	MPI_Send(&(xtal->spg), 1, MPI_INT , destination, 6, comm);
	MPI_Send(&(xtal->wyckoff_position), 1, MPI_INT , destination, 7, comm);
	MPI_Send(&(xtal->num_atoms_in_molecule), 1, MPI_INT , destination, 8, comm);
	MPI_Send(&(xtal->Z), 1, MPI_INT , destination, 9, comm);
	MPI_Send(&(xtal->Zp), 1, MPI_INT , destination, 10, comm);
}

void receive_xtal(MPI_Comm comm, int source, crystal* xtal, int total_atoms)
{
	MPI_Status *status;
	MPI_Recv(xtal->lattice_vectors, 9, MPI_FLOAT, source, 1, comm, status);
	MPI_Recv(xtal->Xcord, total_atoms, MPI_FLOAT, source, 2, comm, status);
	MPI_Recv(xtal->Ycord, total_atoms, MPI_FLOAT, source, 3, comm, status);
	MPI_Recv(xtal->Zcord, total_atoms, MPI_FLOAT, source, 4, comm, status);
	MPI_Recv(xtal->atoms, 2*total_atoms, MPI_CHAR, source, 5, comm, status);
	MPI_Recv(&(xtal->spg), 1, MPI_INT, source, 6, comm, status);
	MPI_Recv(&(xtal->wyckoff_position), 1, MPI_INT, source, 7, comm, status);
	MPI_Recv(&(xtal->num_atoms_in_molecule), 1, MPI_INT, source, 8, comm, status);
	MPI_Recv(&(xtal->Z), 1, MPI_INT, source, 9, comm, status);
	MPI_Recv(&(xtal->Zp), 1, MPI_INT, source, 10, comm, status);


}

/*
void create_mpi_xtal_type(MPI_Datatype* XTAL_TYPE, int total_atoms)
{
	const int struct_len = 10;
	int block_lengths[ ] = {9, total_atoms, total_atoms, total_atoms, 
						   2*total_atoms, 1, 1, 1, 1, 1};
	MPI_Datatype types[ ] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,
							 MPI_CHAR, MPI_INT, MPI_INT, MPI_INT, MPI_INT,
							 MPI_INT};

	//Create a fake crystal to get the memory displacements
	crystal xtal;
	allocate_xtal(&xtal, total_atoms, 1);

	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			xtal.lattice_vectors[i][j] = 1.0;

	for(int i = 0; i < total_atoms; i++)
	{
		xtal.Xcord[i] = 1.0;
		xtal.Ycord[i] = 1.0;
		xtal.Zcord[i] = 1.0;
		xtal.atoms[2*i] ='C';
		xtal.atoms[2*i + 1] ='C';
	}
	xtal.spg = 1;
	xtal.wyckoff_position = 1;
	xtal.num_atoms_in_molecule = 1;
	xtal.Z = total_atoms;
	xtal.Zp = 1;

	//print_crystal(&xtal);
	//Now find the displacements
	MPI_Aint displacements[struct_len];
	MPI_Get_address(&xtal.lattice_vectors, &displacements[0]);
	MPI_Get_address(&xtal.Xcord, &displacements[1]);
	MPI_Get_address(&xtal.Ycord, &displacements[2]);
	MPI_Get_address(&xtal.Zcord, &displacements[3]);
	MPI_Get_address(&xtal.atoms, &displacements[4]);
	MPI_Get_address(&xtal.spg, &displacements[5]);
	MPI_Get_address(&xtal.wyckoff_position, &displacements[6]);
	MPI_Get_address(&xtal.num_atoms_in_molecule, &displacements[7]);
	MPI_Get_address(&xtal.Z, &displacements[8]);
	MPI_Get_address(&xtal.Zp, &displacements[9]);

	//make it relative to the first
	for(int i = 0; i < struct_len; i++)
	{	displacements[i] -= displacements[0];
		printf("%d\n", displacements[i] );
	}

	//create and commit the new MPI datatype
	MPI_Type_create_struct(struct_len, block_lengths, displacements,
						   types, XTAL_TYPE);
    
	free(xtal.atoms);
	free(xtal.Xcord);
	free(xtal.Ycord);
	free(xtal.Zcord);

}
*/

