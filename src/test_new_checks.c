#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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


typedef struct
{ 
	char *atoms;
	float *X;
	float *Y;
	float *Z;
	int num_of_atoms;
}molecule;


float vector3_norm(float a[3])
{
	return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

float det_mat3b3(float a[3][3])
{
	 return a[0][0] * ((a[1][1]*a[2][2]) - (a[2][1]*a[1][2])) -a[0][1] 
	 * (a[1][0] * a[2][2] - a[2][0] * a[1][2]) + a[0][2] * (a[1][0] 
	 * a[2][1] - a[2][0] * a[1][1]);
}

void vector3_subtract(float a[3], float b[3], float diff[3])
{
	diff[0] = a[0] - b[0];
	diff[1] = a[1] - b[1];
	diff[2] = a[2] - b[2];
	return;
}

void print_mat3b3(float mat[3][3])
{
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
			printf("%f   ", mat[i][j]);
		printf("\n");
	}
}

void compute_molecule_COM(crystal xtal, float com[3], int i)
{	/* computes the COM of a molecule (whose first atom is i) in
	a crystal. 
	*/
	com[0] = 0;
	com[1] = 0;
	com[2] = 0;
	
	int N = xtal.num_atoms_in_molecule;
	for(int j = i; j < i + N ; j++)
	{
		com[0] += xtal.Xcord[j]; 
		com[1] += xtal.Ycord[j];
		com[2] += xtal.Zcord[j];
	}
	com[0] /= N;
	com[1] /= N;
	com[2] /= N;
	
	//print_vec3(com);
	
}

void inverse_mat3b3(float Tinv[3][3], float T[3][3])
{
	float det = det_mat3b3(T);

    if(det == 0)
       {printf("error: det = 0, matrix has no inverse \n"); print_mat3b3(T);}
 
	for(int i = 0;i < 3; i++)
	{
    	for(int j = 0; j < 3; j++)
		{
			Tinv[j][i]=((T[(i+1)%3][(j+1)%3] * T[(i+2)%3][(j+2)%3]) - 
				(T[(i+1)%3][(j+2)%3]*T[(i+2)%3][(j+1)%3]))/ det;
		}
    }
}

/*
Compute the length of the first molecule in xtal
*/
float get_molecule_length(crystal xtal)
{
	int N = xtal.num_atoms_in_molecule;
	float com[3];
	compute_molecule_COM( xtal, com, 0);
	float max = 0;
	float dist = 0;
	int first_time = 1;

	for(int i = 0; i < N; i++)
	{
		float atom_vec[3] = {xtal.Xcord[i] , xtal.Ycord[i] , xtal.Zcord[i]};
		float dist_vec[3] = {0 , 0 , 0};
		vector3_subtract(atom_vec, com, dist_vec);
		dist = vector3_norm(dist_vec);

		if (first_time)
		{
			first_time = 0;
			max = dist;
		}
		else if (max < dist)
		{
			max = dist;
		}

	}
	//printf("molecule distance  = %f\n", dist);
	return 2*dist;
}


float pdist_appx(float T[3][3],
			float T_inv[3][3],
			float x1, 
			float x2,
			float x3,
			float y1,
			float y2,
			float y3  );

void pdist_2_appx(float T[3][3],
			float T_inv[3][3],
			float x1, 
			float x2,
			float x3,
			float y1,
			float y2,
			float y3,
			float *p_dist,
			float *p_dist_2  );

//#############check structure###########3
float CONSTANT_TOL = 0.001;
float MOL_SIM_TOL = 0.01;

static inline int int_floor(float x)
{
  int i = (int)x; /* truncate */
  return i - ( i > x ); /* convert trunc to floor */
}



/* The most function of this file. returns the distance between
 * two points under periodic boundary. T and T_inv are lattice
 * vector matrix and its inverse. T is in row major form. x1,x2,x3 are
 * the first point and y's are second in cartesian
 * returns the shortest distance as a float.
 */
float pdist(float T[3][3],
			float T_inv[3][3],
			float x1, 
			float x2,
			float x3,
			float y1,
			float y2,
			float y3  )
{
	float p_dist = pdist_appx(T,
			T_inv,
			x1, 
			x2,
			x3,
			y1,
			y2,
			y3  );
			
	float cartesian_distance[3] = {x1 - y1, x2 - y2, x3 - y3};
    float fractional_distance[3];
	for (int i = 0; i < 3; i++)
	{	
		float sum = 0;
		for (int j = 0; j < 3; j++)
			sum = sum + T_inv[j][i] * cartesian_distance[j];	

		float frac_part = sum - ((long)sum); //if -ve add one
		if (frac_part < 0)
			fractional_distance[i] = frac_part + 1;	
		else
			fractional_distance[i] = frac_part;
	}

	//find cartesian distance vector in the bounding box
	//reduced cartesian distance vector
	float red_cart_distance[3];
	for(int i = 0; i < 3; i++)
		red_cart_distance[i] = T[0][i] * fractional_distance[0] + 
							   T[1][i] * fractional_distance[1] +
							   T[2][i] * fractional_distance[2];
	
	//use pdist as radius of sphere are search
	int search_safety = 2;
	int limx = abs( int_floor ( p_dist/ T[0][0] ) ) + 1 + search_safety;
	int limy = abs( int_floor ( p_dist/ T[1][1] ) ) + 1 + search_safety;
	int limz = abs( int_floor ( p_dist/ T[2][2] ) ) + 1 + search_safety;
	
	float test_dist[3] = {0,0,0};
	
	for (int i = -limx; i <= limx; i++)
	for (int j = -limy; j <= limy; j++)
	for (int k = -limz; k <= limz; k++)
	{
		test_dist[0] = i*T[0][0] + j*T[1][0] + k*T[2][0] - red_cart_distance[0];
		test_dist[1] = j*T[1][1] + k*T[2][1] - red_cart_distance[1];
		test_dist[2] = k*T[2][2] - red_cart_distance[2];
		
		if (vector3_norm(test_dist) < p_dist )
			p_dist = vector3_norm(test_dist);
	}

	return p_dist;
}

/* Almost the same as p_dist, but returns the shortest and the second
 * shortest distance. return value is in argument. used for self image 
 * check
 */
void pdist_2(float T[3][3],
			float T_inv[3][3],
			float x1, 
			float x2,
			float x3,
			float y1,
			float y2,
			float y3,
			float *p_dist,
			float *p_dist_2  )
{
	*p_dist = 0.0;
	*p_dist_2 = 0.0;
	
	pdist_2_appx(T,
				T_inv,
				x1,
				x2,
				x3,
				y1,
				y2,
				y3,
				p_dist,
				p_dist_2);
	
	float cartesian_distance[3] = {x1 - y1, x2 - y2, x3 - y3};
    float fractional_distance[3];
	for (int i = 0; i < 3; i++)
	{	
		float sum = 0;
		for (int j = 0; j < 3; j++)
			sum = sum + T_inv[j][i] * cartesian_distance[j];	

		float frac_part = sum - ((int)sum); //if -ve add one
		if (frac_part < 0)
			fractional_distance[i] = frac_part + 1;	
		else
			fractional_distance[i] = frac_part;
	}

	float red_cart_distance[3];
	for(int i = 0; i < 3; i++)
		red_cart_distance[i] = T[0][i] * fractional_distance[0] + 
							   T[1][i] * fractional_distance[1] +
							   T[2][i] * fractional_distance[2];
	
	//use pdist as radius of sphere are search
	int search_safety = 2;
	int limx = abs( int_floor ( *p_dist/ T[0][0] ) ) + 1 + search_safety;
	int limy = abs( int_floor ( *p_dist/ T[1][1] ) ) + 1 + search_safety;
	int limz = abs( int_floor ( *p_dist/ T[2][2] ) ) + 1 + search_safety;
	
	float test_dist[3] = {0,0,0};
	
	for (int i = -limx; i <= limx; i++)
	for (int j = -limy; j <= limy; j++)
	for (int k = -limz; k <= limz; k++)
	{
		test_dist[0] = i*T[0][0] + j*T[1][0] + k*T[2][0] - red_cart_distance[0];
		test_dist[1] = i*T[0][1] + j*T[1][1] + k*T[2][1] - red_cart_distance[1];
		test_dist[2] = i*T[0][2] + j*T[1][2] + k*T[2][2] - red_cart_distance[2];
		
		float norm = vector3_norm(test_dist);
		if ( norm + CONSTANT_TOL < *p_dist )
		{/*
            printf("ijl = %d %d %d , norm = %f, red_cart_dist=%f\n", i, j, k, norm, vector3_norm(red_cart_distance));
            printf("red_cart_dist = %f %f %f\n", red_cart_distance[0], red_cart_distance[1], red_cart_distance[2]);
            printf("x =  %f, %f, %f \n", x1, x2, x3);
            printf("y =  %f, %f, %f \n", y1, y2, y3);
            printf("fractional dist =  %f %f %f \n", fractional_distance[0], fractional_distance[1], fractional_distance[2] );
           */
           // print_mat3b3(T);
           // print_mat3b3(T_inv);
			*p_dist_2 = *p_dist;
			*p_dist = norm;
		}
		else if ( norm + CONSTANT_TOL < *p_dist_2 && *p_dist > norm + CONSTANT_TOL)
			*p_dist_2 = norm;
	}
	
}



/*aproximate version of pdist()
 */ 
float pdist_appx(float T[3][3],
			float T_inv[3][3],
			float x1, 
			float x2,
			float x3,
			float y1,
			float y2,
			float y3  )
{
    //intialising variables
	float p_dist = 0;

	static float Q[8][3]={ {0,0,0},{1,0,0},{0,1,0},{0,0,1},
						   {1,1,0},{0,1,1},{1,0,1},{1,1,1}  };
		
	float cartesian_distance[3] = {x1 - y1, x2 - y2, x3 - y3};
	//vector3_subtract(x,y,cartesian_distance);

    /*computes tthe distance in fractinal space and reduces it 
     * inside first unit cube
     * #TODO: make this outside pdist().
     *  Needs to be done only once per crystal
     */
    float fractional_distance[3];
	for (int i = 0; i < 3; i++)
	{	
		float sum = 0;
		for (int j = 0; j < 3; j++)
			sum = sum + T_inv[j][i] * cartesian_distance[j];	

		float frac_part = sum - ((long)sum); //if -ve add one
		if (frac_part < 0)
			fractional_distance[i] = frac_part + 1;	
		else
			fractional_distance[i] = frac_part;
	}

	for (int z = 0; z < 8; z++)
	{
		float A[3] = {Q[z][0], Q[z][1], Q[z][2]};
		float dist_corner[3];

        dist_corner[0] = fractional_distance[0] - A[0];
		dist_corner[1] = fractional_distance[1] - A[1];
        dist_corner[2] = fractional_distance[2] - A[2];

		//distance vector to 8 corners in cartesian.
		float dist_z = 0;
		for (int i = 0; i < 3; i++)
		{
			float sum = 0;
			for (int j = 0; j < 3 ;j++)
				sum = sum + dist_corner[j] * T[j][i];	
			dist_z += sum*sum;
		}
        // length odistnce vector to the zth corner 
        // for finding the minimum distance (min dist_z)
		if (z == 0)
			p_dist = dist_z;
		else if (p_dist > dist_z)
			p_dist = dist_z;
	}
	return sqrt(p_dist); 
}

/* Almost the same as p_dist, but returns the shortest and the second
 * shortest distance. return value is in argument. used for self image 
 * check
 */
void pdist_2_appx(float T[3][3],
			float T_inv[3][3],
			float x1, 
			float x2,
			float x3,
			float y1,
			float y2,
			float y3,
			float *p_dist,
			float *p_dist_2  )
{
    //intialising variables
	*p_dist = 0;
	*p_dist_2 = 0;
	static float Q[8][3]={ {0,0,0},{1,0,0},{0,1,0},{0,0,1},
						   {1,1,0},{0,1,1},{1,0,1},{1,1,1}  };
		
	float cartesian_distance[3] = {x1 - y1, x2 - y2, x3 - y3};
	//vector3_subtract(x,y,cartesian_distance);

    /*computes tthe distance in fractinal space and reduces it 
     * inside first unit cube
     * #TODO: make this outside pdist().
     *  Needs to be done only once per crystal
     */
    float fractional_distance[3];
	for (int i = 0; i < 3; i++)
	{	
		float sum = 0;
		for (int j = 0; j < 3; j++)
			sum = sum + T_inv[j][i] * cartesian_distance[j];	

		float frac_part = sum - ((long)sum); //if -ve add one
		if (frac_part < 0)
			fractional_distance[i] = frac_part + 1;	
		else
			fractional_distance[i] = frac_part;
	}

    //computes the distances to all the corners of the unit cube
	for (int z = 0; z < 8; z++)
	{
		float A[3] = {Q[z][0], Q[z][1], Q[z][2]};
		float dist_corner[3];

        dist_corner[0] = fractional_distance[0] - A[0];
		dist_corner[1] = fractional_distance[1] - A[1];
        dist_corner[2] = fractional_distance[2] - A[2];

		//distance vector to 8 corners in cartesian.
		float dist_z = 0;
		for (int i = 0; i < 3; i++)
		{
			float sum = 0;
			for (int j = 0; j < 3 ;j++)
				sum = sum + dist_corner[j] * T[j][i];	
			dist_z += sum*sum;
		}
        // length odistnce vector to the zth corner 
        // for finding the minimum distance (min dist_z)
		if (z == 0)
		{
			*p_dist = dist_z;
		}
		else if (z == 1)
		{
			if (*p_dist > dist_z)
				{*p_dist_2 = *p_dist; *p_dist = dist_z;}
			else
				*p_dist_2 = dist_z;
		}
		else if (*p_dist > dist_z)
		{
			*p_dist_2 = *p_dist;
			*p_dist = dist_z;
		}
		else if (*p_dist_2 > dist_z)
		{
			*p_dist_2 = dist_z;
		}
		
	}
	*p_dist 	= sqrt(*p_dist);
	*p_dist_2	= sqrt(*p_dist_2);

}



/* convert the atoms char array into vdw radii information for structure checking
 * uses Bondii radii 
 */

void convert_atom2atom_vdw(char *atom,float *atom_vdw, int num_atoms)
{

	for (int i = 0; i < num_atoms; i++)
	{
		if      (atom[2*i] == 'C' && atom[2*i+1] == ' ')
		atom_vdw[i]=1.7;
		else if (atom[2*i] == 'H' && atom[2*i+1] == ' ')
		atom_vdw[i]=1.1;
		else if (atom[2*i] == 'N' && atom[2*i+1] == ' ')
		atom_vdw[i] = 1.55 ;
		else if (atom[2*i] == 'O' && atom[2*i+1] == ' ')
		atom_vdw[i] = 1.52;
		else if (atom[2*i] == 'F' && atom[2*i+1] == ' ')
		atom_vdw[i] = 1.47;
		else if (atom[2*i] == 'P' && atom[2*i+1] == ' ')
		atom_vdw[i] = 1.8;
		else if (atom[2*i] == 'S' && atom[2*i+1] == ' ')
		atom_vdw[i] = 1.8;
		else if (atom[2*i] == 'C' && atom[2*i+1] == 'l')
		atom_vdw[i] = 1.75;
		else if (atom[2*i] == 'B' && atom[2*i+1] == 'r')
		atom_vdw[i] = 1.85;
		else if (atom[2*i] == 'I' && atom[2*i+1] == ' ')
		atom_vdw[i] = 1.98;
		else if (atom[2*i] == 'B' && atom[2*i+1] == ' ')
		atom_vdw[i] = 1.92;
		else if (atom[2*i] == 'H' && atom[2*i+1] == 'e')
		atom_vdw[i] = 1.40;
		else if (atom[2*i] == 'N' && atom[2*i+1] == 'e')
		atom_vdw[i] = 1.54;
		else if (atom[2*i] == 'K' && atom[2*i+1] == 'r')
		atom_vdw[i] = 2.02;
		else if (atom[2*i] == 'S' && atom[2*i+1] == 'i')
		atom_vdw[i] = 2.10;		
		else
			{printf("***ERROR: atom2atom_vdw: atom not found -> %c%c\n", atom[2*i], atom[2*i+1]);exit(0);}
	//	printf("%d --> %f \n", i, atom_vdw[i]);
	}
}

/* used for checking self-image overlap. this has the information if the 
 * atoms are bonded.
 */
void calc_bond_length(float *bond_length,float *X, float *Y, float *Z, int N)
{
	
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			*(bond_length+i*N+j) = sqrt( (X[i] - X[j])*(X[i] - X[j]) + 
			(Y[i] - Y[j])*(Y[i] - Y[j]) + (Z[i] - Z[j])*(Z[i] - Z[j]) );
		//printf("%f \t" , *(bond_length+i*N+j));	
		}
		//printf("\n");
	}
}

int fast_screener_vdw(crystal xtal, float *vdw_matrix)
{
									//number of atoms in a molecule
	int N = xtal.num_atoms_in_molecule;	
	int m = xtal.Z;					//number of molecules in a unit cell;
								// numberof atom in a unit cell = N*m
	int total_atoms = N*m;	
	static float mol_len = 0;
    static int first_time = 1;

    if (first_time)
    {
        first_time = 0;
        mol_len = get_molecule_length(xtal);
    }

	float small_number = 2;
	
	for(int i = 0; i < total_atoms; i += N)
	{
		float com1[3];
		compute_molecule_COM( xtal, com1, i);
		for(int j = i + N; j < total_atoms; j += N)
		{
			//check if the molecule COM are far. if they are, dont 
			//bother checking distances
			if (j == i + N)
			{	float com2[3];
				compute_molecule_COM( xtal, com2, j);
				if( sqrt ((com1[0] - com2[0])*(com1[0] - com2[0])+
						  (com1[1] - com2[1])*(com1[1] - com2[1])+
						  (com1[2] - com2[2])*(com1[2] - com2[2]))
					<     (2*mol_len + small_number)                  )
					continue;
			}
			
			//molecule COM are close than molecule length
			else
			{
				for(int k = i; k < i + N; k++)
				{	
					for (int z = j; z < j + N; z++ )
					{
						if(		(xtal.Xcord[k] - xtal.Xcord[z])*
							   (xtal.Xcord[k] - xtal.Xcord[z])+
							   (xtal.Ycord[k] - xtal.Ycord[z])*
							   (xtal.Ycord[k] - xtal.Ycord[z])+
							   (xtal.Zcord[k] - xtal.Zcord[z])*
							   (xtal.Zcord[k] - xtal.Zcord[z])
							   < *(vdw_matrix + total_atoms*k + z) *
								 *(vdw_matrix + total_atoms*k + z)	)
						    return 0;
					}
				}
			}
		}
	}

	return 1;
}


/*checks the distance between two pair of molecules
 */
int check_pair_vdw_tier2( float T[3][3],
				float T_inv[3][3],
				float *X,
				float *Y,
				float *Z,
				float *vdw_matrix,
				int i,
				int j,
				int N,
				int total_atoms)
{
	
	for (int k = i; k < i + N; k++)
	{
		for (int l = j; l < j + N; l++)
		{
		//	printf("k=%d \t l= %d  \n", k, l);
			if (pdist_appx(T, T_inv, X[k], Y[k], Z[k], X[l], Y[l], Z[l]) <
				 *(vdw_matrix + k*total_atoms + l) 					)
			{
				//printf("atom1 = %f %f %f \n", X[k], Y[k], Z[k]);
				//printf("atom2 = %f %f %f \n", X[l], Y[l], Z[l]);
				//printf("%f  %f  %f", sr , atom_vdw[k], atom_vdw[l]);
				//printf("failed at %d atom and %d atom with distance %f,
				// expected %f \n", k+1%N, l+1%N,
				// pdist(T,T_inv, X[k], Y[k], Z[k], X[l], Y[l], Z[l]),
				//sr*(atom_vdw[k]+atom_vdw[l]) );
				return 0;
			}
		}
	}
//	printf("yes");
	return 1;
}

/*check with periodic image of one molecule: approximate but faster
*/
int check_self_vdw_tier2(float T[3][3], float T_inv[3][3],float *X,float *Y,
	float *Z, float *vdw_matrix, int i,int  N,int total_atoms, float *bond_length)
{

	int j,k, modj,modk;
	float pdist_kj, pdist_kj_2;
	for (j = i; j < i+N; j++)
	{
		modj = j % N;
		for (k = j ; k < i + N; k++)
		{
			modk = k % N;
			pdist_2_appx(T,T_inv, X[k], Y[k], Z[k], X[j], Y[j], Z[j], &pdist_kj, &pdist_kj_2);
			
			if (k == j && pdist_kj_2 < *(vdw_matrix + k*total_atoms + j) )
				return 0;
				
			if (pdist_kj - *(bond_length+modj*N+modk) < -0.02)
			{
				if (pdist_kj < *(vdw_matrix + k*total_atoms + j) )
				{
					//printf("failed at %d atom and %d atom with distance %f, expected %f , bondlength of %f\n",
					// k+1, j+1, pdist_kj,sr*(atom_vdw[k]+atom_vdw[j]), *(bond_length+modj*N+modk) );
					return 0 ;
				}
			}
			
			if(pdist_kj_2 < *(vdw_matrix + k*total_atoms + j) )
				return 0;
			
		}
	}
	return 1;
}


int check_pair_vdw_tier3( float T[3][3],
				float T_inv[3][3],
				float *X,
				float *Y,
				float *Z,
				float *vdw_matrix,
				int i,
				int j,
				int N,
				int total_atoms)
{
	
	for (int k = i; k < i + N; k++)
	{
		for (int l = j; l < j + N; l++)
		{
		//	printf("k=%d \t l= %d  \n", k, l);
			if (pdist(T, T_inv, X[k], Y[k], Z[k], X[l], Y[l], Z[l]) <
				 *(vdw_matrix + k*total_atoms + l) 					)
			{
				//printf("atom1 = %f %f %f \n", X[k], Y[k], Z[k]);
				//printf("atom2 = %f %f %f \n", X[l], Y[l], Z[l]);
				//printf("%f  %f  %f", sr , atom_vdw[k], atom_vdw[l]);
				//printf("failed at %d atom and %d atom with distance %f,
				// expected %f \n", k+1%N, l+1%N,
				// pdist(T,T_inv, X[k], Y[k], Z[k], X[l], Y[l], Z[l]),
				//sr*(atom_vdw[k]+atom_vdw[l]) );
				return 0;
			}
		}
	}
//	printf("yes");
	return 1;
}

/*check with periodic image of one molecule: rigourous and slow
*/
int check_self_vdw_tier3(float T[3][3], float T_inv[3][3],float *X,float *Y,
	float *Z, float *vdw_matrix, int i,int  N,int total_atoms, float *bond_length)
{

	int j,k, modj,modk;
	float pdist_kj, pdist_kj_2;
	for (j = i; j < i+N; j++)
	{
		modj = j % N;
		for (k = j ; k < i + N; k++)
		{
			modk = k % N;
			pdist_2(T,T_inv, X[k], Y[k], Z[k], X[j], Y[j], Z[j], &pdist_kj, &pdist_kj_2);
			
			if (k == j && pdist_kj_2 < *(vdw_matrix + k*total_atoms + j) )
				return 0;
				
			if (pdist_kj - *(bond_length+modj*N+modk) < -0.02)
			{
				if (pdist_kj < *(vdw_matrix + k*total_atoms + j) )
				{
					//printf("failed at %d atom and %d atom with distance %f, expected %f , bondlength of %f\n", k+1, j+1, pdist_kj,sr*(atom_vdw[k]+atom_vdw[j]), *(bond_length+modj*N+modk) );
					return 0 ;
				}
			}
			
			if(pdist_kj_2 < *(vdw_matrix + k*total_atoms + j) )
				return 0;
			
		}
	}
	return 1;
}


int check_structure_with_vdw_matrix(crystal random_crystal,
	float *vdw_matrix,
	int dim1,
	int dim2)
{
	int N = random_crystal.num_atoms_in_molecule;	//number of atoms in a molecule
	int m = random_crystal.Z;						//number of molecules in a unit cell;
													// numberof atom in a unit cell = N*m	
	int total_atoms = m*N;
	
	if (dim1 != dim2)
		{printf("ERROR:matrix not square\n"); exit(0);}

	if (dim1 != total_atoms)
		{printf("matrix size doesnt match total atoms\n"); exit(0);}

	if( fast_screener_vdw(random_crystal, vdw_matrix) == 0)
		{return 0;}
	
	int final_verdict =1;
	int check_val=1;
	float T_inv[3][3];
	int i,j;
	inverse_mat3b3(T_inv, random_crystal.lattice_vectors);

	//tier2 checks
	//start checking each pair of molecule
	for (i = 0; i < total_atoms; i = i + N)
	{
		for(j= i + N; j < total_atoms; j = j + N)
		{
			check_val = check_pair_vdw_tier2(random_crystal.lattice_vectors,
									   T_inv,
									   random_crystal.Xcord,
									   random_crystal.Ycord,
									   random_crystal.Zcord,
									   vdw_matrix,
									   i, 
									   j,
									   N,
									   total_atoms);
			
			//if(i == 0)
			//printf("check val = %d \n", check_val );
			if (check_val == 0) // 0 if check failed , 1 if it is ok
			{
				i = N*m + 1; //to break out of both loops if the check failed
				j = N*m + 1;
				final_verdict = 0;
			}
		}
	}
	
	//check self-image
	float *bond_length = malloc (N*N*sizeof(float));
	 calc_bond_length(bond_length,
					  random_crystal.Xcord,
					  random_crystal.Ycord,
					  random_crystal.Zcord,
					  N);
		
    if (final_verdict == 1)
    {
		for (i = 0; i < total_atoms; i = i + N)
		{
			check_val = check_self_vdw_tier2(random_crystal.lattice_vectors,
									T_inv,
									random_crystal.Xcord,
								   random_crystal.Ycord,
								   random_crystal.Zcord,
								   vdw_matrix,
								   i,
								   N,
								   total_atoms,
								   bond_length);
				
			if (check_val == 0)
			{
				i = total_atoms+ 1;
				//break;
			}
		}
	}
	if (check_val == 0)
		final_verdict = 0;

														//start tier-3 check
	
	if (final_verdict == 1)
	{
		for (i = 0; i < total_atoms; i = i + N)
		{
			for(j= i + N; j < total_atoms; j = j + N)
			{
				check_val = check_pair_vdw_tier3(random_crystal.lattice_vectors,
									   T_inv,
									   random_crystal.Xcord,
									   random_crystal.Ycord,
									   random_crystal.Zcord,
									   vdw_matrix,
									   i, 
									   j,
									   N,
									   total_atoms);
				
				if (check_val == 0) // 0 if check failed , 1 if it is ok
				{
					i = N*m + 1; //to break out of both loops if the check failed
					j = N*m + 1;
					final_verdict = 0;
				}
			}
		}
	} 
	
    if (final_verdict == 1)
    {
		for (i = 0; i < total_atoms; i = i + N)
		{
			check_val = check_self_vdw_tier3(random_crystal.lattice_vectors,
									T_inv,
									random_crystal.Xcord,
								   random_crystal.Ycord,
								   random_crystal.Zcord,
								   vdw_matrix,
								   i,
								   N,
								   total_atoms,
								   bond_length);
				
			if (check_val == 0)
			{
				i = total_atoms+ 1;
				final_verdict = 0;
			}
		}
	}
														//end tier-3

	free(bond_length);

	return final_verdict;
} 


//#######################################
void print_crystal(crystal* xtal)
{
	
	
	 int N = xtal->num_atoms_in_molecule;
	 int m = xtal->Z;
	 
	 printf("#Z = %d \n", xtal->Z);
	 printf("#napm = %d \n", xtal->num_atoms_in_molecule);
	
	for(int i = 0; i < 3; i++)
	{
		printf("lattice_vector %12f %12f %12f \n",
			xtal->lattice_vectors[i][0], xtal->lattice_vectors[i][1],
			xtal->lattice_vectors[i][2]);
	}
	
	for(int i = 0; i < N*m; i++)
	{
		printf("atom %12f %12f %12f %4c \n", xtal->Xcord[i],
			xtal->Ycord[i],  xtal->Zcord[i],  xtal->atoms[2*i]);
	}
}


void read_xtal_from_file(char *fname, crystal* xtal)
{
	FILE *fileptr;
	size_t len = 0;
	char *line = NULL;
	char *sub_line = NULL;
	int read;
	fileptr = fopen(fname,"r");

	if(!fileptr)
	{
		printf("***ERROR: xtal file not found \n");
		exit(0);
	}

	int total_atoms = 0;
	while ((read = getline(&line, &len, fileptr)) != -1)
	{
        if (strstr(line, "#") != NULL)
            continue;

		sub_line=strtok(line," ");
        //printf("%s \n" , sub_line);
        if(strcmp(sub_line, "atom") == 0)
            total_atoms++;
        else
            continue;
    }
    fclose(fileptr);

    int N = total_atoms;
    printf("total atoms = %d\n", total_atoms);
    xtal->Xcord = (float *)malloc(N*sizeof(float));
    
    xtal->Ycord = (float *)malloc(N*sizeof(float));
    xtal->Zcord = (float *)malloc(N*sizeof(float));
    xtal->atoms = (char *)malloc(2*N*sizeof(char));

    fileptr = fopen(fname,"r");
    int i = 0; 
    total_atoms = 0;

    while ((read = getline(&line, &len, fileptr)) != -1)
	{
        if (strstr(line, "#") != NULL)
            continue;

		sub_line=strtok(line," ");
        
        if(strcmp(sub_line, "atom") == 0)
         {
        	sub_line=strtok(NULL," ");
			(*xtal).Xcord[total_atoms]=atof(sub_line);
			sub_line=strtok(NULL," ");
			(*xtal).Ycord[total_atoms]=atof(sub_line);
			sub_line=strtok(NULL," ");
			(*xtal).Zcord[total_atoms]=atof(sub_line);
			sub_line=strtok(NULL," ");
			(*xtal).atoms[2*total_atoms]=*sub_line;
        	if(*(sub_line+1) == '\n' || *(sub_line+1) == ' ' || *(sub_line+1) == '\0' )
            	(*xtal).atoms[2*total_atoms+1]=' ';            
        	else
            	(*xtal).atoms[2*total_atoms+1]=*(sub_line+1);

            total_atoms++;
        }
        
        if (strcmp(sub_line, "lattice_vector") == 0)
        {
        	sub_line=strtok(NULL," ");
        	xtal->lattice_vectors[i][0] = atof(sub_line);
        	sub_line=strtok(NULL," ");
        	xtal->lattice_vectors[i][1] = atof(sub_line);
        	sub_line=strtok(NULL," ");
        	xtal->lattice_vectors[i][2] = atof(sub_line);

        	i++;
        } 

        else
            continue;
   	}
    fclose(fileptr);

    xtal->Z = 0;
    xtal->spg = 0;

}

void read_geometry(molecule* mol)
{
    FILE *fileptr;
	size_t len = 0;
    int read;	
	char *line = NULL;
	char *sub_line = NULL;
	int i = 0;
    int atom_count = 0;

    //find_number_of atoms
	fileptr = fopen("geometry.in","r");
	//check if file exits
	if(!fileptr)
	{
		printf("***ERROR: no geometry.in file \n");
		exit(0);
	}

     while ((read = getline(&line, &len, fileptr)) != -1)
	{
        if (strstr(line, "#") != NULL)
            continue;

		sub_line=strtok(line," ");
        //printf("%s \n" , sub_line);
        if(strcmp(sub_line, "atom") == 0)
            atom_count++;
        else
            continue;
    }
    fclose(fileptr);

    //printf("Total number of atoms in molecule = %d\n", atom_count);
    int N = atom_count;
	//memory allocation
	(*mol).atoms = (char *)malloc(2*N*sizeof(char));
	(*mol).X = (float *)malloc(N*sizeof(float));
	(*mol).Y = (float *)malloc(N*sizeof(float));
	(*mol).Z = (float *)malloc(N*sizeof(float));
	
	//creates a file pointer and opens geometry.in
	    
    fileptr = fopen("geometry.in","r");
    while ((read = getline(&line, &len, fileptr)) != -1)
	{
		//printf("Retrieved line of length %zu :\n", read);
		//printf("%s", line);
        if (strstr(line, "#") != NULL)
            continue;

		sub_line=strtok(line," ");
        //printf("%s \n" , sub_line);
        if(strcmp(sub_line, "atom") == 0)
            atom_count++;
        else
            continue;
		sub_line=strtok(NULL," ");
		(*mol).X[i]=atof(sub_line);
		//	printf("%f \t",X[i-3]);	
		sub_line=strtok(NULL," ");
		(*mol).Y[i]=atof(sub_line);
		//	printf("%f \t",Y[i-3]);
		sub_line=strtok(NULL," ");
		(*mol).Z[i]=atof(sub_line);
		//	printf("%f \t",Z[i-3]);	
		sub_line=strtok(NULL," ");
		(*mol).atoms[2*i]=*sub_line;
        if(*(sub_line+1) == '\n' || *(sub_line+1) == ' ' || *(sub_line+1) == '\0' )
            (*mol).atoms[2*i+1]=' ';            
        else
            (*mol).atoms[2*i+1]=*(sub_line+1);
			//printf("char = %c \n",mol->atoms[2*i+1]); exit(0);
            //printf("string = %c l %c \n", *sub_line, *(sub_line+1));	
		i++;

   	}
 

	fclose(fileptr);
   // printf("atoms = %d \n", atom_count);	
	mol->num_of_atoms = N;
	

}

void create_vdw_matrix_from_sr( molecule *mol,
								float *vdw_matrix,
								float sr,
								int Z);
int structure_checker(float L[3][3],
	float *X ,
	int dim1,
	float *Y,
	int dim2,
	float *Z,
	int dim3,
	float *vdw_cutoff,
	int dim4, 
	int dim5,
	int *mol_id,
	int num_mols
	);

int read_geometry_out(char* file_name, crystal *xtal);

int main()
{
    crystal xtal;
    float sr = 0.85;
    molecule *mol = (molecule*)malloc(sizeof(molecule));
    read_geometry(mol);
    float vdw_matrix[3600];
    create_vdw_matrix_from_sr(mol, vdw_matrix, sr, 2); 
    read_xtal_from_file("test.in", &xtal);
    xtal.num_atoms_in_molecule = 30;
    xtal.Z = 2;
    //print_crystal(  &xtal);
    int mol_id[] = {0, 30};
    int num_mols = 2;
    /*
    int verdict = structure_checker(xtal.lattice_vectors, xtal.Xcord, 60, xtal.Ycord, 60, xtal.Zcord, 60, vdw_matrix, 60, 60, mol_id, num_mols);
    printf("verdict = %d\n", verdict );
    */

    //test reaf_geometry.out
    /*
    
    crystal txtal;
    read_geometry_out(fname, &txtal);
    printf("read file\n");
    print_crystal( &txtal);
    read_geometry_out(fname, &txtal);
    print_crystal( &txtal);
    */
    char fname[] = "good_geometry.out";
    //char fname[] = "junk.in";

    int ret_val = read_geometry_out(fname, &xtal);
    int count = 1;
    while ( ret_val )
    {
    	int verdict = structure_checker(xtal.lattice_vectors, xtal.Xcord, 60, xtal.Ycord, 60, xtal.Zcord, 60, vdw_matrix, 60, 60, mol_id, num_mols); 

        int verdict2 = check_structure_with_vdw_matrix(xtal, vdw_matrix, 60, 60);

        if (verdict != verdict2)
        {
        	printf("disagree - verdict = %d, verditct2 = %d count = %d \n",verdict, verdict2, count);

        }

		
        if (!verdict)
        	printf("failed new checks %d\n", count);

    	ret_val = read_geometry_out(fname, &xtal);
    	count++;
    }

}


#define CONSTANT_TOL 0.001
#define MAXVDW 2.1
#define NORM2(a,b) ( (a[0]-b[0])*(a[0]-b[0]) +\
 					 (a[1]-b[1])*(a[1]-b[1]) +\
  					 (a[2]-b[2])*(a[2]-b[2]) )

#define NORMSQR(a) ( a[0]*a[0] + a[1]*a[1] + a[2]*a[2] )
#define SQR(x) (x*x)

float find_mol_len(float *X, float *Y, float *Z, int len);
int check_pairwise_if_mol_close(float L[3][3], float *vdw_matrix, int dim1,\
 float *X, float *Y, float *Z, int id1, int N1, int id2, int N2,\
 float com1[3], float com2[3], float sum_len);
int check_pairwise(float L[3][3], float *vdw_matrix, int dim, float *X, float *Y, float *Z, \
	int id1, int id2, int N1, int N2, int i, int j, int k);



void create_vdw_matrix_from_sr( molecule *mol,
								float *vdw_matrix,
								float sr,
								int Z)
{	int num_atoms_in_molecule = mol->num_of_atoms;
	float *atom_vdw_vector = (float *) malloc( num_atoms_in_molecule*
							sizeof(float) );
	//create vector with vdw radii of size num of atoms in mol
	convert_atom2atom_vdw( mol->atoms, atom_vdw_vector,
							num_atoms_in_molecule);
							
	//now create matrix
	int dim_vdw_matrix = num_atoms_in_molecule * Z ;
	for(int i = 0; i < dim_vdw_matrix; i++)
		for(int j = 0; j < dim_vdw_matrix; j++)
			{
				*(vdw_matrix + i*dim_vdw_matrix +j)
					=  (*(atom_vdw_vector + i%num_atoms_in_molecule) +
						*(atom_vdw_vector + j%num_atoms_in_molecule) ) *
						sr;
			}
	free(atom_vdw_vector);
	
	
}

/* for pygenarris API only
* Assume lattice to be lower triangular in row major form
*/
int structure_checker(float L[3][3],
	float *X ,
	int dim1,
	float *Y,
	int dim2,
	float *Z,
	int dim3,
	float *vdw_cutoff,
	int dim4, 
	int dim5,
	int *mol_id,
	int num_mols
	)
{
	//preliminary checks
	if(dim1 != dim2 || dim2 != dim3 || dim3 != dim1)
	{
		printf("***ERROR: length of coordinate arrays don't match\n");
		return -1;
	}

	if(dim4 != dim5)
	{
		printf("***ERROR: vdw cutoff matrix should be a square matrix\n");
		return -1;
	}

	if(!check_lower_triangular(L))
	{
		printf("***ERROR: lattice_vectors should be a lower triangular matrix in row major form\n");
		return -1;
	}	

	//find number of atoms in each molecule
	int num_atoms_in_molecule[num_mols];
	for(int i = 0; i < num_mols -1; i++)
	{
		num_atoms_in_molecule[i] = mol_id[i+1] - mol_id[i];
		//printf("%s = %d\n", "num_atoms_molecule = ", num_atoms_in_molecule[i]);
	}
	//for last mol:
	num_atoms_in_molecule[num_mols-1] = dim1 - mol_id[num_mols-1];

	// estimate molecule lengths
	float mol_len[num_mols];
	for(int i = 0 ; i < num_mols; i++)
	{
		mol_len[i] = find_mol_len(X + mol_id[i],
								 Y + mol_id[i],
								 Z + mol_id[i],
								 num_atoms_in_molecule[i]);
		//printf("molLen = %f\n", mol_len[i]);
	}

	// find com
	float com[num_mols * 3];
	for(int i = 0; i < num_mols; i++)
	{
		//printf("num_atoms_in_molecule = %d\n", num_atoms_in_molecule[i]);
		find_mol_com(X + mol_id[i], Y +  mol_id[i], Z +  mol_id[i], num_atoms_in_molecule[i], com + 3*i);
		//printf("i = %d , com  === %f, %f, %f\n",i, com[3*i], com[3*i +1], com[3*i+2] );
	}

	for(int i = 0; i < num_mols; i++)
	{
		for(int j = i; j < num_mols; j++)
		{
			//check pair of selcted molecules.
			//printf("i, j = %d, %d \n", i, j );
			float sum_len = mol_len[i] + mol_len[j];
			int verdict = check_pairwise_if_mol_close(L, vdw_cutoff, dim1, X, Y, Z, mol_id[i], num_atoms_in_molecule[i], mol_id[j], num_atoms_in_molecule[j], com+3*i, com+3*j, sum_len);
			if (!verdict)
			{
				return 0 ;
			}
		}
	}

	return 1;
}


int check_pairwise_if_mol_close(float L[3][3], float *vdw_matrix, int dim, float *X, float *Y, float *Z, int id1, int N1, int id2, int N2, float com1[3], float com2[3], float sum_len)
{
	int xmax = fabs(sum_len/L[0][0]) + 2; // can't ignore nearest cells +2 ; loop goes until xmax-1
	int ymax = (sum_len + xmax*fabs(L[1][0]/L[1][1]) ) + 2;
	int zmax = (sum_len + xmax*fabs(L[2][0]) + ymax*fabs(L[2][0]/L[2][2]) ) + 2;
	//printf("xmax, ymax, zmax = %d %d %d\n", xmax, ymax, zmax);

	for(int i = -xmax; i < xmax; i++)
	for(int j = -ymax; j < ymax; j++)
	for(int k = -zmax; k < zmax; k++)
	{
		// Skip checking the molecule with itself
		if (id1 == id2 && i == 0 && j == 0 && k == 0)
			continue;

		float dist[3] = {com1[0] - com2[0] + i*L[0][0] + j*L[1][0] + k*L[2][0], 
						 com1[1] - com2[1]             + j*L[1][1] + k*L[2][1],
						 com1[2] - com2[2]                         + k*L[2][2]};

		
		if (NORMSQR(dist) < (sum_len/2 + MAXVDW)*(sum_len/2 + MAXVDW))
		{
			/*
			printf("i, j, k = %d, %d, %d dist = %f %f %f\n",i, j, k, dist[0], dist[1], dist[2] );
			printf("sumlen = %f\n", sum_len);
			printf("dist = %f\n", NORMSQR(dist));
			printf("close; ijk = %d %d %d\n", i,j,k );
			*/
			int ret_val = check_pairwise(L, vdw_matrix, dim, X, Y, Z, id1, id2, N1, N2, i, j, k);
			if(!ret_val)
			{

				return 0;
			}
		}
		
	}
	return 1;
}

int check_pairwise(float L[3][3], float *vdw_matrix, int dim, float *X, float *Y, float *Z, int id1, int id2, int N1, int N2, int i, int j, int k)
{
	for(int i1 = 0; i1 < N1; i1++)
	for(int i2 = 0; i2 < N2; i2++)
	{
		int atom1 = i1 + id1;
		int atom2 = i2 + id2;

		float dist[3] = {X[atom1] - X[atom2] + i*L[0][0] + j*L[1][0] + k*L[2][0],
						 Y[atom1] - Y[atom2] + j*L[1][1] + k*L[2][1],
						 Z[atom1] - Z[atom2] + k*L[2][2]
						};

		if (NORMSQR(dist) < SQR(*(vdw_matrix + atom1*dim + atom2)) )
		{
			/*printf("failedat i=%d, j=%d, k=%d, atom1 = %d, atom2=%d \n", i, j, k, atom1, atom2);
			printf("cutoff = %f, dist = %f\n", *(vdw_matrix + atom1*dim + atom2),  sqrt(NORMSQR(dist)) );
			printf("atom1 = %f %f %f, atom2 = %f %f %f\n", X[atom1], Y[atom1], Z[atom1], X[atom2], Y[atom2], Z[atom2]);
			*/
			return 0;
		}
	}

	return 1;
}

int check_lower_triangular(float L[3][3])
{
	if (L[0][1] < CONSTANT_TOL &&
	    L[0][2] < CONSTANT_TOL &&
	    L[1][2] < CONSTANT_TOL
	    )
	    return 1;

	return 0;
}

float find_mol_len(float *X, float *Y, float *Z, int len)
{
	float com[3] = {0, 0, 0};
	find_mol_com(X, Y, Z, len, com);
	//printf("com = %f %f %f\n", com[0], com[1], com[2]);
	float first_atom[3] = {X[0], Y[0], Z[0]};
	float max = NORM2(com, first_atom);
	float dist_sqr = 0;

	for(int i = 0; i < len; i++)
	{
		float atom_vec[3] = {X[i], Y[i], Z[i]};
		dist_sqr = NORM2(atom_vec, com);

		if (max < dist_sqr)
			max = dist_sqr;
	}

	return 2*sqrt(max);
}



void find_mol_com(float *X, float *Y, float *Z, int len, float com[3])
{
	com[0] = 0;
	com[1] = 0;
	com[2] = 0;

	for(int i = 0; i < len; i++)
	{
		com[0] += X[i];
		com[1] += Y[i];
		com[2] += Z[i];
	}

	com[0] /= len;
	com[1] /= len;
	com[2] /= len;
}


int read_geometry_out(char* file_name, crystal *xtal)
{
  static FILE *fptr = NULL;
  char * line=NULL;
  char * sub_line=NULL;
  int num_atom=0;
  size_t len=0;
  int read;
  int spg;
  int Z;
  int wyckoff_position = 0;
  int Zp=0;
  int i;
  int counter;
  int num_atoms_in_molecule;

  static int first_time = 1;
  if (first_time)
  {
  	fptr= fopen(file_name,"r");
  	first_time = 0;
  }

  if (!fptr)
  {
       printf("***ERROR: Cannot open the file %s\n", file_name);
       exit(0);
  }

	while (read = getline(&line,&len,fptr) != -1)
 	{
    	if (strstr(line,"#") != NULL || *line=='\n')
	    {	
			if (strstr(line,"attempted_spacegroup")!=NULL)
		    {
			     sub_line = strtok(line," ");
			     sub_line = strtok(NULL," ");
			     sub_line = strtok(NULL," ");
			     spg = atof(sub_line);
		    }

			else if (strstr(line,"Z")!=NULL)
			{
			    sub_line = strtok(line," ");
			    sub_line = strtok(NULL," ");
			    sub_line = strtok(NULL," ");
			    Z = atof(sub_line);
		        i = 0;
		        counter = 0;
		    }

		    else if (strstr(line,"number_of_atoms_in_molecule") != NULL)
		    {
		    	sub_line = strtok(line," ");
		    	sub_line = strtok(NULL," ");
			    sub_line = strtok(NULL," ");
			    num_atoms_in_molecule = atof(sub_line);
		    }

			else if (strstr(line,"END") != NULL)
		    {
			     xtal->spg = spg;
			     xtal->wyckoff_position = wyckoff_position;
			     xtal->num_atoms_in_molecule = num_atoms_in_molecule;
			     xtal->Z = Z;
			     xtal->Zp = Zp;
			     return 1;
		    }

 		    continue;
	    }

	    sub_line = strtok(line,"    ");
	    if (strcmp(sub_line,"lattice_vector") == 0)
		{
		    sub_line = strtok(NULL,"    ");	
		    xtal->lattice_vectors[i][0] = atof(sub_line);
		    sub_line = strtok(NULL,"    ");
		    xtal->lattice_vectors[i][1] = atof(sub_line);
		    sub_line = strtok(NULL,"    ");
		    xtal->lattice_vectors[i][2] = atof(sub_line);
		    if (!i)
		    { 
		     	num_atom = num_atoms_in_molecule *Z;  
		        xtal->Xcord = (float *) malloc (num_atom*sizeof(float));
	     		xtal->Ycord = (float *) malloc (num_atom*sizeof(float));
	     		xtal->Zcord = (float *) malloc (num_atom*sizeof(float));
	     		xtal->atoms = (char *) malloc (2*num_atom*sizeof(char));
		     }
		    i++;
		}

	    if (strcmp(sub_line,"atom")==0)
		{
		    xtal->Xcord[counter] = atof(strtok(NULL,"    "));
		    xtal->Ycord[counter] = atof(strtok(NULL,"    "));
		    xtal->Zcord[counter] = atof(strtok(NULL,"    "));

		    sub_line = strtok(NULL," ");
		    xtal->atoms[2*counter] = *sub_line;
            if(*(sub_line+1) == '\n' || *(sub_line+1) == ' ' || *(sub_line+1) == '\0' )
            	xtal->atoms[2*counter+1] = ' ';            
        	else
            	xtal->atoms[2*counter+1] = *(sub_line+1);
		    counter++;
		}

    }

    fclose(fptr);
    return 0;
    
}