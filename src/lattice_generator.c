#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "lattice_generator.h"
#include "randomgen.h"
#include "niggli.h"

#define LOWB 2.5  //lower bound for length of lattice vector
#define PI 3.141592653

//****
//
//***

//angles are in radians

void gen_triclinic_lattice(float lattice_vector[3][3], 
	float target_volume, float max_angle, float min_angle)
{
	int attempt = 0;
    float obliqness = 0.6;
    float random01;
    float x, y, z;
    do {x = normal_dist_ab(1, obliqness); } while(x < 0);
    do {y = normal_dist_ab(1, obliqness); } while(y < 0);
    do {z = normal_dist_ab(1, obliqness); } while(z < 0);
    float factor =  cbrt (target_volume / (x*y*z) );
    float ax = x*factor;
    float by = y*factor;
    float cz = z*factor;
	//select principal components
/*
	float random01 = uniform_dist_01();
	float ax = random01*(target_volume/(LOWB*LOWB) - LOWB) + LOWB;
	random01 = uniform_dist_01();
	float by = random01*(target_volume/(ax*LOWB) - LOWB) +LOWB;
	float cz = target_volume/(ax*by); 
*/


	//generate random angles
	random01 = uniform_dist_01();;
	float beta = random01*(max_angle - min_angle) + min_angle; 
	random01 = uniform_dist_01();
	float gamma = random01*(max_angle - min_angle) + min_angle;
	float alpha;
	float cosbeta2;
	float cosbeta3;
	do
	{
		random01 = uniform_dist_01();;
		alpha = random01*(max_angle - min_angle) + min_angle;
		cosbeta2 = (cos(alpha) - cos(gamma)*cos(beta))/sin(gamma);
		cosbeta3 = sqrt(1 - cos(beta)*cos(beta) - cosbeta2*cosbeta2);
		attempt++;
		if (attempt > 1000)
		{	
			printf("*WARNING: bad input*\n");	
		}
	} 
	while(abs(cosbeta2) >= 1 ||
		(cos(beta)*cos(beta) + cosbeta2*cosbeta2) > 1 );

	float bx = by/tan(gamma);
	float modc = cz/cosbeta3;
	float cx = modc*cos(beta);
	float cy = modc*cosbeta2;

	lattice_vector[0][0] = ax;
	lattice_vector[1][1] = by;
	lattice_vector[2][2] = cz;

	lattice_vector[0][1] = 0;
	lattice_vector[0][2] = 0;
	lattice_vector[1][2] = 0;

	lattice_vector[1][0] = bx;
	lattice_vector[2][0] = cx;
	lattice_vector[2][1] = cy;
}


void gen_monoclinic_lattice(float lattice_vector[3][3], 
	float target_volume, float max_angle, float min_angle)
{
	float random01 =uniform_dist_01();
	float ax = random01*(target_volume/(LOWB*LOWB) - LOWB) + LOWB;
	random01 = uniform_dist_01();
	float by = random01*(target_volume/(ax*LOWB) - LOWB) +LOWB;
	float cz = target_volume/(ax*by); 
	
	random01 = uniform_dist_01();;
	float beta = random01*(max_angle - min_angle) + min_angle;

	float cx = cz / tan(beta); 

	lattice_vector[0][0] = ax;
	lattice_vector[1][1] = by;
	lattice_vector[2][2] = cz;

	lattice_vector[0][1] = 0;
	lattice_vector[0][2] = 0;
	lattice_vector[1][2] = 0;

	lattice_vector[1][0] = 0;
	lattice_vector[2][0] = cx;
	lattice_vector[2][1] = 0;
}

void gen_orthorhombic_lattice(float lattice_vector[3][3], 
	float target_volume)
{
	float random01 = uniform_dist_01();
	float ax = random01*(target_volume/(LOWB*LOWB) - LOWB) + LOWB;
	random01 = uniform_dist_01();
	float by = random01*(target_volume/(ax*LOWB) - LOWB) +LOWB;
	float cz = target_volume/(ax*by); 

	lattice_vector[0][0] = ax;
	lattice_vector[1][1] = by;
	lattice_vector[2][2] = cz;

	lattice_vector[0][1] = 0;
	lattice_vector[0][2] = 0;
	lattice_vector[1][2] = 0;

	lattice_vector[1][0] = 0;
	lattice_vector[2][0] = 0;
	lattice_vector[2][1] = 0;
}


void gen_tetragonal_lattice(float lattice_vector[3][3],
	float target_volume)
{
	float random01 = uniform_dist_01();
	float cz = random01*(target_volume/(LOWB*LOWB) - LOWB) + LOWB;
	float by = sqrt (target_volume/cz) ; 
	float ax = by; 

	lattice_vector[0][0] = ax;
	lattice_vector[1][1] = by;
	lattice_vector[2][2] = cz;

	lattice_vector[0][1] = 0;
	lattice_vector[0][2] = 0;
	lattice_vector[1][2] = 0;

	lattice_vector[1][0] = 0;
	lattice_vector[2][0] = 0;
	lattice_vector[2][1] = 0;
}

void gen_hexagonal_lattice(float lattice_vector[3][3], 
	float target_volume)
{
	float random01 = uniform_dist_01();
	float cz = random01*(target_volume/(LOWB*LOWB) - LOWB) + LOWB;
	float ax = sqrt (target_volume/cz) ; 
	float gamma = 120* PI/180;
	float by = ax*sin(gamma);
	float bx = ax*cos(gamma);

	lattice_vector[0][0] = ax;
	lattice_vector[1][1] = by;
	lattice_vector[2][2] = cz;

	lattice_vector[0][1] = 0;
	lattice_vector[0][2] = 0;
	lattice_vector[1][2] = 0;

	lattice_vector[1][0] = bx;
	lattice_vector[2][0] = 0;
	lattice_vector[2][1] = 0;
}

void gen_cubic_lattice(float lattice_vector[3][3],
	float target_volume)
{
	float a = cbrt(target_volume); 

	lattice_vector[0][0] = a;
	lattice_vector[1][1] = a;
	lattice_vector[2][2] = a;

	lattice_vector[0][1] = 0;
	lattice_vector[0][2] = 0;
	lattice_vector[1][2] = 0;

	lattice_vector[1][0] = 0;
	lattice_vector[2][0] = 0;
	lattice_vector[2][1] = 0;
}


void generate_lattice(float lattice_vector[3][3], int spg,
 float max_angle, float min_angle, float target_volume)
{
	
	if(spg < 1 || spg > 230)
	printf("***ERROR: generate_lattice: spg out of bounds***");
	
	else if (spg <= 2)
	gen_triclinic_lattice(lattice_vector, target_volume, max_angle, min_angle);

	else if (spg <= 15)
	gen_monoclinic_lattice(lattice_vector, target_volume, max_angle, min_angle);

	else if (spg <= 74)
	gen_orthorhombic_lattice(lattice_vector, target_volume);

	else if (spg <= 142)
	gen_tetragonal_lattice(lattice_vector, target_volume);

	else if (spg <= 167)
	gen_hexagonal_lattice(lattice_vector, target_volume);
	//same as heaxagonal?
	
	else if (spg <= 194) 
	gen_hexagonal_lattice(lattice_vector, target_volume);

	else if (spg <= 230)
	gen_cubic_lattice(lattice_vector, target_volume);
	
	//standardise_lattice(lattice_vector, spg);

	return;
}

//create a large volume lattice for testing compatiility
void generate_fake_lattice(float lattice_vector[3][3], int spg)
{
	const float ax = 15;
	
	
	if(spg < 1 || spg > 230)
	printf("***ERROR: generate_lattice: spg out of bounds***");
	
	else if (spg <= 2)
	{	
		lattice_vector[0][0] = ax;
		lattice_vector[1][1] = 0.8*ax;
		lattice_vector[2][2] = 0.5*ax;

		lattice_vector[0][1] = 0;
		lattice_vector[0][2] = 0;
		lattice_vector[1][2] = 0;

		lattice_vector[1][0] = 0.8*ax*tan(10*PI/180);
		lattice_vector[2][0] = 0.5*ax*cos(85*PI/180);
		lattice_vector[2][1] = 0.5*ax*cos(70*PI/180);
	}
	
	else if ( spg <= 15 )
	{
		lattice_vector[0][0] = ax;
		lattice_vector[1][1] = 0.8*ax;
		lattice_vector[2][2] = 0.5*ax;

		lattice_vector[0][1] = 0;
		lattice_vector[0][2] = 0;
		lattice_vector[1][2] = 0;

		lattice_vector[1][0] = 0;
		lattice_vector[2][0] = 0.5*ax/tan(70*PI/180);
		lattice_vector[2][1] = 0;
	}
	
	else if (spg <= 74)
	{
		lattice_vector[0][0] = ax;
		lattice_vector[1][1] = 0.8*ax;
		lattice_vector[2][2] = 0.5*ax;

		lattice_vector[0][1] = 0;
		lattice_vector[0][2] = 0;
		lattice_vector[1][2] = 0;

		lattice_vector[1][0] = 0;
		lattice_vector[2][0] = 0;
		lattice_vector[2][1] = 0;
	}
	
	else if (spg <= 142)
	{
		lattice_vector[0][0] = ax;
		lattice_vector[1][1] = ax;
		lattice_vector[2][2] = 0.5*ax;

		lattice_vector[0][1] = 0;
		lattice_vector[0][2] = 0;
		lattice_vector[1][2] = 0;

		lattice_vector[1][0] = 0;
		lattice_vector[2][0] = 0;
		lattice_vector[2][1] = 0;
	}
	
	else if (spg <= 194)
	{
		lattice_vector[0][0] = ax;
		lattice_vector[1][1] = ax*sin(120*PI/180);
		lattice_vector[2][2] = 0.5*ax;

		lattice_vector[0][1] = 0;
		lattice_vector[0][2] = 0;
		lattice_vector[1][2] = 0;

		lattice_vector[1][0] = ax*cos(120*PI/180);
		lattice_vector[2][0] = 0;
		lattice_vector[2][1] = 0;
	}
	
	else if (spg <= 230)
	{
		lattice_vector[0][0] = ax;
		lattice_vector[1][1] = ax;
		lattice_vector[2][2] = ax;

		lattice_vector[0][1] = 0;
		lattice_vector[0][2] = 0;
		lattice_vector[1][2] = 0;

		lattice_vector[1][0] = 0;
		lattice_vector[2][0] = 0;
		lattice_vector[2][1] = 0;
	}
}

static inline float fmodulo (float n, float d)
{
	long q =n/d;
	float r = n - q*d ;
	return r;	
}

void standardise_lattice( float lattice[3][3], int spg)
{
/*
    double nlattice[9] = {lattice[0][0], lattice[0][1], lattice[0][2],
                         lattice[1][0], lattice[1][1], lattice[1][2],
                         lattice[2][0], lattice[2][1], lattice[2][2]
                        };
    
    int status = niggli_reduce(nlattice, 0.0001);

    if(status)
    {
        lattice[0][0] = nlattice[0];
        lattice[0][1] = nlattice[1];
        lattice[0][2] = nlattice[2];
        
        lattice[1][0] = nlattice[3];
        lattice[1][1] = nlattice[4];
        lattice[1][2] = nlattice[5];

        lattice[2][0] = nlattice[6];
        lattice[2][1] = nlattice[7];
        lattice[2][2] = nlattice[8];
    }
*/

    
	//triclinic
	if (spg == 1 || spg == 2)
	{
		float by = lattice[1][1];
		float cy = lattice[2][1];
        float cx = lattice[2][0];
        float bx = lattice[1][0];
        long q = cy/by;
		float cy_new = cy - q * by;
        float cx_new = cx - q * bx;
		lattice[2][1] = cy_new;
        lattice[2][0] = cx_new;
		
		bx = lattice[1][0];
		float ax = lattice[0][0];
		float bx_new = fmodulo (bx, ax);
		lattice[1][0] = bx_new;
		
		cx = lattice[2][0];
		cx_new = fmodulo(cx, ax);
		lattice[2][0] = cx_new;
	}
	
	else if (spg > 2 && spg < 16)
	{
		float bx = lattice[1][0];
		float ax = lattice[0][0];
		float bx_new = fmodulo (bx, ax);
		lattice[1][0] = bx_new;
		
		float cx = lattice[2][0];
		float cx_new = fmodulo(cx, ax);
		lattice[2][0] = cx_new;
		
	}
	
	else if (spg <= 0 || spg > 230)
		printf("***ERROR: lattice_generator: standardise_lattice invalid spg");
	
}












