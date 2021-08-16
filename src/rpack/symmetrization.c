
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "../spglib.h"
#include "../algebra.h"
#include "symmetrization.h"

static void quaternion2matrix(float *quat, float *mat, const int inv);
static void matrix2quaternion(float *rot, float *quat, int *inv);
static void symmetrize_matrix(float *mat, int dim, int spg);
static void symmetrize_vector(float *vec, int dim, float lattice[3][3], int spg);

static void test_symmetrize_state();

/*
  Symmetrizes a vector wrt symmetry operations of a spg.
  Can be used for symmetrizing positions and position gradients.
  vec, size -> input vector and its dimensions. Shape of vec is (dim, 3).
               dim should be equal to number of symm operation or Z for now
	       This needs to be generalized for cocrystals
  lattice   -> lattice vectors in row major form
  spg       -> spacegroup which defines the symm operations.
*/
static void symmetrize_vector(float *vec, int dim, float lattice[3][3], int spg)
{
    // Get symm operations
    double translations[192][3];
    int rotations[192][3][3];
    int hall_number = hall_number_from_spg(spg);
    int Z = spg_get_symmetry_from_database(rotations,
					   translations,
					   hall_number);
    assert(dim == Z);
    
    // Get inverse lattice
    float inverse_lattice[3][3], inverse_lattice_t[3][3], lattice_t[3][3];
    inverse_mat3b3(inverse_lattice, lattice);

    mat3b3_transpose(inverse_lattice_t, inverse_lattice);
    mat3b3_transpose(lattice_t, lattice);

    // Convert to fractional space
    float vec_frac[3*dim];
    for(int i = 0; i < dim; i++)
	vector3_mat3b3_multiply(inverse_lattice_t, vec + 3*i, vec_frac +3*i);

    // Array to compute symm average
    float vec_symm[3] = {0};
    
    for(int op = 0; op < Z; op++)
    {
	int (*rot_i)[3] = rotations[op];
	float rot[3][3] = {{rot_i[0][0], rot_i[0][1],rot_i[0][2]},
			   {rot_i[1][0], rot_i[1][1],rot_i[1][2]},
			   {rot_i[2][0], rot_i[2][1],rot_i[2][2]}};

	float inv_rot[3][3];
	inverse_mat3b3(inv_rot, rot);

	float temp[3];
	vector3_mat3b3_multiply(inv_rot, vec_frac + 3*op, temp);
	vector3_add(temp, vec_symm, vec_symm);
    }

    // Take mean over total operations
    for(int i = 0; i < 3; i++)
	vec_symm[i] /= Z;
    
    // Apply symmetry operation to regenerate the entire vec
    for(int i = 0; i < Z; i++)
    {
	int (*rot)[3] = rotations[i];
	vector3_intmat3b3_multiply(rot, vec_symm, vec_frac + 3*i);
    }

    // Convert back to cartesian space
    for(int i = 0; i < dim; i++)
	vector3_mat3b3_multiply(lattice_t, vec_frac + 3*i, vec + 3*i);
}

void test_symmetrize_vector()
{
    int dim = 2;
    float vec[] = {5, 4, 6, -2, 1, -6};
    float lattice[3][3] = {{5, 0, 0}, {0.5, 3, 0}, {1, 0, 2}};
    for(int i = 0; i < dim; i++)
	printf("%f %f %f\n", vec[0 + 3*i], vec[1 + 3*i], vec[2 + 3*i]);

    symmetrize_vector(vec, dim, lattice, 2);
    for(int i = 0; i < dim; i++)
	printf("%f %f %f\n", vec[0 + 3*i], vec[1 + 3*i], vec[2 + 3*i]);
}

void test_symmetrize_matrix()
{
    float mat[] = { 1.02, 0, 0, 0,  1, 0, 0, 0,  1,
		    -1, 0, 0, 0, -1.05, 0, 0, 0, -1};
    int dim = 2;
    int spg = 2;

   symmetrize_matrix(mat, dim, spg);

    for(int i = 0; i < dim; i++)
    for(int j = 0; j < 3; j++)
	printf("%f, %f, %f \n",
	       mat[9*i + 0 + 3*j],
	       mat[9*i + 1 + 3*j],
	       mat[9*i + 2 + 3*j]
	       );
}

void test_rot_quat_conversion()
{
    float mat[9] = {-1, 0, 0, 0, -1, 0, 0, 0, -1};
    float quat[4];
    int inv;
    for(int j = 0; j < 3; j++)
	printf("%f, %f, %f \n",
	       mat[0 + 3*j],
	       mat[1 + 3*j],
	       mat[2 + 3*j]
	       );
    
    matrix2quaternion(mat, quat, &inv);
    printf("quat = %f %f %f %f %d\n", quat[0], quat[1], quat[2], quat[3], inv);

    quaternion2matrix(quat, mat, inv);
    for(int j = 0; j < 3; j++)
	printf("%f, %f, %f \n",
	       mat[0 + 3*j],
	       mat[1 + 3*j],
	       mat[2 + 3*j]
	       );
 }

int main()
{
    test_symmetrize_vector();
    test_symmetrize_matrix();
    //test_rot_quat_conversion();
    test_symmetrize_state();
}

/*
  Symmetrizes a matrix wrt symmetry operations of a spg.
  Useful for rotation matrices.
  mat, size -> input square matrix and its dimension. Shape = (dim, 3, 3)
  lattice   -> lattice vectors in row major form.
  spg       -> spacegroup which defines the symm operations.

  NOTE : Average rotation cannot be calculated using simple mean.
  Assumption is that elements of the input matrices are very close to
  each other.
*/

void symmetrize_matrix(float *mat, int dim, int spg)
{
    // Get symm operations
    double translations[192][3];
    int rotations[192][3][3];
    int hall_number = hall_number_from_spg(spg);
    int Z = spg_get_symmetry_from_database(rotations,
					   translations,
					   hall_number);
    assert(dim == Z);

    // Matrix to compute symm average
    float mat_symm[3][3] = {0};
    float temp[3][3] = {0};
    for(int op = 0; op < Z; op++)
    {
	int (*rot_i)[3] = rotations[op];
	float rot[3][3] = {{rot_i[0][0], rot_i[0][1], rot_i[0][2]},
			   {rot_i[1][0], rot_i[1][1], rot_i[1][2]},
			   {rot_i[2][0], rot_i[2][1], rot_i[2][2]}};

	float inv_rot[3][3]; 
	inverse_mat3b3(inv_rot, rot);

	float (*mat_i)[3] = (float (*)[3]) (mat + 3*3*op);
	mat3b3_mat3b3_multiply(inv_rot, mat_i, temp);
	mat3b3_add(temp, mat_symm, mat_symm);
    }

    // Take average
    for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
	mat_symm[i][j] /= Z;

    // For safety, ensure det(R) = 1
    float det = det_mat3b3(mat_symm);
    float renorm = fabsf(1.0 / cbrtf(det));
    for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
	mat_symm[i][j] *= renorm;
    
    // Apply symm operation to get full mat
    for(int op = 0; op < Z; op++)
    {
	int (*rot_i)[3] = rotations[op];
	float rot[3][3] = {{rot_i[0][0], rot_i[0][1], rot_i[0][2]},
			   {rot_i[1][0], rot_i[1][1], rot_i[1][2]},
			   {rot_i[2][0], rot_i[2][1], rot_i[2][2]}};

	float (*mat_i)[3] = (float (*)[3]) (mat + 3*3*op);
	mat3b3_mat3b3_multiply(rot, mat_symm, mat_i);
    }
}

// https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
// inv takes -1 or +1 and determines if the rotation is proper or improper
static void quaternion2matrix(float *quat, float *mat, const int inv)
{
    float s = 1.0/(quat[0]*quat[0] +
		   quat[1]*quat[1] +
		   quat[2]*quat[2] +
		   quat[3]*quat[3]);

    mat[0] = 1 - 2 * s * (quat[2]*quat[2] + quat[3]*quat[3]);
    mat[4] = 1 - 2 * s * (quat[1]*quat[1] + quat[3]*quat[3]);
    mat[8] = 1 - 2 * s * (quat[2]*quat[2] + quat[1]*quat[1]);

    mat[1] = 2 * s * (quat[1]*quat[2] - quat[3]*quat[0]);
    mat[2] = 2 * s * (quat[1]*quat[3] + quat[2]*quat[0]);

    mat[3] = 2 * s * (quat[1]*quat[2] + quat[3]*quat[0]);
    mat[5] = 2 * s * (quat[2]*quat[3] - quat[1]*quat[0]);
    
    mat[6] = 2 * s * (quat[1]*quat[3] - quat[2]*quat[0]);
    mat[7] = 2 * s * (quat[2]*quat[3] + quat[1]*quat[0]);

    for(int i = 0; i < 9; i++)
	mat[i] *= inv;
}

// inv -> +1 or -1 depending on proper/improper rotation
static void matrix2quaternion(float *rot, float *quat, int *inv)
{

    float det = rot[0+0*3]*rot[1+1*3]*rot[2+2*3]
	-rot[0+0*3]*rot[1+2*3]*rot[2+1*3]
	+rot[0+1*3]*rot[1+2*3]*rot[2+0*3]
	-rot[0+1*3]*rot[1+0*3]*rot[2+2*3]
	+rot[0+2*3]*rot[1+0*3]*rot[2+1*3]
	-rot[0+2*3]*rot[1+1*3]*rot[2+0*3];

    if(det < 0)
	*inv = -1;
    else
	*inv = +1;

    for(int i = 0; i < 9; i++)
	rot[i] *= *inv;
    
    float t;
    if (rot[2+2*3] < 0.0)
    {
        if (rot[0+0*3] > rot[1+1*3])
        {
            t = 1.0 + rot[0+0*3] - rot[1+1*3] - rot[2+2*3];
            quat[0] = rot[2+1*3]-rot[1+2*3];
            quat[1] = t;
            quat[2] = rot[1+0*3]+rot[0+1*3];
            quat[3] = rot[2+0*3]+rot[0+2*3];
        }
        else
        {
            t = 1.0 - rot[0+0*3] + rot[1+1*3] - rot[2+2*3];
            quat[0] = rot[0+2*3]-rot[2+0*3];
            quat[1] = rot[1+0*3]+rot[0+1*3];
            quat[2] = t;
            quat[3] = rot[2+1*3]+rot[1+2*3];
        }
    }
    else
    {
        if (rot[0+0*3] < -rot[1+1*3])
        {
            t = 1.0 - rot[0+0*3] - rot[1+1*3] + rot[2+2*3];
            quat[0] = rot[1+0*3]-rot[0+1*3];
            quat[1] = rot[0+2*3]+rot[2+0*3];
            quat[2] = rot[2+1*3]+rot[1+2*3];
            quat[3] = t;
        }
        else
        {
            t = 1.0 + rot[0+0*3] + rot[1+1*3] + rot[2+2*3];
            quat[0] = t;
            quat[1] = rot[2+1*3]-rot[1+2*3];
            quat[2] = rot[0+2*3]-rot[2+0*3];
            quat[3] = rot[1+0*3]-rot[0+1*3];
        }
    }
    t = 0.5/sqrt(t);
    quat[0] *= t;
    quat[1] *= t;
    quat[2] *= t;
    quat[3] *= t;
    
}

void test_symmetrize_state()
{
    double state[6 + 2*7] = {1, 0, 1, 0, 0, 1,
			     5, 4, 6,
			     1, 0, 0, 0,   
			     -2, 1, -6,
			     1, 0, 0, 0};   
    int nmol = 2;
    int invert[2] = {1, -1};
    int spg = 2;
    symmetrize_state(state, invert, nmol, spg);
    for(int i = 0; i < 7*nmol + 6; i++)
	printf("%f ", state[i]);
}

void symmetrize_state(double *state, int *invert, const int nmol, const int spg)
{
    // Collect position vectors
    float pos[nmol * 3];
    int dim = 0;

    float lattice[3][3] = {state[0],        0,       0,
			   state[1], state[2],       0,
			   state[3], state[4], state[5]};

    // each mol is represented by 7 numbers; first 3 are positions
    for(int i = 6; i < 7*nmol + 6; i += 7)
    {
	pos[3*dim + 0] = state[i + 0];
	pos[3*dim + 1] = state[i + 1];
	pos[3*dim + 2] = state[i + 2];
	dim++;
    }
    symmetrize_vector(pos, dim, lattice, spg);

    // Collect orientational quaternions and convert to rotation matrices
    float quat[4];
    float mat[nmol*3*3];
    dim = 0;
    for(int i = 9; i < 7*nmol + 6; i += 7)
    {
	quat[0] = state[i + 0];
	quat[1] = state[i + 1];
	quat[2] = state[i + 2];
	quat[3] = state[i + 3];
	quaternion2matrix(quat, mat + 3*3*dim, invert[dim]);
	dim++;
    }
    symmetrize_matrix(mat, dim, spg);

    // Convert back to quaternions
    dim = 0;
    float quat_all[nmol*4];
    for(int i = 9; i < 7*nmol + 6; i += 7)
    {
	matrix2quaternion(mat + 3*3*dim, quat_all + 4*dim, invert + dim);
	dim++;
    }

    // Fill in state vector using the results
    dim = 0;
    for(int i = 6; i < 7*nmol + 6; i += 7)
    {
	state[i + 0] = pos[3*dim + 0];
	state[i + 1] = pos[3*dim + 1];
	state[i + 2] = pos[3*dim + 2];
	dim ++;
    }
    dim = 0;
    for(int i = 9; i < 7*nmol + 6; i += 7)
    {
        state[i + 0] = quat_all[4*dim + 0];
        state[i + 1] = quat_all[4*dim + 1];
        state[i + 2] = quat_all[4*dim + 2];
        state[i + 3] = quat_all[4*dim + 3];
	dim++;
    }
    
}

