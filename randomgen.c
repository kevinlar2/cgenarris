# include <math.h>
# include <stdint.h>
# include <stdio.h>
# include <stdlib.h>
# include <time.h>
 #include <stdbool.h>

# include "randomgen.h"

extern int *seed;
#pragma omp threadprivate(seed)
extern int *seed2;
#pragma omp threadprivate(seed2)

float uniform_dist_01 ()

/******************************************************************************/
/*
  Purpose:

    R4_UNIFORM_01 returns a unit pseudorandom R4.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      r4_uniform_01 = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

    If the initial seed is 12345, then the first three computations are

      Input     Output      R4_UNIFORM_01
      SEED      SEED

         12345   207482415  0.096616
     207482415  1790989824  0.833995
    1790989824  2035175616  0.947702

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2004

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has been updated.

    Output, float R4_UNIFORM_01, a new pseudorandom variate, strictly between
    0 and 1.
*/
{
  const int i4_huge = 2147483647;
  int k;
  float value;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R4_UNIFORM_01 - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }
  value = ( float ) ( *seed ) * 4.656612875E-10;

  return value;
}


float normal_dist_01 ()

/******************************************************************************/
/*
  Purpose:

    R4_NORMAL_01 returns a unit pseudonormal R4.

  Discussion:

    The standard normal probability distribution function (PDF) has 
    mean 0 and standard deviation 1.

    The Box-Muller method is used, which is efficient, but 
    generates two values at a time.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 June 2013

  Author:

    John Burkardt

  Parameters:

    Input/output, int *SEED, a seed for the random number generator.

    Output, float R4_NORMAL_01, a normally distributed random value.
*/
{
  float r1;
  float r2;
  const double r4_pi = 3.141592653589793;
  float x;

  r1 = uniform_dist_01 ( );
  r2 = uniform_dist_01 ( );
  x = sqrt ( - 2.0 * log ( r1 ) ) * cos ( 2.0 * r4_pi * r2 );

  return x;
}

float normal_dist_ab ( float mean, float stdev )

/******************************************************************************/
/*
  Purpose:

    R4_NORMAL_AB returns a scaled pseudonormal R4.

  Discussion:

    The normal probability distribution function (PDF) is sampled,
    with mean A and standard deviation B.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 June 2006

  Author:

    John Burkardt

  Parameters:

    Input, float A, the mean of the PDF.

    Input, float B, the standard deviation of the PDF.

    Input/output, int *SEED, a seed for the random number generator.

    Output, float R4_NORMAL_AB, a sample of the normal PDF.
*/
{
  float value;

  value = mean + stdev * normal_dist_01 (  );

  return value;
}

