#include <math.h>
#include "srng.h"
#include <iostream>
#include "gasdev.h"
#include <cstdlib>
#include <cmath>
#include <limits>

#define NT_OPENMP 8
#define PAD 16

// Adapted from Numerical Recipes p. 289
// Date adapted September 3rd, 2009.
// In the original simulation, gasdev was the Box Muller routine from Press et al, Numerical Recipes in C (2nd Edition).  Due to copyright reasons, we have replaced the original code with a different implementation adapted from wikipedia (Box Muller transform).  Although both implementations are derived from the same Box Muller method of producing Gaussian distributed random variables from uniform distributed ones,  please note that the authors have not used this implementation in any of the simulations and have not tested the suitability of this implementation in our sims.  (May, 2016)


float gasdev(long *idum, int* iset_in, float* gset_in,  int* iset_out, float* gset_out, long* idum2_in, long* iy_in, long* iv_in, long* idum2_out, long* iy_out, long *iv_out, int& done_init_flag, int& done_init_ran2_flag, int& done_dump_flag, int& done_dump_ran2_flag, const int& thread){


	const float epsilon = std::numeric_limits<float>::min();
	//const double epsilon = 1e-32;

	//const double two_pi = 2.0*3.14159265358979323846;
	const float two_pi = 6.2831853;

	static float z0, z1;
	static bool generate;
	float mu = 0;
	float sigma = 1;

	generate = !generate;

	if (!generate)
	   return z1*sigma+mu;

	float u1, u2;
	do
	 {
	   u1 = ran2(idum);
	   u2 = ran2(idum);
	 }
	while (u1<=epsilon);

	z0 = sqrt(-2.0*log(u1))*cos(two_pi*u2);
	z1 = sqrt(-2.0*log(u1))*sin(two_pi*u2);
	return z0 * sigma + mu;
}


float gasdev(long *idum, int* iset_in, float* gset_in,  int* iset_out, float* gset_out, long* idum2_in, long* iy_in, long* iv_in, long* idum2_out, long* iy_out, long *iv_out, int& done_init_flag, int& done_init_ran2_flag, int& done_dump_flag, int& done_dump_ran2_flag){
 return gasdev(idum, iset_in, gset_in,  iset_out, gset_out, idum2_in, iy_in, iv_in, idum2_out, iy_out, iv_out, done_init_flag, done_init_ran2_flag, done_dump_flag, done_dump_ran2_flag, int(0));
}



//subroutine to initialize or dump the internal state of the random number generator 
void gasdev_initdump(int* iset_in, int* iset, float* gset_in, float* gset, int& done_flag)
{
  if (done_flag==1 || gset_in==NULL || gset==NULL){return;}
  *(iset) = *(iset_in);
  *(gset) = *(gset_in);
  done_flag = 1;
  return;
}

float gasdev(long *idum)
{
  int dummy = 1;
  return gasdev(idum, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,dummy,dummy,dummy,dummy,int(0));
}

float gasdev(long *idum, const int& thread)
{
  int dummy = 1;
  return gasdev(idum, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,dummy,dummy,dummy,dummy,thread);
}
