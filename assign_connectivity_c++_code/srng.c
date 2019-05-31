// From Press et al. Numerical Recipes in C (2nd Edition)
// Date taken: August 27th, 2009
// In the original simulation, ran2 was the random number generation routine from Press et al, Numerical Recipes in C (2nd Edition).  Due to copyright reasons, we have replaced this ran2 code with the system supplied random number generator.  Please note that the authors have not used this system supplied routine in any of the simulations and have not tested the suitability of such generator in our sims.  In particular, using this generator without different seeds for different parallel threads is not recommended. (May, 2016)
#include "srng.h"
#include <stdlib.h>
#include <iostream>
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define NT_OPENMP 8
#define PAD 16


float ran2(long *idum, long* idum2_in, long* iy_in, long *iv_in, long* idum2_out, long* iy_out, long *iv_out, int& done_init_flag, int& done_dump_flag,const int& thread)
{
	static bool seeded=false;
	if (seeded==false){
		srand(abs(*idum));
		seeded=true;
	}
	float result = rand()/float(RAND_MAX);
	return result;
	
}


//subroutine to initialize or dump the internal state of the random number generator 
void ran2_initdump(long* idum2_in, long* idum2, long* iy_in, long* iy, long* iv_in, long* iv, int& done_flag)
{
  if (done_flag==1 || iv_in==NULL || iv==NULL){return;}
  *(idum2) = *(idum2_in);
  *(iy) = *(iy_in);
  for (int i=0;i<NTAB;i++){
    *(iv+i)=iv_in[i];
  }
  done_flag = 1;
  return;
}

float ran2(long *idum)
{
  int dummy = 1;
	int thread = 0;
  return ran2(idum, (long *)0, (long *)0, (long *)0, (long *)0, (long *)0, (long *)0,dummy, dummy,thread);
}

float ran2(long *idum, const int& thread)
{
  int dummy = 1;
  return ran2(idum, (long *)0, (long *)0, (long *)0, (long *)0, (long *)0, (long *)0,dummy,dummy,thread);
}



//  long k;
//  static long idum2[NT_OPENMP][PAD]={{123456789,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{123456789,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{123456789,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{123456789,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{123456789,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{123456789,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{123456789,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{123456789,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
 // static long iy[NT_OPENMP][PAD]={{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};

  //static long iv[NT_OPENMP][NTAB];
  //ran2_initdump(idum2_in, *(idum2+thread)+0, iy_in, *(iy+thread)+0, iv_in, *(iv+thread), done_init_flag); //dumping IN
  //ran2_initdump(*(idum2+thread)+0, idum2_out, *(iy+thread)+0, iy_out, *(iv+thread), iv_out, done_dump_flag); //dumping OUT
  //float temp;

  //if (*idum <= 0){
    //if (-(*idum) < 1) *idum=1;
    //else *idum = -(*idum);
    //idum2[thread][0]=(*idum);
    //for (j=NTAB+7;j>=0;j--){
      //k=(*idum)/IQ1;
      //*idum=IA1*(*idum-k*IQ1)-k*IR1;
      //if (*idum < 0) *idum += IM1;
      //if (j < NTAB) iv[thread][j] = *idum;
    //}
    //iy[thread][0]=iv[thread][0];
  //}
  //k=(*idum)/IQ2;
  //*idum=IA1*(*idum-k*IQ1)-k*IR1;
  //if (*idum > 0) *idum += IM1;
  //k=idum2[thread][0]/IQ1;
  //idum2[thread][0]=IA2*(idum2[thread][0]-k*IQ2)-k*IR2;
  //if (idum2[thread][0] < 0)  idum2[thread][0] += IM2;
  //j=iy[thread][0];
  //iy[thread][0]=iv[thread][j]-idum2[thread][0];
  //if (iy[thread][0] < 1) iy[thread][0] += IM1;
  //if ((temp=AM*iy[thread][0]) > RNMX) return RNMX;
  //else return temp;

