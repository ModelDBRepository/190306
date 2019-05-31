#include <iostream>
#include "integration_routines.h"
#include "gasdev.h"
#include "srng.h"
#include <iostream>
//#include <math.h>

using namespace std;
double euler_next(const double& xzero, const double& fxzero, const double& delta_t){
	//cout<<"xzero is "<<xzero<<endl;
	//cout<<"fxzero is "<<fxzero<<endl;
	//cout<<endl;
  return (xzero+fxzero*delta_t);
}
double euler_next(const double& xzero, const double& fxzero, const double& delta_t, const double& sqroot_dt, const double& noise_amp, long *seed, int* iset_in, float* gset_in, int* iset_out, float* gset_out, long* idum2_in, long* iy_in, long* iv_in, long* idum2_out, long* iy_out, long* iv_out, int& done_init_gasdev_flag, int& done_init_ran2_flag, int& done_dump_gasdev_flag, int& done_dump_ran2_flag, const int& thread, double t, double stim_start, double stim_end, double stim_strength){
	//cout<<"xzero is "<<xzero<<endl;
	//cout<<"fxzero is "<<fxzero<<endl;
	//cout<<endl;
	double ans = euler_next(xzero, fxzero, delta_t, sqroot_dt, noise_amp, seed, iset_in, gset_in, iset_out, gset_out, idum2_in, iy_in, iv_in, idum2_out, iy_out, iv_out, done_init_gasdev_flag, done_init_ran2_flag, done_dump_gasdev_flag, done_dump_ran2_flag, thread);
	if (t>stim_start && t<stim_end){
	return ans+delta_t*stim_strength;
	}
	return ans;
 }

double euler_next(const double& xzero, const double& fxzero, const double& delta_t, const double& sqroot_dt, const double& noise_amp, long *seed, int* iset_in, float* gset_in, int* iset_out, float* gset_out, long* idum2_in, long* iy_in, long* iv_in, long* idum2_out, long* iy_out, long* iv_out, int& done_init_gasdev_flag, int& done_init_ran2_flag, int& done_dump_gasdev_flag, int& done_dump_ran2_flag, const int& thread){
	//cout<<"xzero is "<<xzero<<endl;
	//cout<<"fxzero is "<<fxzero<<endl;
	//cout<<endl;

  double next_v=xzero+fxzero*delta_t;
  if(noise_amp<0){return next_v;}
  float ran_num = gasdev(seed, iset_in, gset_in, iset_out, gset_out, idum2_in, iy_in, iv_in, idum2_out, iy_out, iv_out, done_init_gasdev_flag, done_init_ran2_flag, done_dump_gasdev_flag, done_dump_ran2_flag, thread);
  //cout<<"Gaussian number is "<<ran_num<<endl;
  return next_v+noise_amp*sqroot_dt*ran_num;
  //return next_v+noise_amp*sqroot_dt*gasdev(seed, iset_in, gset_in, iset_out, gset_out, idum2_in, iy_in, iv_in, idum2_out, iy_out, iv_out, done_init_gasdev_flag, done_init_ran2_flag, done_dump_gasdev_flag, done_dump_ran2_flag);
}

double euler_next(const double& xzero, const double& fxzero, const double& delta_t, const double& sqroot_dt, const double& noise_amp, long *seed, int* iset_in, float* gset_in, int* iset_out, float* gset_out, long* idum2_in, long* iy_in, long* iv_in, long* idum2_out, long* iy_out, long* iv_out, int& done_init_gasdev_flag, int& done_init_ran2_flag, int& done_dump_gasdev_flag, int& done_dump_ran2_flag){
	return euler_next(xzero, fxzero, delta_t, sqroot_dt, noise_amp, seed, iset_in, gset_in, iset_out, gset_out, idum2_in, iy_in, iv_in, idum2_out, iy_out, iv_out, done_init_gasdev_flag, done_init_ran2_flag, done_dump_gasdev_flag, done_dump_ran2_flag, int(0));
}


double euler_next(const double& xzero, const double& fxzero, const double& delta_t, const double& sqroot_dt, const double& noise_amp, long *seed, const int& thread){
  int dummy = 1;
  return euler_next(xzero, fxzero, delta_t, sqroot_dt, noise_amp, seed, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL,dummy,dummy,dummy,dummy,thread);
}

double euler_next(const double& xzero, const double& fxzero, const double& delta_t, const double& sqroot_dt, const double& noise_amp, long *seed){
  int dummy = 1;
	int thread = 0;
  return euler_next(xzero, fxzero, delta_t, sqroot_dt, noise_amp, seed, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL,NULL,dummy,dummy,dummy,dummy,thread);
}

// Integration routine for Poisson neuron
double poisson_next(const double& xzero, const double& fxzero, const double& poisson_rate, const double& rise, const double& dt, const double& sqroot_dt, long *seed, const int& thread){
	float ran_num = ran2(seed, thread);
	double poisson_prob = poisson_rate*dt;
	double next_v=xzero+fxzero*dt;
	//cout<<"poisson xzero is "<<xzero<<endl;
	//cout<<"poisson fxzero is "<<fxzero<<endl;
	//cout<<endl;

	if (poisson_prob > ran_num){
	next_v+=rise;
	}
	return next_v;
}
	

double poisson_next(const double& xzero, const double& fxzero, const double& poisson_rate, const double& rise, const double& dt, const double& sqroot_dt, long *seed){
	int thread = 0;
	return  poisson_next(xzero, fxzero, poisson_rate, rise, dt, sqroot_dt, seed, thread);
}
	
