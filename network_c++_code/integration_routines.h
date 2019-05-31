#include <iostream>
//#ifndef INTEGRATION_ROUTINES_H
//#define INTEGRATION_ROUTINES_H
//using namespace std;
double euler_next(const double&, const double&, const double&);
double euler_next(const double& xzero, const double& fxzero, const double& delta_t, const double& sqroot_dt, const double& noise_amp, long *seed, int* iset_in, float* gset_in, int* iset_out, float* gset_out, long* idum2_in, long* iy_in, long* iv_in, long* idum2_out, long* iy_out, long* iv_out, int& done_init_gasdev_flag, int& done_init_ran2_flag, int& done_dump_gasdev_flag, int& done_dump_ran2_flag, const int& thread, double t, double stim_start, double stim_end, double stim_strength);
double euler_next(const double&, const double&, const double&, const double&, const double&, long*,int*, float*, int*, float*, long*, long*, long*, long*, long*, long*, int&, int&, int&, int&);
double euler_next(const double& xzero, const double& fxzero, const double& delta_t, const double& sqroot_dt, const double& noise_amp, long *seed, int* iset_in, float* gset_in, int* iset_out, float* gset_out, long* idum2_in, long* iy_in, long* iv_in, long* idum2_out, long* iy_out, long* iv_out, int& done_init_gasdev_flag, int& done_init_ran2_flag, int& done_dump_gasdev_flag, int& done_dump_ran2_flag, const int& thread);
double euler_next(const double&, const double&, const double&, const double&, const double&, long*);
double euler_next(const double& xzero, const double& fxzero, const double& delta_t, const double& sqroot_dt, const double& noise_amp, long *seed, const int& thread);
double poisson_next(const double& xzero, const double& fxzero, const double& poisson_rate, const double& rise, const double& dt, const double& sqroot_dt, long *seed);
double poisson_next(const double& xzero, const double& fxzero, const double& poisson_rate, const double& rise, const double& dt, const double& sqroot_dt, long *seed, const int& thread);
//#endif
