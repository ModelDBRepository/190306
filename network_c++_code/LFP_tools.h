#include "LFP.h"
#include "cortical_assignment.h"

using namespace std;
void get_all_LFP_points(int* pts_per_dim, double& radius, int n_per_dim);
void initialize_LFP(LFP* lfp_objs, int* pts_per_dim, double radius, int n_per_dim);
void initialize_Ca(LFP* lfp_objs, int lfp_objs_size);
void calculate_and_print_LFP_results(LFP* lfp_objs, int n_per_dim_sq, float* potentials, ostream& out);
void LFP_jobs_partition(int* start, int* end, int* displ, int* rcounts, int n_workers, int LFParray_length);
