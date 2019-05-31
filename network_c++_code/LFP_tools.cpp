#include "LFP.h"
#include "cortical_assignment.h"
#include <math.h>
#include <assert.h>

void get_all_LFP_points(int* pts_per_dim, double& radius, int n_per_dim){
	assert(n_per_dim>1);
	int separation = int(floor(XDIM/float(n_per_dim-1)));
	radius = separation*WIDTH/2.0;
	for (int i=0; i<n_per_dim; i++){
		pts_per_dim[i] = i*separation;
	}
	return;
}

void initialize_LFP(LFP* lfp_objs, int* pts_per_dim, double radius, int n_per_dim){
	int count = 0;
	for (int i=0; i<n_per_dim; i++){
		for (int j=0; j<n_per_dim; j++){
		LFP temp(pts_per_dim[i], pts_per_dim[j], radius);
		cout<<pts_per_dim[i]<<"\t"<<pts_per_dim[j]<<endl;
		*(lfp_objs + count) = temp;
		(lfp_objs + count)->count_and_fetch_neighbours();
		(lfp_objs + count)->pick_all_start_indicies();
		count++;	
		}	
	}
	return;
}

void initialize_Ca(LFP* lfp_objs, int lfp_objs_size){
	int count = 0;
	for (int i=0; i<XDIM; i++){
		for (int j=0; j<YDIM; j++){
		assert(count<lfp_objs_size);
		LFP temp(i, j, 0.01);
		//cout<<pts_per_dim[i]<<"\t"<<pts_per_dim[j]<<endl;
		*(lfp_objs + count) = temp;
		(lfp_objs + count)->count_and_fetch_neighbours();
		(lfp_objs + count)->pick_all_start_indicies();
		count++;	
		}	
	}
	return;
}


void calculate_and_print_LFP_results(LFP* lfp_objs, int n_per_dim_sq, float* potentials, ostream& out){
	for (int i=0; i<n_per_dim_sq; i++){
		//float res = (lfp_objs + i)->calculate_LFP_L23(potentials);
		//out<<res<<"\t";
	}
	out<<"\n";
	return;
}

void LFP_jobs_partition(int* start, int* end, int* displ, int* rcount, int n_workers, int LFParray_length){
	int quo = int(ceil(LFParray_length/double(n_workers)));
	for (int i=0; i<n_workers; i++){
		start[i] = i*quo;
		end[i] = (i+1)*quo;
		if (start[i]>=LFParray_length){
			start[i] = LFParray_length;
			end[i] = LFParray_length;
		}
		else if (end[i]>=LFParray_length){
			end[i] = LFParray_length;
		}
		rcount[i] = end[i] - start[i];
		(i==0)?displ[i] = 0:displ[i] = displ[i-1]+rcount[i-1];
	}
	return;
}
