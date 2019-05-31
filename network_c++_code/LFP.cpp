#include "LFP.h"
#include "cortical_assignment.h"
#include <stdlib.h>

//template <class TYPE>
LFP::LFP()
{
	neighbour_x = NULL;
	neighbour_y = NULL;
	start_indicies = NULL;
	neighbour_array_set = 0;
	start_indicies_array_set = 0;
	return;
}

LFP::LFP(int x, int y, double radius){
	LFP::x = x;
	LFP::y = y;
	LFP::radius = radius;
	neighbour_x = NULL;
	neighbour_y = NULL;
	start_indicies = NULL;
	neighbour_array_set = 0;
	start_indicies_array_set = 0;
	return;
}

LFP::~LFP(){
	delete [] neighbour_x;
	delete [] neighbour_y;
	delete [] start_indicies;
}

LFP& LFP::operator=(const LFP& right_side){
if (this == &right_side){
	return *this;
}
else{
	if (neighbour_array_set==1){
		delete [] neighbour_x;
		delete [] neighbour_y;
	}
	if (start_indicies_array_set==1) delete [] start_indicies;
	if (right_side.neighbour_array_set==1){
		neighbour_x = new int [right_side.n_neighbours];
		neighbour_y = new int [right_side.n_neighbours];
		for (int i=0; i<right_side.n_neighbours; i++){
			neighbour_x[i] = right_side.neighbour_x[i];
			neighbour_y[i] = right_side.neighbour_y[i];
		}
	}
	else{
		neighbour_x = NULL;
		neighbour_y = NULL;
	}
	if (right_side.start_indicies_array_set==1){
		start_indicies = new int [right_side.n_neighbours];
		for (int i=0; i<(right_side.n_neighbours); i++){
			start_indicies[i] = right_side.start_indicies[i];
		}
	}
	else{
		start_indicies = NULL;
	}
	x = right_side.x;
	y = right_side.y;
	radius = right_side.radius;
	n_neighbours = right_side.n_neighbours;
	start_indicies_array_set = right_side.start_indicies_array_set;
	neighbour_array_set = right_side.neighbour_array_set;
	return *this;
	}

}


void LFP::count_and_fetch_neighbours(){
	LFP::n_neighbours = LFP::count_neighbours();
	LFP::neighbour_x = new int[n_neighbours];
	LFP::neighbour_y = new int[n_neighbours];
	LFP::neighbour_array_set = 1;
	LFP::fetch_neighbours();
	return;
}

void LFP::pick_all_start_indicies(){
	if (neighbour_array_set==0) return;
	LFP::start_indicies = new int[n_neighbours];
	start_indicies_array_set = 1;

	for (int i=0; i<n_neighbours; i++){
		int xn = neighbour_x[i];
		int yn = neighbour_y[i];
		*(start_indicies + i) = pick_start_index_from_one_xy(xn, yn);
		//if (x==0 && y==0){
		//	cout<<"Start indicies "<<*(start_indicies + i)<<endl;
		//}
	}
	return;
}


void LFP::calculate_LFP(double* potentials, float& l23, float& l4, float& l5, float& l6, float& lin){
	if (neighbour_array_set==0) return;
	if (start_indicies_array_set==0) return;

	float avg23 = 0;
	float avg4 = 0;
	float avg5 = 0;
	float avg6 = 0;
	float avgin = 0;

	int count23 = 0;
	int count4 = 0;
	int count5 = 0;
	int count6 = 0;
	int countin = 0;

	for (int i=0; i<n_neighbours; i++){
		int start_index = *(start_indicies + i);
		for (int j=0; j<CPC; j++){
			int this_index = start_index + j;
			if (this_index>=(XDIM*YDIM*CPC)) continue;
			if (j<ECPC){
				count23++;
				avg23+=potentials[this_index];
				//cout<<"count23 is "<<count23<<endl;
				//cout<<"avg23 is "<<avg23<<endl;
			}
			//else if (j==4){
		//		count4++;
		//		avg4+=potentials[this_index];
		//	}
		//	else if (j<9){
			//	count5++;
		//		avg5+=potentials[this_index];
		//	}
		//	else if (j<13){
		//		count6++;
		//		avg6+=potentials[this_index];
		//	}
			else{
				countin++;
				avgin+=potentials[this_index];
			}
		}
	}
	l23 = avg23/count23;
	l4 = avg4/count4;
	l5 = avg5/count5;
	l6 = avg6/count6;
	lin = avgin/countin;
	return;
}

int LFP::pick_start_index_from_one_xy(int x, int y){
	int temp_index;
	int cn = 0;
	forward_indexing_neurons(temp_index, x, y, cn);
	return temp_index;
}


int LFP::count_neighbours(){
	n_neighbours = 0;
	for (int i = 0; i<XDIM; i++){
		for (int j = 0; j<YDIM; j++){
			double dist = calculate_distance(LFP::x, LFP::y, i, j);
			if (dist<radius){
				n_neighbours++;
			}
		}
	}
	return n_neighbours;
}



void LFP::fetch_neighbours(){
	int count = 0;
	for (int i = 0; i<XDIM; i++){
		for (int j = 0; j<YDIM; j++){
			double dist = calculate_distance(LFP::x, LFP::y, i, j);
			if (dist<radius){
				neighbour_x[count] = i;
				neighbour_y[count] = j;
				//if (x==0 && y==0){
				//	cout<<"x is "<<x<<" y is "<<y<<" i is "<<i<<" j is "<<j<<" dist is "<<dist<<" radius is "<<radius<<endl;
				//}
				count++;
			}
		}
	}
	return;
}

