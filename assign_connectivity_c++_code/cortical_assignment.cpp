#include "cortical_assignment.h"
#include <math.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "srng.h"
#include <string>
#include <cstring>

using namespace std;
static long s_seed =-800094;//default -80

//			"0" => "EXT",
//			"1" => "INH"


int this_processor(int neuron_index, int num_proc){
	int m_size = XDIM*YDIM*CPC;
	int quotient = m_size/num_proc;
	assert(m_size%quotient==0);
	return int(neuron_index/quotient);
}

int this_processor(int x, int y, int cn, int num_proc){
	int index;
	forward_indexing_neurons(index, x, y, cn);
	return this_processor(index, num_proc);
}

int local_neuron_index(int neuron_index, int this_proc, int num_proc){
	int m_size = XDIM*YDIM*CPC;
	int quotient = m_size/num_proc;
	return (neuron_index-quotient*this_proc);
}
	

double calculate_distance(int x1, int y1, int x2, int y2){
	double dx = x2-x1;
	double dy = y2-y1;
	return WIDTH*sqrt(dx*dx + dy*dy);
}

void assign_distance_matrix(double**** distance_matrix){
	double dist = 0;
	for (int x1=0; x1<XDIM; x1++){
		for (int y1=0; y1<YDIM; y1++){
			for (int x2=0; x2<=x1; x2++){
				for (int y2=0; y2<=y1; y2++){
					dist = calculate_distance(x1, y1, x2, y2);
					distance_matrix[x1][y1][x2][y2] = dist;
					distance_matrix[x2][y1][x1][y2] = dist;
					distance_matrix[x1][y2][x2][y1] = dist;
					distance_matrix[x2][y2][x1][y1] = dist;
				}
			}
		}
	}
	return;
}

string replace_conn_m_with_symbol(int cnf, int cnt){
	if (cnf<ECPC && cnt<ECPC){
		return "EtoE";
	}
	else if (cnf<ECPC && cnt>=ECPC){
		return "EtoI";
	}
	else if (cnf>=ECPC && cnt<ECPC){
		return "ItoE";
	}
	else if (cnf>=ECPC && cnt>=ECPC){
		return "ItoI";
	}
	
}

int check_distance_within(int xf, int yf, int xt, int yt, double**** distance_matrix, double dist){
		if (distance_matrix[xf][yf][xt][yt] < dist){
			return 0;
		}
		else{
			return 1;
		}
}
void connect_e_e(double& conn, double distance, int xf, int yf, int cnf, int xt, int yt, int cnt, int& count){
	//if (xf==xt && yf==yt){
	//conn = 0.004*assign_connect_internal_by_probability(xf, yf, xt, yt,distance, 0, 0.45, count); // all-to-all internal connection (Aug 14th, 2013)
	//}
	//else{
	//conn = 0.004*assign_connect_by_probability(xf, yf, xt, yt,distance, 10000.00, 0.1, count); // June 11th, 2013; boost L23->L23 connectivitity by a factor of 3 (x3) GLOBAL
	
	// conn = 0.004*assign_connect_by_probability(xf, yf, xt, yt,distance, 450.00, 0.23, count); // June 11th, 2013; boost L23->L23 connectivitity by a factor of 3 (x3) LOCAL 450
	//conn = 0.004*assign_connect_by_probability(xf, yf, xt, yt,distance, 250.00, 0.48, count); // June 11th, 2013; boost L23->L23 connectivitity by a factor of 3 (x3) LOCAL 450
//conn = 0.004*assign_connect_by_probability(xf, yf, xt, yt,distance, 150.00, 0.65, count); // June 11th, 2013; boost L23->L23 connectivitity by a factor of 3 (x3) LOCAL 450
//conn = 0.004*assign_connect_by_probability(xf, yf, xt, yt,distance, 75.00, 0.85, count); // June 11th, 2013; boost L23->L23 connectivitity by a factor of 3 (x3) LOCAL 450
conn = 0.004*assign_connect_by_probability_no_prune(xf, yf, xt, yt,distance, 10000.00, 0.15, count); // June 11th, 2013; boost L23->L23 connectivitity by a factor of 3 (x3) LOCAL 450

	//revert to original (July 26th, 2013)
	//default probability 0.05
	//}
	return;
}
void connect_e_i(double& conn, double distance, int xf, int yf, int cnf, int xt, int yt, int cnt, int& count){
	//if (xf==xt && yf==yt){
	//conn = 0.004*assign_connect_internal_by_probability(xf, yf, xt, yt,distance, 0, 0.45, count); // all-to-all internal connection (Aug 14th, 2013)
	//}
	//else{
	//conn = 0.004*assign_connect_by_probability(xf, yf, xt, yt,distance, 10000.00, 0.1, count); // June 11th, 2013; boost L23->L23 connectivitity by a factor of 3 (x3)
//conn = 0.004*assign_connect_by_probability(xf, yf, xt, yt,distance, 450.00, 0.23, count); 
		//conn = 0.004*assign_connect_by_probability(xf, yf, xt, yt,distance, 250.00, 0.48, count);
//conn = 0.004*assign_connect_by_probability(xf, yf, xt, yt,distance, 150.00, 0.65, count);
//conn = 0.004*assign_connect_by_probability(xf, yf, xt, yt,distance, 75.00, 0.85, count);
conn = 0.004*assign_connect_by_probability_no_prune(xf, yf, xt, yt,distance, 75000.00, 0.15, count);

	//revert to original (July 26th, 2013)
	//default probability 0.05
	//}
	return;
}
void connect_i_e(double& conn, double distance, int xf, int yf, int cnf, int xt, int yt, int cnt, int& count){
	//if (xf==xt && yf==yt){
	//conn = 0.03*assign_connect_internal_by_probability(xf, yf, xt, yt,distance, 0, 0.35, count); // all-to-all internal connection (Aug 14th, 2013)
	//}
	//else{
//	conn = 0.03*assign_connect_by_probability(xf, yf, xt, yt,distance, 10000.00, 0.1, count); // June 11th, 2013; boost L23->L23 connectivitity by a factor of 3 (x3)
	// conn = 0.03*assign_connect_by_probability_no_prune(xf, yf, xt, yt,distance, 450.00, 0.20, count); 
	// conn = 0.03*assign_connect_by_probability_no_prune(xf, yf, xt, yt,distance, 250.00, 0.35, count); 
//conn = 0.03*assign_connect_by_probability_no_prune(xf, yf, xt, yt,distance, 150.00, 0.45, count); 
//conn = 0.03*assign_connect_by_probability_no_prune(xf, yf, xt, yt,distance, 75.00, 0.6, count); 
conn = 0.03*assign_connect_by_probability_no_prune(xf, yf, xt, yt,distance, 75000.00, 0.25, count); 

	//revert to original (July 26th, 2013)
	//default probability 0.05
//	}
	return;
}
void connect_i_i(double& conn, double distance, int xf, int yf, int cnf, int xt, int yt, int cnt, int& count){
	//if (xf==xt && yf==yt){
	//conn = 0.03*assign_connect_internal_by_probability(xf, yf, xt, yt,distance, 0, 0.35, count); // all-to-all internal connection (Aug 14th, 2013)
	//}
	//else{
	//conn = 0.03*assign_connect_by_probability(xf, yf, xt, yt,distance, 10000.00, 0.1, count); // June 11th, 2013; boost L23->L23 connectivitity by a factor of 3 (x3)
	// conn = 0.03*assign_connect_by_probability_no_prune(xf, yf, xt, yt,distance, 450.00, 0.20, count); 
	//conn = 0.03*assign_connect_by_probability_no_prune(xf, yf, xt, yt,distance, 250.00, 0.35, count); 
//conn = 0.03*assign_connect_by_probability_no_prune(xf, yf, xt, yt,distance, 150.00, 0.45, count); 
//conn = 0.03*assign_connect_by_probability_no_prune(xf, yf, xt, yt,distance, 75.00, 0.6, count);
	conn = 0.03*assign_connect_by_probability_no_prune(xf, yf, xt, yt,distance, 75000.00, 0.25, count); 
	//revert to original (July 26th, 2013)
	//default probability 0.05
	//}
	return;
}



double prep_a_b_dying_probability(double b_density, double range, double target_connection_to_b_for_one_a){
	return target_connection_to_b_for_one_a/(5*PI*range*range*b_density);
}

double assign_connect_internal_by_probability(int xf, int yf, int xt, int yt, double distance, double max_dist, double prob, int& count){
		if (xt!=xf || yt!=yf) return 0; 
		double r = ran2(&s_seed);
		prob*=impose_boundary_conditions(xt, yt);

		if (r<prob){
			count++;
			return 1+0.4*(ran2(&s_seed)-0.5);
		}
		else{
			return 0;
		}	
}


double assign_connect_by_probability(int xf, int yf, int xt, int yt, double distance, double max_dist, double prob, int& count){
		//if (distance_matrix[xf][yf][xt][yt]>max_dist) return 0;
		if (distance>max_dist) return 0;
		double r = ran2(&s_seed);
		prob*=impose_boundary_conditions(xt, yt);

		if (r<prob){
			count++;
			return 1+0.4*(ran2(&s_seed)-0.5);
		}
		else{
			return 0;
		}	
}
double assign_connect_by_probability_no_prune(int xf, int yf, int xt, int yt, double distance, double max_dist, double prob, int& count){
		//if (distance_matrix[xf][yf][xt][yt]>max_dist) return 0;
		if (distance>max_dist) return 0;
		double r = ran2(&s_seed);
		//prob*=impose_boundary_conditions(xt, yt);

		if (r<prob){
			count++;
			return 1+0.4*(ran2(&s_seed)-0.5);
		}
		else{
			return 0;
		}	
}


double assign_connect_by_dying_probability(int xf, int yf, int xt, int yt, double distance, double range, double prob, int& count){
	//double dist = distance_matrix[xf][yf][xt][yt];
	double r = ran2(&s_seed);
	prob*=impose_boundary_conditions(xt, yt);
	if (distance<range){
		if (r<prob){
			count++;
			return 1+0.4*(ran2(&s_seed)-0.5);
		}
	}
	else if (r<prob*exp(-(distance/range)/range)){
		count++;
		return 1;
	}


	return 0;
}

double impose_boundary_conditions(int x, int y){
// a la Anderson!
double reduction = 1;
int y1 = YDIM -1;
int x1 = XDIM -1;

if ((y<CRITDEX) && (x>=CRITDEX) && (x<=(x1-CRITDEX))){
	//reduction*=0.625;
	reduction*=0.25;

}
if ((y>(y1-CRITDEX)) && (x>=CRITDEX) && (x<=(x1-CRITDEX))){
	//reduction*=0.625;
	reduction*=0.25;

}
if ((x<CRITDEX) && (y>=CRITDEX) && (y<=(y1-CRITDEX))){
	//reduction*=0.625;
	reduction*=0.25;

}
if ((x>(x1-CRITDEX)) && (y>=CRITDEX) && (y<=(y1-CRITDEX))){
	//reduction*=0.625;
	reduction*=0.25;

}
if ((y<CRITDEX) && (x<CRITDEX)){
//	reduction*=0.375;
	reduction*=0.0625;

}
if ((y<CRITDEX) && (x>(x1-CRITDEX))){
	//reduction*=0.375;
	reduction*=0.0625;

}
if((y>(y1-CRITDEX)) && (x<CRITDEX)){
//	reduction*=0.375;
	reduction*=0.0625;

}
//Added May 11th,2014
if((y>(y1-CRITDEX)) && (x>(x1-CRITDEX))){
	reduction*=0.0625;
}


return reduction;
}


int calculate_total_number_of_cells(){
	return XDIM*YDIM*CPC;
}

double calculate_density(){
	int num = calculate_total_number_of_cells();
	double area = XDIM*YDIM*WIDTH*WIDTH;
	return num/area;
}


int calculate_total_e(){
	return XDIM*YDIM*ECPC;
}

double calculate_e_density(){
	int num = calculate_total_e();
	double area = XDIM*YDIM*WIDTH*WIDTH;
	return num/area;
}


int calculate_total_i(){
	return XDIM*YDIM;
}

double calculate_i_density(){
	int num = calculate_total_i();
	double area = XDIM*YDIM*WIDTH*WIDTH;
	return num/area;
}


void forward_indexing_neurons(int& indexed_number, int x, int y, int cn){
	assert(cn<CPC);
	assert(XDIM==YDIM);
	assert(x<XDIM);
	indexed_number = CPC*(x + y*YDIM) + cn;
	return;
}

void reverse_indexing_neurons(int indexed_number, int& x, int& y, int& cn){
	cn = indexed_number%CPC;
	int temp = (indexed_number - cn)/CPC; // should be a whole number
	y = temp/YDIM;
	x = temp%YDIM;
	return;
}

void custom_convert_lcc_to_cn(int layer, int cell_type, int cell_number, int& cn){
	cn = -1;
	//cout<<layer<<endl;
	//cout<<cell_type<<endl;
	//cout<<cell_number<<endl;
	//cout<<endl;
if (layer==0 && cell_type==0 && cell_number==0){ // excit
	cn = 0;
}
else if (layer==0 && cell_type==0 && cell_number==1){ // excit
	cn = 1;
}
else if (layer==0 && cell_type==0 && cell_number==2){ // excit
	cn = 2;
}
else if (layer==0 && cell_type==0 && cell_number==3){ // excit
	cn = 3;
}
else if (layer==0 && cell_type==0 && cell_number==4){ // excit

	cn = 4;
}
else if (layer==0 && cell_type==0 && cell_number==5){ // excit
	cn = 5;
}
else if (layer==0 && cell_type==0 && cell_number==6){ // excit
	cn = 6;
}
else if (layer==0 && cell_type==0 && cell_number==7){ // excit
	cn = 7;
}
else if (layer==0 && cell_type==0 && cell_number==8){ // excit
	cn = 8;
}
else if (layer==0 && cell_type==0 && cell_number==9){ // excit
	cn = 9;
}
else if (layer==0 && cell_type==0 && cell_number==10){ // excit
	cn = 10;
}
else if (layer==0 && cell_type==0 && cell_number==11){ // excit
	cn = 11;
}
else if (layer==0 && cell_type==0 && cell_number==12){ // inhib
	cn = 12;
}
else if (layer==0 && cell_type==1 && cell_number==0){ // inhib
	cn = 13;
}
else if (layer==0 && cell_type==1 && cell_number==1){ // inhib
	cn = 14;
}
else if (layer==0 && cell_type==1 && cell_number==2){ // inhib
	cn = 15;
}

assert(cn!=-1);
return;
}

void custom_convert_cn_to_lcc(int& layer, int& cell_type, int& cell_number, int cn){
	switch(cn){
		case 0:
			layer = 0; cell_type = 0; cell_number = 0;
			break;
		case 1:
			layer = 0; cell_type = 0; cell_number = 1;
			break;
		case 2:
			layer = 0; cell_type = 0; cell_number = 2;
			break;
		case 3:
			layer = 0; cell_type = 0; cell_number = 3;
			break;
		case 4:
			layer = 0; cell_type = 0; cell_number = 4;
			break;
		case 5:
			layer = 0; cell_type = 0; cell_number = 5;
			break;
		case 6:
			layer = 0; cell_type = 0; cell_number = 6;
			break;
		case 7:
			layer = 0; cell_type = 0; cell_number = 7;
			break;
		case 8:
			layer = 0; cell_type = 0; cell_number = 8;
			break;
		case 9:
			layer = 0; cell_type = 0; cell_number = 9;
			break;
		case 10:
			layer = 0; cell_type = 0; cell_number = 10;
			break;
		case 11:
			layer = 0; cell_type = 0; cell_number = 11;
			break;
		case 12:
			layer = 0; cell_type = 0; cell_number = 12;
			break;
		case 13:
			layer = 0; cell_type = 1; cell_number = 0;
			break;
		case 14:
			layer = 0; cell_type = 1; cell_number = 1;
			break;
		case 15:
			layer = 0; cell_type = 1; cell_number = 2;
			break;

			default:
			break;
	}
	return;
}

void print_num_neurons(ostream& out){
	int num_neurons = XDIM*YDIM*CPC;
	out.write((char*)&(num_neurons), sizeof(int));
	return;
}

void read_num_neurons(istream& in, int* num_neurons){
	in.read((char*)num_neurons, sizeof(int));
	return;
}

void print_seed(ostream& out, int this_proc){
	long seed = -(this_proc+1)*8396;//default 100
	//long seed = -8396; // Same seed for all MPI processes (wrong!)
	out.write((char*)&seed, sizeof(long));
	return;
}

void read_seed(istream& in, long* seed){
	in.read((char*)seed, sizeof(long));
	return;
}

void print_connectivity_matrix_testing(double** conn_m, ostream& out, int this_proc, int num_proc){
	assert(num_proc==2);
	assert(this_proc<2);
	int min_i = 2*this_proc;
	int max_i = 2*(this_proc +1);
	int num_neurons = 4;
	assert(num_neurons==4);
	out.write((char*)&num_neurons, sizeof(int));
	print_seed(out, this_proc);
	for (int i = min_i; i<max_i; i++){
		for (int j = 0; j<num_neurons; j++){
		out.write((char*)&(conn_m[i][j]), sizeof(double));
			}
		//out<<"\n";
		}
	return;
}
void print_connectivity_matrix_custom(double****** conn_m, ostream& out){
	int m_size = XDIM*YDIM*CPC;
	int xf, yf, cnf, xt, yt, cnt;
	int typef, numf, layerf, typet, numt, layert;
	for (size_t i=0; i<m_size; i++){
		for (size_t j=0; j<m_size; j++){
				reverse_indexing_neurons(i, xf, yf, cnf);
				reverse_indexing_neurons(j, xt, yt, cnt);
				print_connectivity_matrix_custom_per_row(i, xf, yf, cnf, j, xt, yt, cnt, conn_m[xf][yf][cnf][xt][yt][cnt], out);
				//custom_convert_cn_to_lcc(layerf, typef, numf, cnf);
				//custom_convert_cn_to_lcc(layert, typet, numt, cnt);
				//if (conn_m[xf][yf][cnf][xt][yt][cnt]>1e-50) out<<xf<<"\t"<<yf<<"\t"<<layerf<<"\t"<<typef<<"\t"<<numf<<"\t"<<xt<<"\t"<<yt<<"\t"<<layert<<"\t"<<typet<<"\t"<<numt<<"\t"<<conn_m[xf][yf][cnf][xt][yt][cnt]<<endl;		
		}
	}
	return;
}
template <class	TYPE> void print_connectivity_matrix_custom_per_row(int i, int xf, int yf, int cnf, int j, int xt, int yt, int cnt, TYPE conn, ostream&out){
	//if (conn<1e-50) return;
	int typef, numf, layerf, typet, numt, layert;
	custom_convert_cn_to_lcc(layerf, typef, numf, cnf);
	custom_convert_cn_to_lcc(layert, typet, numt, cnt);
	out<<xf<<"\t"<<yf<<"\t"<<layerf<<"\t"<<typef<<"\t"<<numf<<"\t"<<xt<<"\t"<<yt<<"\t"<<layert<<"\t"<<typet<<"\t"<<numt<<"\t"<<conn<<endl;	
	return;
}
template void print_connectivity_matrix_custom_per_row<double>(int i, int xf, int yf, int cnf, int j, int xt, int yt, int cnt, double conn, ostream&out);
template void print_connectivity_matrix_custom_per_row<string>(int i, int xf, int yf, int cnf, int j, int xt, int yt, int cnt, string conn, ostream&out);

template <class	TYPE> void print_connectivity_matrix_custom_per_row(int xf, int yf, int cnf, int xt, int yt, int cnt, TYPE conn, ostream&out){
	int i=0;
	int j=0;
	print_connectivity_matrix_custom_per_row(i, xf, yf, cnf, j, xt, yt, cnt, conn, out);
	return;
}
template void print_connectivity_matrix_custom_per_row<double>(int xf, int yf, int cnf, int xt, int yt, int cnt, double conn, ostream&out);
template void print_connectivity_matrix_custom_per_row<string>(int xf, int yf, int cnf, int xt, int yt, int cnt, string conn, ostream&out);
	
template <class	TYPE> void print_connectivity_matrix_custom_per_row_sorted(int xf, int yf, int cnf, int xt, int yt, int cnt, TYPE conn, int syn_count[], int this_proc, int num_proc, ostream&out){
	//if (conn<1e-50) return;
	int neuron_numf, neuron_numt;
	forward_indexing_neurons(neuron_numf, xf, yf, cnf);
	forward_indexing_neurons(neuron_numt, xt, yt, cnt);
	int local_neuron_numt = local_neuron_index(neuron_numt, this_proc, num_proc);
	syn_count[neuron_numt]++;
	out << neuron_numf<<"\t"<<local_neuron_numt<<"\t"<<conn<<endl;
	return;
}
template void print_connectivity_matrix_custom_per_row_sorted<double>(int xf, int yf, int cnf, int xt, int yt, int cnt, double conn, int syn_count[], int this_proc, int num_proc, ostream&out);
template void print_connectivity_matrix_custom_per_row_sorted<string>(int xf, int yf, int cnf, int xt, int yt, int cnt, string conn, int syn_count[], int this_proc, int num_proc, ostream&out);

void print_connectivity_matrix(double****** conn_m, ostream& out, int this_proc, int num_proc){	
	int m_size = XDIM*YDIM*CPC;
	int quotient = m_size/num_proc;
	int xf, yf, cnf, xt, yt, cnt;
	int min_i = quotient*this_proc;
	int max_i = quotient*(this_proc +1);
	assert(m_size%num_proc==0);
	assert(this_proc<num_proc);
	assert(this_proc>=0);
	print_num_neurons(out);
	print_seed(out, this_proc);

	for (size_t i=min_i; i<max_i; i++){
		for (size_t j=0; j<m_size; j++){
			reverse_indexing_neurons(i, xf, yf, cnf);
			reverse_indexing_neurons(j, xt, yt, cnt);
			//out<<conn_m[xf][yf][cnf][xt][yt][cnt]<<"\t";
			out.write((char*)&(conn_m[xf][yf][cnf][xt][yt][cnt]), sizeof(double));
		}
		//out<<"\n";
	}
	return;
}
void print_connectivity_matrix(ostream& out, int this_proc, int num_proc){	
	int m_size = XDIM*YDIM*CPC;
	int quotient = m_size/num_proc;
	int xf, yf, cnf, xt, yt, cnt;
	int min_i = quotient*this_proc;
	int max_i = quotient*(this_proc +1);
	assert(m_size%num_proc==0);
	assert(this_proc<num_proc);
	assert(this_proc>=0);
	print_num_neurons(out);
	print_seed(out, this_proc);

	return;
}

void read_connectivity_matrix(double **syn, istream& in, int num_neurons_this_proc, int m_size, int num_proc, int this_proc){
	assert(num_proc*num_neurons_this_proc==m_size);
	size_t start_neuron = this_proc*num_neurons_this_proc;
	size_t end_neuron = (this_proc+1)*num_neurons_this_proc;

	for (size_t i=start_neuron; i<end_neuron; i++){
		for (size_t j=0; j<m_size; j++){
			in.read((char*)&(syn[i-start_neuron][j]), sizeof(double));
		}
	}
	return;
}

void read_connectivity_matrix_custom(double **syn_sorted, int **syn_location, int *syn_counter, int* syn_size, istream& in_conn){
	int n_indexf, n_indext_local;
	double conn;
	while (in_conn>>n_indexf>>n_indext_local>>conn){
		//in_conn>>n_indexf>>n_indext_local>>conn;
		int count = syn_counter[n_indext_local];
		assert(count<syn_size[n_indext_local]);
		syn_sorted[n_indext_local][count] = conn;
		syn_location[n_indext_local][count] = n_indexf;
		syn_counter[n_indext_local]++;
		//cout<<"count is "<<count<<" size is "<<syn_size[n_indext_local]<<" n_indexf is "<<n_indexf<<" n_indext_local is "<<n_indext_local<<endl;
	}
	return;
}

void read_syn_size(int *syn_size, int size_of_array, istream& in_size){
	int count = 0;
	for (int i=0; i<size_of_array; i++){
		int size;
		in_size>>size;
		syn_size[i] = size;
		count++;
	}
assert(count==size_of_array);
return;
}
	

void read_connectivity_matrix_custom_per_row(int& neuron_numf, int& neuron_numt, double& conn_value, istream& in){
	int xf, yf, layerf, typef, numf, cnf;
	int xt, yt, layert, typet, numt, cnt;
	in>>xf>>yf>>layerf>>typef>>numf>>xt>>yt>>layert>>typet>>numt>>conn_value;
	custom_convert_lcc_to_cn(layerf, typef, numf, cnf);
	custom_convert_lcc_to_cn(layert, typet, numt, cnt);
	forward_indexing_neurons(neuron_numf, xf, yf, cnf);
	forward_indexing_neurons(neuron_numt, xt, yt, cnt);

	return;
}

void assign_neurons_testing(Fmap* var_vec, Fmap* par_vec, int num_neurons_this_proc, int total_neurons, int num_proc, int this_proc){
assign_neurons_testing(var_vec, par_vec, num_neurons_this_proc, total_neurons, num_proc, this_proc, -0.08, -0.07, 0.11);
return;
}

void assign_neurons(Fmap* var_vec, Fmap* par_vec, int num_neurons_this_proc, int total_neurons, int num_proc, int this_proc){
assign_neurons(var_vec, par_vec, num_neurons_this_proc, total_neurons, num_proc, this_proc, -0.08, -0.07, 0.11);
return;
}


void assign_neurons_testing(Fmap* var_vec, Fmap* par_vec, int num_neurons_this_proc, int total_neurons, int num_proc, int this_proc, double rs_iext, double ib_iext, double fs_iext){
	assert(num_proc==2);
	assert(this_proc<num_proc);
	assert(total_neurons==4);
	assert(num_neurons_this_proc==2);

	size_t start_neuron = this_proc*num_neurons_this_proc;
	size_t end_neuron = (this_proc+1)*num_neurons_this_proc;

	for (size_t i = start_neuron; i<end_neuron; i++){

			if (i==0){
			*(par_vec+i-start_neuron) = regular_spiking_excitatory_par(rs_iext);
			*(var_vec+i-start_neuron) = anderson_Ca_var();
		}
		else if (i==1){
			*(par_vec+i-start_neuron) = regular_spiking_inhibitory_par(rs_iext);
			*(var_vec+i-start_neuron) = anderson_Ca_var();

		}
		else if (i==2){
			*(par_vec+i-start_neuron) = intrinsic_bursting_excitatory_par(ib_iext);
			*(var_vec+i-start_neuron) = anderson_Ca_var();

		}
		else if (i==3){
			*(par_vec+i-start_neuron) = fast_spiking_inhibitory_par(fs_iext);
			*(var_vec+i-start_neuron) = anderson_regular_var();

		}
	}
	return;
}

void assign_neurons(Fmap* var_vec, Fmap* par_vec, int num_neurons_this_proc, int total_neurons, int num_proc, int this_proc, double rs_iext, double ib_iext, double fs_iext){
assign_neurons(var_vec, par_vec, num_neurons_this_proc, total_neurons, num_proc, this_proc, rs_iext, ib_iext, fs_iext, 20, 15);
}
//void assign_neurons_ahmed_modified(Fmap* var_vec, Fmap* par_vec, int num_neurons_this_proc, int total_neurons, int num_proc, int this_proc, double rs_iext, double ib_iext, double fs_iext){
//assign_neurons_ahmed_modified(var_vec, par_vec, num_neurons_this_proc, total_neurons, num_proc, this_proc, rs_iext, ib_iext, fs_iext, 20, 15);
//}

void assign_neurons(Fmap* var_vec, Fmap* par_vec, int num_neurons_this_proc, int total_neurons, int num_proc, int this_proc, double rs_iext, double ib_iext, double fs_iext, double delay, double delay_var){
	assert(num_proc*num_neurons_this_proc==total_neurons);
	size_t start_neuron = this_proc*num_neurons_this_proc;
	size_t end_neuron = (this_proc+1)*num_neurons_this_proc;
	for (size_t i=start_neuron; i<end_neuron; i++){
		int x, y, cn;
		reverse_indexing_neurons(i, x, y, cn);
		assert(CPC>ECPC);
		if (cn<ECPC){
		//	*(par_vec+i-start_neuron) = if_excitatory_par(rs_iext, delay, delay_var);
		//		*(var_vec+i-start_neuron) = if_var();
				*(par_vec+i-start_neuron) = fast_spiking_inhibitory_par(rs_iext, delay, delay_var);
				*(var_vec+i-start_neuron) = anderson_regular_var();

		}
		else if (cn<CPC){
			//	*(par_vec+i-start_neuron) = if_inhibitory_par(fs_iext, delay, delay_var);
			//	*(var_vec+i-start_neuron) = if_var();	
					*(par_vec+i-start_neuron) = fast_spiking_inhibitory_par(rs_iext, delay, delay_var);
				*(var_vec+i-start_neuron) = anderson_regular_var();

		}
		}
	return;
}

void assign_neurons_ahmed_modified(Fmap* var_vec, Fmap* par_vec, int num_neurons_this_proc, int total_neurons, int num_proc, int this_proc, double rs_iext, double ib_iext, double fs_iext, double delay, double delay_var){
	assert(num_proc*num_neurons_this_proc==total_neurons);
	size_t start_neuron = this_proc*num_neurons_this_proc;
	size_t end_neuron = (this_proc+1)*num_neurons_this_proc;
	for (size_t i=start_neuron; i<end_neuron; i++){
		int x, y, cn;
		reverse_indexing_neurons(i, x, y, cn);
			assert(CPC>ECPC);
		if (cn<ECPC){
			*(par_vec+i-start_neuron) = regular_spiking_excitatory_ahmed_modified_par(rs_iext, delay, delay_var);
				*(var_vec+i-start_neuron) = ahmed_modified_var();
		}
		else if (cn<CPC){
				*(par_vec+i-start_neuron) = fast_spiking_inhibitory_ahmed_modified_par(fs_iext, delay, delay_var);
			*(var_vec+i-start_neuron) = ahmed_modified_var();

		}
		}
	return;
}


void append_file_with_proc_number(string file_name_str, int proc_number, char*& result_file_name, string name_type){

	char rank_char[10];
	sprintf(rank_char, "%d", proc_number);
	std::string rank_string(rank_char);
	rank_string = name_type + rank_string;

	file_name_str=file_name_str.insert(file_name_str.length()-4, rank_string);
	//assert(result_file_name==NULL);
	result_file_name = new char[file_name_str.size()+1];
	strcpy(result_file_name, file_name_str.c_str());
	cout<<"result_file_name is "<<result_file_name<<endl;
	return;
}

Fmap if_excitatory_par(double iext, double delay, double delay_var){
	Fmap spec;
	spec["Vtheta"]=20.0;
	spec["external_current"]=iext;
	spec["membrane_tau"]=20.0;
	spec["Vreset"]=10.0;
	spec["refractory_time"]=1.5;
	spec["tau_i"]=0.5;
	spec["synapse_i_max"]=1; //0.01 is dt to accomodate the delta function for instanteous rise for IF. Hack.
	spec["Vmax"]=30;
	spec["reversal_potential_e"]=10;//Not used here
	spec["reversal_potential_i"]=10;//Not used here
	spec["synapse_delay"] = delay + delay_var*(ran2(&s_seed)-0.5); //Aug 19th, 2013

	return spec;
}
Fmap if_inhibitory_par(double iext, double delay, double delay_var){
	Fmap spec;
	spec = if_excitatory_par(iext, delay, delay_var);
	spec["ei"]=-1;
	return spec;
}

Fmap if_var(){
	Fmap spec;
	spec["membrane_potential"]=10;
	return spec;
}

Fmap regular_spiking_excitatory_par(double iext){
	return regular_spiking_excitatory_par(iext, 20, 15);
}

Fmap regular_spiking_excitatory_ahmed_modified_par(double iext){
	return regular_spiking_excitatory_ahmed_modified_par(iext, 20, 15);
}


Fmap regular_spiking_excitatory_par(double iext, double delay, double delay_var){// becomes RS
	Fmap spec;
	spec["aA"] = 0.02;
	spec["vhalfA"] = -20.0;
	spec["aX"] = 0.2;//(2);
	spec["vhalfX"] = -45.0;
	spec["aB"]= -0.1;
	spec["vhalfB"] = -70.0;
	spec["aW"] = 0.055;
	spec["vhalfW"] = -35.0;
	spec["am"] = 0.065;
	spec["vhalfm"] = -31.0;
	spec["lambda"] = 0.08;
	spec["tauX"] = 3.0;
	spec["tauB"] = 1.0;
	spec["VNa"] = 55.0;
	spec["VK"] = -72.0;
	spec["VL"] = -50.0;
	spec["VCa"] = 124.0;
	spec["gNa"] = 120.0;//*(1+0.1*(ran2(&s_seed)-0.5));
	spec["gK"] = 15.0;
	spec["gL"] = 0.3;
	spec["gA"] = 12.5;
	spec["gKCa"] = 3.5;
	spec["Kc"] = 2.0;
	spec["Kd"] = 0.5;
	spec["PCa"] = 0.12;
	spec["C0"] = 2.0;
	spec["Rp"] = 3.0;
	spec["Kp"] = 0.5;
	spec["b"] = 0.1;
	spec["f"] = 0.3;
	spec["B0"] = 60.0;
	spec["D"] = 2e-9;//6e-9;
	spec["taud"] = 3;
	spec["tau_i"] = 3;//*(1+0.1*(ran2(&s_seed)-0.5));;
	//spec["synapse_delay"] = 8.0*(1+0.5*(ran2(&s_seed)-0.5));
	spec["synapse_delay"] = delay + delay_var*(ran2(&s_seed)-0.5); //Aug 19th, 2013
	//spec["Esyn"] = -10.0;
	spec["reversal_potential_e"] = -10.0;
	spec["reversal_potential_i"] = -72.0;
	spec["synapse_i_max"] = 0.55;  // default 0.55
	spec["poisson_rate"] = 0.0;
	spec["external_current"] = iext;//-0.09;//0.2; //-0.09 is the threshold for firing
	spec["nt"] = 1; // Anderson
	spec["ei"] = 1; // Excitatory
	return spec;
}

Fmap regular_spiking_excitatory_ahmed_modified_par(double iext, double delay, double delay_var){// becomes RS

	Fmap spec;
	spec["a"] = 0.02; // inverse time constant of adaptation current
	spec["k"] = 0.1;
	spec["k_after"] = 0.2;
	spec["b"] = -0.001;

	spec["Vreset"] = -65.0;
	spec["Vtshift"] = 0.0;
	spec["Vmax"] = 30.0;
	spec["Vchangek"] = -20.0;
	spec["Vr"]=-65.0;

	spec["Rreset"] = 2;
	spec["reversal_potential_e"] = -10.0;
	spec["reversal_potential_i"] = -72.0;
	spec["synapse_i_max"] = 0.55;  // default 0.55
	spec["poisson_rate"] = 0.0;
	spec["ei"] = 1; // excitatory
	spec["d_step"]= 1.0;
	spec["external_current"] = iext;
	spec["synapse_delay"] = delay + delay_var*(ran2(&s_seed)-0.5); //Aug 19th, 2013
	spec["tau_i"]= 3;
	return spec;
}


Fmap fast_spiking_inhibitory_par(double iext){
	return fast_spiking_inhibitory_par(iext, 20, 15);
}

Fmap fast_spiking_inhibitory_ahmed_modified_par(double iext){
	return fast_spiking_inhibitory_ahmed_modified_par(iext, 20, 15);
}


Fmap fast_spiking_inhibitory_par(double iext, double delay, double delay_var){//becomes FS
	Fmap spec;
	spec["aA"] = 0.02;
	spec["vhalfA"] = -20;
	spec["aX"] = 0.2;//(2);
	spec["vhalfX"] = -45;
	spec["aB"] =  -0.1;
	spec["vhalfB"] = -70;
	spec["aW"] = 0.055; 
	spec["vhalfW"] = -35;
	spec["am"] = 0.065;
	spec["vhalfm"] = -31;
	spec["lambda"] = 0.08;
	spec["tauX"] = 12.0;
	spec["tauB"] = 10.0;
	spec["VNa"] = 55;
	spec["VK"] = -72;
	spec["VL"] = -53;
	spec["VCa"] = 124;
	spec["gNa"] = 120;//*(1+0.1*(ran2(&s_seed)-0.5));
	spec["gK"] = 15;
	spec["gL"] = 0.3;
	spec["gA"] = 9;
	spec["gKCa"] = 0;
	spec["Kc"] = 0;
	spec["Kd"] = 0;
	spec["PCa"] = 0;
	spec["C0"] = 0;
	spec["Rp"] = 0;
	spec["Kp"] = 0;
	spec["b"] = 0;
	spec["f"] = 0;
	spec["B0"] = 0;
	spec["D"] = 0;
	spec["taud"] = 3;
	spec["tauo"] = 0.5;
	spec["tau_i"] = 5.0;//used to be 3//*(1+0.1*(ran2(&s_seed)-0.5));
	spec["tauDelay"] = 2;
	//spec["synapse_delay"] =6.0*(1+0.5*(ran2(&s_seed)-0.5));
	//spec["synapse_delay"] =20.0*(1+0.5*(ran2(&s_seed)-0.5)); // Aug 19th, 2013
	spec["synapse_delay"] = delay + delay_var*(ran2(&s_seed)-0.5); //Aug 19th, 2013
	spec["reversal_potential_e"] = -10;
	spec["reversal_potential_i"] = -72;
	//Hack
	//spec["synapse_i_max"] = 0.55;  // default 0.55
	spec["synapse_i_max"]=1;
	//Hack to accomodate IF
	spec["refractory_time"]=2;
	spec["poisson_rate"] = 0;
	spec["external_current"] = iext;//0.2;
	spec["nt"] = -1; // sans Ca
	spec["ei"] = -1; // inh

			return spec;

}

Fmap fast_spiking_inhibitory_ahmed_modified_par(double iext, double delay, double delay_var){//becomes FS
	Fmap spec;
	spec["a"] = 0.02; // inverse time constant of adaptation current
	spec["k"] = 0.1;
	spec["k_after"] = 0.1;
	spec["b"] = -0.00;

	spec["Vreset"] = -65.0;
	spec["Vtshift"] = 0.0;
	spec["Vmax"] = 30.0;
	spec["Vchangek"] = -20.0;
	spec["Vr"]=-65.0;

	spec["Rreset"] = 2;
	spec["reversal_potential_e"] = -10.0;
	spec["reversal_potential_i"] = -72.0;
	spec["synapse_i_max"] = 0.55;  // default 0.55
	spec["poisson_rate"] = 0.0;
	spec["ei"] = -1; // excitatory
	spec["d_step"]= 0.0;
	spec["external_current"] = iext;
	spec["synapse_delay"] = delay + delay_var*(ran2(&s_seed)-0.5); //Aug 19th, 2013
	spec["tau_i"]=5.0;
	return spec;
}


Fmap regular_spiking_inhibitory_par(double iext){
	return regular_spiking_inhibitory_par(iext, 20, 15);
}
Fmap regular_spiking_inhibitory_ahmed_modified_par(double iext){
	return regular_spiking_inhibitory_ahmed_modified_par(iext, 20, 15);
}

Fmap regular_spiking_inhibitory_par(double iext, double delay, double delay_var){
	Fmap spec = regular_spiking_excitatory_par(iext, delay, delay_var);
	spec["ei"] = -1;
	spec["tau_i"]=5.0;
	return spec;
}
Fmap regular_spiking_inhibitory_ahmed_modified_par(double iext, double delay, double delay_var){
	Fmap spec = regular_spiking_excitatory_ahmed_modified_par(iext, delay, delay_var);
	spec["ei"] = -1;
	spec["tau_i"]=5.0;
	return spec;
}




Fmap intrinsic_bursting_excitatory_par(double iext){
	return intrinsic_bursting_excitatory_par(iext, 20, 15);
}
Fmap intrinsic_bursting_excitatory_ahmed_modified_par(double iext){
	return intrinsic_bursting_excitatory_ahmed_modified_par(iext, 20, 15);
}



Fmap intrinsic_bursting_excitatory_par(double iext, double delay, double delay_var){//becomes IB
	Fmap s;
	s["aA"] = 0.02;
	s["vhalfA"] = -20;
	s["aX"] = 0.2;//(2)
	s["vhalfX"] = -45;
	s["aB"] =  -0.1;
	s["vhalfB"] = -70;
	s["aW"] = 0.055;
	s["vhalfW"] = -35;
	s["am"] = 0.065;
	s["vhalfm"] = -31;
	s["lambda"] = 0.08;
	s["tauX"] = 15.0;
	s["tauB"] = 0.2;
	s["VNa"] = 55;
	s["VK"] = -72;
	s["VL"] = -50;
	s["VCa"] = 124;
	s["gNa"] = 120;//*(1+0.1*(ran2(&s_seed)-0.5));
	s["gK"] = 15;
	s["gL"] = 0.3;
	s["gA"] = 12.5;
	s["gKCa"] = 3.0;
	s["Kc"] = 2.0;
	s["Kd"] = 0.5;
	s["PCa"] = 0.15;
	s["C0"] = 2.0;
	s["Rp"] = 3.0;
	s["Kp"] = 0.75;
	s["b"] = 0.1;
	s["f"] = 0.3;
	s["B0"] = 60;
	s["D"] = 2e-9;//6e-9;
	s["taud"] = 3;
	s["tauo"] = 0.5;
	s["tau_i"] = 3;//*(1+0.1*(ran2(&s_seed)-0.5));
	//s["synapse_delay"] = 7.0*(1+0.5*(ran2(&s_seed)-0.5));
	//s["synapse_delay"] = 20.0*(1+0.5*(ran2(&s_seed)-0.5));//Aug 19th, 2013
	s["synapse_delay"] = delay + delay_var*(ran2(&s_seed)-0.5); //Aug 19th, 2013

	s["reversal_potential_e"] = -10;
	s["reversal_potential_i"] = -72;
	s["synapse_i_max"] = 0.55;  // default 0.55
	s["poisson_rate"] = 0;
	s["external_current"] = iext;//-0.08 is the threshold;
		s["nt"]=1; // anderson
	s["ei"]=1; // exc

	return s;
}

Fmap intrinsic_bursting_excitatory_ahmed_modified_par(double iext, double delay, double delay_var){//becomes IB
	Fmap spec;
	spec["a"] = 0.02; // inverse time constant of adaptation current
	spec["k"] = 0.1;
	spec["k_after"] = 2;
	spec["b"] = -0.004;

	spec["Vreset"] = -61.0;
	spec["Vtshift"] = 0.0;
	spec["Vmax"] = 30.0;
	spec["Vchangek"] = -20.0;
	spec["Vr"]=-65.0;

	spec["Rreset"] = 2;
	spec["reversal_potential_e"] = -10.0;
	spec["reversal_potential_i"] = -72.0;
	spec["synapse_i_max"] = 0.55;  // default 0.55
	spec["poisson_rate"] = 0.0;
	spec["ei"] = 1; // excitatory
	spec["d_step"]= 1.0;
	spec["external_current"] = iext;
	spec["synapse_delay"] = delay + delay_var*(ran2(&s_seed)-0.5); //Aug 19th, 2013
	spec["tau_i"]=3;
	return spec;
}


Fmap ahmed_modified_var(){
	Fmap s;
	s["membrane_potential"] = -60.02;
	s["synapse_i"] = 0;
	s["synapse_i_delayed"] = 0;
	s["synapse_e"] = 0;
	s["potassium_channel"] = 0;
	return s;
}




Fmap anderson_regular_var(){
	Fmap s;
	s["membrane_potential"] = -60.02;//*(1+0.15*(ran2(&s_seed)-0.5));
	s["synapse_i"] = 0;
	s["synapse_i_delayed"] = 0;
	s["synapse_e"] = 0;
	s["X"] = 0.00245*(1+0.1*(ran2(&s_seed)-0.5));
	s["W"] = 0.0599;
	s["B"] = 0.11967;
	return s;
}

Fmap anderson_Ca_var(){
	Fmap s = anderson_regular_var();
	s["ICa"] = 0;
	s["ICa_next"] = 0;
	s["Ca_C0"] = 0.1;
	s["Ca_C1"] = 0.1;
	s["Ca_C2"] = 0.1;
	s["Ca_C3"] = 0.1;
	s["Ca_C4"] = 0.1;
	s["Ca_C5"] = 0.1;
	s["Ca_C6"] = 0.1;
	s["Ca_U0"] = 0;
	s["Ca_U1"] = 0;
	s["Ca_U2"] = 0;
	s["Ca_U3"] = 0;
	s["Ca_U4"] = 0;
	s["Ca_U5"] = 0;
	s["Ca_U6"] = 0;
	return s;
}
