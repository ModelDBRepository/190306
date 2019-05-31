#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <unordered_map>
#include <string>
#define XDIM 9 //default 28
#define YDIM 9
#define WIDTH 25
#define CPC 16
#define ECPC 12 // number of exciatory cells per column (original 12; changed to 16 for ex only--Dec 3, 2014)
#define TC 2 // Types of cells
#define ET 1 //Types of e cells
#define PI 3.1415927
#define DIV_E 1.0
//#define CRITDEX (125.00/WIDTH)//default 300.00/WIDTH but when XDIM is <= 16, we use 75
#define CRITDEX (100.00/WIDTH)
using namespace std;
//template <typename TYPE>;
typedef unordered_map<string, double> Fmap;
int this_processor(int neuron_index, int num_proc);
int this_processor(int x, int y, int cn, int num_proc);
int local_neuron_index(int neuron_index, int this_proc, int num_proc);
double calculate_distance(int x1, int y1, int x2, int y2);
string replace_conn_m_with_symbol(int cnf, int cnt);
void assign_distance_matrix(double**** distance_matrix);
int check_distance_within(int xf, int yf, int xt, int yt, double**** distance_matrix, double dist);
void connect_e_e(double& conn, double distance, int xf, int yf, int cnf, int xt, int yt, int cnt, int& count);
void connect_e_i(double& conn, double distance, int xf, int yf, int cnf, int xt, int yt, int cnt, int& count);
void connect_i_e(double& conn, double distance, int xf, int yf, int cnf, int xt, int yt, int cnt, int& count);
void connect_i_i(double& conn, double distance, int xf, int yf, int cnf, int xt, int yt, int cnt, int& count);

void connect_doub_chan(double& conn, double**** distance_matrix, int xf, int yf, int cnf, int xt, int yt, int cnt, int& count);

double prep_a_b_dying_probability(double b_density, double range, double target_connection_to_b_for_one_a);
double prep_e_i_conn_prob(double total_e_e_conn, double e_density, double i_density, double range, double proportion);
double assign_connect_internal_by_probability(int xf, int yf, int xt, int yt, double distance, double max_dist, double prob, int& count);
double assign_connect_by_probability(int xf, int yf, int xt, int yt, double distance , double max_dist, double prob, int& count);
double assign_connect_by_probability_no_prune(int xf, int yf, int xt, int yt, double distance , double max_dist, double prob, int& count);
double assign_connect_by_dying_probability(int xf, int yf, int xt, int yt, double distance, double range, double prob, int& count);
double impose_boundary_conditions(int x, int y);
int calculate_total_number_of_cells();
double calculate_density();
void forward_indexing_neurons(int& indexed_number, int x, int y, int cn);
void reverse_indexing_neurons(int indexed_number, int& x, int& y, int& cn);
void custom_convert_lcc_to_cn(int layer, int cell_type, int cell_number, int& cn);
void custom_convert_cn_to_lcc(int& layer, int& cell_type, int& cell_number, int cn);
void print_num_neurons(ostream& out);
void read_num_neurons(istream& in, int* num_neurons);
void print_seed(ostream& out, int this_proc);
void print_connectivity_matrix_custom(double****** conn_m, ostream& out);
template <class	TYPE> void print_connectivity_matrix_custom_per_row(int i, int xf, int yf, int cnf, int j, int xt, int yt, int cnt, TYPE conn, ostream&out);
template <class	TYPE> void print_connectivity_matrix_custom_per_row(int xf, int yf, int cnf, int xt, int yt, int cnt, TYPE conn, ostream&out);
template <class	TYPE> void print_connectivity_matrix_custom_per_row_sorted(int xf, int yf, int cnf, int xt, int yt, int cnt, TYPE conn, int syn_count[], int this_proc, int num_proc, ostream&out);
void print_connectivity_matrix_testing(double** conn_m, ostream& out, int this_proc, int num_proc);
void print_connectivity_matrix(double****** conn_m, ostream& out, int this_proc, int num_proc);
void print_connectivity_matrix(ostream& out, int this_proc, int num_proc); //for custom conn matrix setup
void read_seed(istream& in, long* seed);
void read_syn_size(int *syn_size, int size_of_array, istream& in_size);
void read_connectivity_matrix(double **syn, istream& in, int num_neurons_this_proc, int m_size, int num_proc, int this_proc);
//void read_connectivity_matrix_custom(double **syn, istream& in, int num_neurons_this_proc, int m_size, int num_proc, int this_proc);
void read_connectivity_matrix_custom(double **syn_sorted, int **syn_location, int *syn_counter, int* syn_size, istream& in_conn);
void read_connectivity_matrix_custom_per_row(int& neuron_numf, int& neuron_numt, double& conn_value, istream& in);
void assign_neurons_testing(Fmap* var_vec, Fmap* par_vec, int num_neurons_this_proc, int total_neurons, int num_proc, int this_proc);
void assign_neurons(Fmap* var_vec, Fmap* par_vec, int num_neurons_this_proc, int total_neurons, int num_proc, int this_proc);
void assign_neurons_testing(Fmap* var_vec, Fmap* par_vec, int num_neurons_this_proc, int total_neurons, int num_proc, int this_proc, double, double, double);
void assign_neurons(Fmap* var_vec, Fmap* par_vec, int num_neurons_this_proc, int total_neurons, int num_proc, int this_proc, double, double, double);
void assign_neurons_ahmed_modified(Fmap* var_vec, Fmap* par_vec, int num_neurons_this_proc, int total_neurons, int num_proc, int this_proc, double rs_iext, double ib_iext, double fs_iext);


void assign_neurons(Fmap* var_vec, Fmap* par_vec, int num_neurons_this_proc, int total_neurons, int num_proc, int this_proc, double, double, double, double, double);
void assign_neurons_ahmed_modified(Fmap* var_vec, Fmap* par_vec, int num_neurons_this_proc, int total_neurons, int num_proc, int this_proc, double rs_iext, double ib_iext, double fs_iext, double delay, double delay_var);

void append_file_with_proc_number(string file_name_str, int proc_number, char*& result_file_name, string name_type);

Fmap if_excitatory_par(double, double, double);
Fmap if_inhibitory_par(double, double, double);

Fmap intrinsic_bursting_excitatory_par(double);
Fmap intrinsic_bursting_excitatory_ahmed_modified_par(double);

Fmap intrinsic_bursting_excitatory_par(double, double, double);
Fmap intrinsic_bursting_excitatory_ahmed_modified_par(double, double, double);

Fmap regular_spiking_inhibitory_par(double);
Fmap regular_spiking_inhibitory_ahmed_modified_par(double);

Fmap regular_spiking_inhibitory_par(double, double, double);
Fmap regular_spiking_inhibitory_ahmed_modified_par(double, double, double);

Fmap regular_spiking_excitatory_par(double);
Fmap regular_spiking_excitatory_ahmed_modified_par(double);

Fmap regular_spiking_excitatory_par(double, double, double);
Fmap regular_spiking_excitatory_ahmed_modified_par(double, double, double);

Fmap fast_spiking_inhibitory_par(double);
Fmap fast_spiking_inhibitory_ahmed_modified_par(double);

Fmap fast_spiking_inhibitory_par(double, double, double);
Fmap fast_spiking_inhibitory_ahmed_modified_par(double, double, double);

Fmap ahmed_modified_var();
Fmap anderson_regular_var();
Fmap anderson_Ca_var();
Fmap if_var();

