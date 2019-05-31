//#include "Neuron.h"
#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>

//#include "cellVariables.h"
//#include "cellParameters.h"
using namespace std;
typedef unordered_map<string, double> Fmap;

void read_neuron_data(istream&, Fmap*);
void write_neuron_data(ostream&, Neuron&);
void read_synapse(istream&, int*, double**);
void write_synapse(ostream&, int*, double**);
void read_other_info(istream& ins, int* N, long* seed, double* s_time, int* iset, float* gset, long* idum2, long* iv, long* iy, int iy_size);
void read_other_info(istream& ins, int* N, long* seed, double* s_time);
void write_other_info(ostream &outs, int *N, long *seed, double* s_time,  int* iset, float* gset, long* idum2, long* iv, long* iy, int iy_size);
void parse_checkpoint_other_and_synapse_info(ifstream& ins, int* N, long* seed, double* s_time,  int* iset, float* gset, long* idum2, long* iv, long* iy, int iy_size, double** syn);
void write_checkpoint_data(ostream &outs, int* N, long* seed, double* s_time,  int* iset, float* gset, long* idum2, long* iv, long* iy, int iy_size, double** syn, Neuron* neu_array);
