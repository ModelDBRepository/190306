#include "Neuron.h"
#include <sstream>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <unordered_map>

// Utility functions for manipulating the I/O of Neuron class

using namespace std;
typedef unordered_map<string, double> Fmap;

void read_neuron_data(istream& ins, Fmap* this_map){
  string s;
  while (getline(ins, s) && !s.empty()){
    stringstream ss(s);
    string key;
    double value;
    ss>>key>>value;
    (*this_map)[ key ]=value;
  }
  assert(s.empty());
  return;
}

void write_neuron_data(ostream& outs, Neuron& neu){
  Fmap var_hash;
  Fmap par_hash;
  par_hash = neu.cellParameters::get_hash();
  var_hash = neu.cellVariables::get_hash();
  for(Fmap::const_iterator it = par_hash.begin(); it != par_hash.end(); ++it)
    {
      outs<<it->first<<"\t"<<it->second<<endl;
    }
  outs<<"\n";
  for(Fmap::const_iterator it = var_hash.begin(); it != var_hash.end(); ++it)
    {
      outs<<it->first<<"\t"<<it->second<<endl;
    }
  outs<<"\n";
  return;
}

void read_synapse(istream &ins, int* N, double** syn){
  double syn_val;
  string s;
  for (int i=0; i<(*N); i++){
    //cout<<"read synapse i is "<<i<<endl;
    getline(ins, s);
		cout<<s<<endl;
    assert(!s.empty());
    stringstream ss(s);
    for (int j=0; j<(*N); j++){
      ss>>syn_val;
      //cout<<"syn_val "<<syn_val<<endl;
      syn[i][j]=syn_val;      
    }
  }
  getline(ins,s);
  assert(s.empty());
  return;
}

void write_synapse(ostream &outs, int* N, double** syn){
  assert((*N)>0);
  for(int i=0; i<(*N); i++){
    for (int j=0; j<(*N); j++){
      outs<<syn[i][j];
      if (j==((*N)-1)){
	outs<<"\n";
      }
      else{
	outs<<"\t";
      }
    }
  }
  outs<<"\n";
  return;
}

  
void read_other_info(istream& ins, int* N, long* seed, double* s_time, int* iset, float* gset, long* idum2, long* iy, long* iv, int iv_size){
  string s;
  int val1;
  long val2;
  double val3;
  long val4;
  long val5;
  int val6;
  float val7;
  std::getline(ins, s);
  //std::getline(ins, s);
  //cout<<"s is "<<&s<<endl;
  //cout<<"ins address "<<&ins<<endl;
  stringstream ss(s);
  //cout<<"ss is "<<ss<<endl;
  ss>>val1>>val2>>val3>>val6>>val7>>val4>>val5;
  //cout<<"val1 is "<<val1<<endl;
  //cout<<"val2 is "<<val2<<endl;
  //cout<<"val3 is "<<val3<<endl;
  (*N) = val1;
  (*seed) = val2;
  (*s_time) = val3;
  (*iset) = val6;
  (*gset) = val7;
  (*idum2) = val4;
  (*iy) = val5;
  for (int i=0; i<iv_size; i++){
    ss>>(*(iv+i));
  }
  getline(ins,s);
  assert(s.empty());
  return;
}


void read_other_info(istream& ins, int* N, long* seed, double* s_time){
	string s;
	int val1;
	long val2;
	double val3;
	std::getline(ins, s);
	stringstream ss(s);
	ss>>val1>>val2>>val3;
	(*N) = val1;
	(*seed) = val2;
	(*s_time) = val3;
	getline(ins,s);
	assert(s.empty());
	return;
}

void write_other_info(ostream &outs, int *N, long *seed, double* s_time,  int* iset, float* gset, long* idum2, long* iy, long* iv, int iv_size){
  outs<<*N<<"\t"<<*seed<<"\t"<<*s_time<<"\t"<<*iset<<"\t"<<*gset<<"\t"<<*idum2<<"\t"<<*iy<<"\t";
  for (int i=0; i<iv_size; i++){
    outs<<(*(iv+i));
    if (i<(iv_size-1)){outs<<"\t";}
  }
  outs<<endl;
  outs<<endl;
  return;
}
  
void parse_checkpoint_other_and_synapse_info(ifstream& ins, int* N, long* seed, double* s_time,  int* iset, float* gset, long* idum2, long* iy, long* iv, int iv_size, double** syn){
  //read N, seed and start time
  read_other_info(ins, N, seed, s_time, iset, gset, idum2, iy, iv, iv_size);
  //read synaptic values
  read_synapse(ins, N, syn);
  //read start cell par and var
  //create pars and vars vectors
  //assert(pars_vec==NULL);
  //assert(vars_vec==NULL);
  //cout<<"parse *N is "<<*N<<" "<<N<<endl;
  assert((*N)>0);
  //cout<<"Before vars_vec is "<<vars_vec<<endl;
  //vars_vec = new Fmap[(*N)];
  //pars_vec = new Fmap[(*N)];
  //cout<<"Inside function vars_vec is "<<vars_vec<<endl;
  //read neuron data
  //for (int i=0; i<(*N); i++){
  //read_neuron_data(ins, (pars_vec+i));
    //string s;
    //getline(ins, s);
    //assert(s.empty());
    //read_neuron_data(ins, (vars_vec+i));
  //}
  return;
}

void write_checkpoint_data(ostream &outs, int* N, long* seed, double* s_time,  int* iset, float* gset, long* idum2, long* iy, long* iv, int iv_size, double** syn, Neuron* neu_array){
  //write N, seed and start time
  //cout<<"Are we here"<<endl;
  write_other_info(outs, N, seed, s_time, iset, gset, idum2, iy, iv, iv_size);
  //write synaptic values
  write_synapse(outs, N, syn);
  //write neuron data
  for (int i=0; i<(*N); i++){
    write_neuron_data(outs, *(neu_array+i));
  }
  return;
}
