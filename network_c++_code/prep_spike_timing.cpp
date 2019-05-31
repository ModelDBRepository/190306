#include <stdio.h>
#include "cortical_assignment.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>
#include <string>
#include <assert.h>
#include <sstream>

int main(int argc, char* argv[]){
	int total_num_proc = atoi(argv[1]);
	string result_file_to_read = argv[2];
	string out_file = argv[3];

	// File handlers
	ifstream* in = new ifstream[total_num_proc];
	ofstream out;

	// Variable for location of neurons
	int x, y, cn, neuron_number;
	int layer, cell_type, cell_number;
	double time;

	char** file_to_read = new char*[total_num_proc];
	char* file_to_write;
	string s;

	// Create the name for the output file
	file_to_write = new char[out_file.size()+1];
	strcpy(file_to_write, out_file.c_str());
	out.open(file_to_write);

	cout<<"Reading result files and processing."<<endl;
	for (int i=0; i<total_num_proc; i++){
		append_file_with_proc_number(result_file_to_read, i, file_to_read[i], "_proc_");
			
		in[i].open(file_to_read[i]);
		while (getline(in[i], s) && !s.empty()){
			stringstream ss(s);
			ss>>neuron_number;
			ss>>time;
			reverse_indexing_neurons(neuron_number, x, y, cn);
			custom_convert_cn_to_lcc(layer, cell_type, cell_number, cn);
			out<<time<<"\t"<<neuron_number<<"\t"<<x<<"\t"<<y<<"\t"<<layer<<"\t"<<cell_type<<"\t"<<cell_number<<endl;
		}
		in[i].close();
	}
	out.close();
	cout<<"End of programme"<<endl;
	return 0;
}
	




					


