#include <stdio.h>
#include "cortical_assignment.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>
#include <string>
#include <assert.h>
#include "srng.h"

using namespace std;
int main(int argc, char* argv[]){
	// Warning: for 64*64 needs 48G of RAM.
	// get the prelims
	//
	//
	//
	//
	// CPU time calculation
	time_t thistime, starttime;
	starttime=time(&thistime);

	int num_proc = atoi(argv[1]);  // NUmber of MPI jobs
	ifstream infile_stream; // stream for input file
	string outfile_str = argv[2];
	string outfile_refresh_str = outfile_str;
	ofstream outfile_stream;
	char rank_char[10];

		
	cout<<"Writing total number of neurons and seed info to files"<<endl;
	for (int i=0; i<num_proc; i++){
		sprintf(rank_char, "%d", i);
		std::string rank_string(rank_char);
		rank_string = std::string("_proc_") + rank_string;
		outfile_str = outfile_refresh_str;
		outfile_str.insert(outfile_str.length()-4, rank_string);
		char* outfile = new char[outfile_str.size()+1];
		strcpy(outfile, outfile_str.c_str());
		outfile_stream.open(outfile, ios::binary);
		cout<<"File is: "<<outfile<<endl;
		//print_num_neurons(outfile_stream);
		//print_seed(outfile_stream, i);
		print_connectivity_matrix(outfile_stream, i, num_proc);
		delete [] outfile;
		outfile_stream.close();
	}



cout<<"Programme ends"<<endl;
cout<<"Time taken "<<time(&thistime)-starttime<<" seconds."<<endl;
return 0;
}





