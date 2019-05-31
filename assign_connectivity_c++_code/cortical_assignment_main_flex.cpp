#include <stdio.h>
#include "cortical_assignment.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>
#include <string>
#include "srng.h"

using namespace std;
int main(int argc, char* argv[]){
	// Warning: for 64*64 needs 48G of RAM.
	// get the prelims
	time_t thistime, starttime;
	starttime=time(&thistime);
	long seed = -35;
	string conn_string;

	//double ******conn_m = new double*****[XDIM];
	//double ****d_m = new double***[XDIM];
	int syn_count[XDIM*YDIM*CPC];
	for (int i=0; i<XDIM*YDIM*CPC; i++) syn_count[i]=0;
	int conn_count[TC][TC]; // count of total number of connection between different types of cells
	int total_ee_conn_count = 0;
	int total_ei_conn_count = 0;
	int this_proc = 0;

	//double total_l23_density = calculate_l23_density();
	//double total_l5_density = calculate_l5_density();
	//double total_l4st_density = calculate_l4st_density();
	//double total_l6_density = calculate_l6_density();
	//double total_bask_density = calculate_bask_density();
	//double total_chan_density = calculate_chan_density();
	//double total_doub_density = calculate_doub_density();

	//double total_e_density = total_l23_density + total_l5_density + total_l4st_density + total_l6_density; // Density of E cells 
	//double total_i_density = total_bask_density + total_chan_density + total_doub_density; // Density of I cells
	//cout<<"total e density is "<<total_e_density<<endl;
	//cout<<"total i density is "<<total_i_density<<endl;

	// For file output
	string outfile_str = argv[2];
	string outfilesyn_str = argv[3];
	string outfilesize_str = argv[4];
	int num_proc = atoi(argv[1]);  // NUmber of MPI jobs
	ofstream outfile_stream;
	ofstream* outfilesize_stream = new ofstream[num_proc];
	ofstream* outfilesyn_stream = new ofstream[num_proc];
	char* outfilesyn;
	char* outfilesize;

	//For file that stores all the connections
	//
	//char* outfile = new char[outfile_str.size()+1];
	//strcpy(outfile, outfile_str.c_str());
	outfile_stream.open(outfile_str.c_str());

	//For files that store connections for each processor
	for (int i=0; i<num_proc; i++){
		append_file_with_proc_number(outfilesyn_str, i, outfilesyn, "_proc_");
		append_file_with_proc_number(outfilesize_str, i, outfilesize, "_proc_");
		outfilesize_stream[i].open(outfilesize);
		outfilesyn_stream[i].open(outfilesyn);
		delete [] outfilesyn;
		delete [] outfilesize;
	}


	
	// variable to store distance number
	double d_m;

	// variable to store connection strength
	double conn_m;


	
	// Initialize conn_cont
	for (int i=0; i<TC; i++){
		for (int j=0; j<TC; j++){
			conn_count[i][j] = 0;
		}
	}

	//cout<<"Assign distance matrix"<<endl;
	//assign_distance_matrix(d_m);
	//cout<<"Random check:"<<endl;
	//for (int n=0;  n<10; n++){
	//	int i = int((XDIM-1)*ran2(&seed, 2));
	//	int j = int((YDIM-1)*ran2(&seed, 2));
	//	int k = int((XDIM-1)*ran2(&seed, 2));
	//	int l = int((YDIM-1)*ran2(&seed, 2));
	//	cout<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<d_m[i][j][k][l]<<endl;
	//}


	cout<<"First pass"<<endl;

	// Do the first pass for e->e

	for (int xf=0; xf<XDIM; xf++){
		for (int yf=0; yf<YDIM; yf++){
				for (int xt=0; xt<XDIM; xt++){
					for (int yt=0; yt<YDIM; yt++){
						d_m = calculate_distance(xf, yf, xt, yt);
						for (int cnf=0; cnf<ECPC; cnf++){
							for (int cnt=0; cnt<ECPC; cnt++){
							//if (d_m>450) continue;
							connect_e_e(conn_m, d_m, xf, yf, cnf, xt, yt, cnt, conn_count[0][0]); //L23->L23
							if (conn_m<1e-50) continue;
							conn_string = replace_conn_m_with_symbol(cnf, cnt);
							print_connectivity_matrix_custom_per_row(xf, yf, cnf, xt, yt, cnt, conn_string, outfile_stream);
							this_proc = this_processor(xt, yt, cnt, num_proc);
							print_connectivity_matrix_custom_per_row_sorted(xf, yf, cnf, xt, yt, cnt, conn_string, syn_count, this_proc, num_proc, outfilesyn_stream[this_proc]);
						}	
					}
				}
			}
		}
	}

	// Compute total ee count
	for (int i=0; i<ET; i++){
		for (int j=0; j<ET; j++){
				total_ee_conn_count+=conn_count[i][j];
		}	
	}

// Compute probability of each e cell connected to an i cell
//double conn_prob_e_bask_chan_doub = prep_e_i_conn_prob(total_ee_conn_count, total_e_density, total_i_density, 300, 0.05); // Total e->i conn = 5% of total e->e conn
//cout<<"e_i probability is "<<conn_prob_e_bask_chan_doub<<endl;

cout<<"Second pass"<<endl;

// Second pass: Now do e->i connection
for (int xf=0; xf<XDIM; xf++){
	for (int yf=0; yf<YDIM; yf++){
			for (int xt=0; xt<XDIM; xt++){
				for (int yt=0; yt<YDIM; yt++){
					d_m = calculate_distance(xf, yf, xt, yt);
					for (int cnf=0; cnf<ECPC; cnf++){
						for (int cnt=ECPC; cnt<CPC; cnt++){
						//if (d_m>450) continue;
						connect_e_i(conn_m, d_m, xf, yf, cnf, xt, yt, cnt, conn_count[0][1]); //L23->BASK (Aug 14th, 2013)
						if (conn_m<1e-50) continue;
						conn_string = replace_conn_m_with_symbol(cnf, cnt);
						print_connectivity_matrix_custom_per_row(xf, yf, cnf, xt, yt, cnt, conn_string, outfile_stream);
						this_proc = this_processor(xt, yt, cnt, num_proc);
						print_connectivity_matrix_custom_per_row_sorted(xf, yf, cnf, xt, yt, cnt, conn_string, syn_count, this_proc, num_proc, outfilesyn_stream[this_proc]);

					}
				}
			}
		}
	}
}

cout<<"Last pass"<<endl;
// Now do the rest i->e, i->i
for (int xf=0; xf<XDIM; xf++){
	for (int yf=0; yf<YDIM; yf++){
			for (int xt=0; xt<XDIM; xt++){
				for (int yt=0; yt<YDIM; yt++){
					for (int cnf=ECPC; cnf<CPC; cnf++){
						for (int cnt=0; cnt<CPC; cnt++){
						d_m = calculate_distance(xf, yf, xt, yt);
						//if (d_m>450) continue;
						if (cnt<ECPC) connect_i_e(conn_m, d_m, xf, yf, cnf, xt, yt, cnt, conn_count[1][0]); //BASK->L23
						else connect_i_i(conn_m, d_m, xf, yf, cnf, xt, yt, cnt, conn_count[1][1]); //BASK->L23
						if (conn_m<1e-50) continue;
						conn_string = replace_conn_m_with_symbol(cnf, cnt);

						print_connectivity_matrix_custom_per_row(xf, yf, cnf, xt, yt, cnt, conn_string, outfile_stream);
						this_proc = this_processor(xt, yt, cnt, num_proc);
						print_connectivity_matrix_custom_per_row_sorted(xf, yf, cnf, xt, yt, cnt, conn_string, syn_count, this_proc, num_proc, outfilesyn_stream[this_proc]);

					}
				}
			}
		}
	}
}

//Compute total e->i conn count

for (int i=0; i<ET; i++){
	for (int j=ET; j<TC; j++){
		total_ei_conn_count+=conn_count[i][j];
	}
}

cout<<"List conn_count entries"<<endl;
for (int i=0; i<TC; i++){
	for (int j=0; j<TC; j++){
		cout<<i<<" "<<j<<" "<<conn_count[i][j]<<endl;
	}
}

for (int j=0; j<TC; j++){
	cout<<"Total e and i synapses going to cell type "<<j<<" (respectively) ";
	int total_e=0;
	int total_i=0;
	for (int i=0; i<TC; i++){
		if (i<ET) total_e+=conn_count[i][j];
		else total_i+=conn_count[i][j];
	}
	cout<<total_e<<"\t"<<total_i<<endl;
}
cout<<endl;
		


cout<<"Total e->e count is "<<total_ee_conn_count<<endl;
cout<<"Total e->i count is "<<total_ei_conn_count<<endl;

cout<<"Writing conectivity matrix size to files"<<endl;
for (int i=0; i<(XDIM*YDIM*CPC); i++){
		int this_proc = this_processor(i, num_proc);
		outfilesize_stream[this_proc] << syn_count[i]<<endl;
}

//for (int i=0; i<num_proc; i++){
//	sprintf(rank_char, "%d", i);
//	std::string rank_string(rank_char);
//	rank_string = std::string("_proc_") + rank_string;
//	outfile_str = outfile_refresh_str;
//	outfile_str.insert(outfile_str.length()-4, rank_string);
//	char* outfile = new char[outfile_str.size()+1];
//	strcpy(outfile, outfile_str.c_str());
//	outfile_stream.open(outfile, ios::binary);
//	cout<<"File is: "<<outfile<<endl;
	//print_num_neurons(outfile_stream);
	//print_seed(outfile_stream, i);
//	print_connectivity_matrix(conn_m, outfile_stream, i, num_proc);
//	delete [] outfile;
//	outfile_stream.close();
//}

outfile_stream.close();
for (int i=0; i<num_proc; i++){
		outfilesyn_stream[i].close();
		outfilesize_stream[i].close();
	}



cout<<"Programme ends"<<endl;
cout<<"Time taken "<<time(&thistime)-starttime<<" seconds."<<endl;
return 0;
}





