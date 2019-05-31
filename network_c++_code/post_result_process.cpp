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
	int total_neurons = atoi(argv[2]);
	string result_file_to_read = argv[3];
	string intermediate_file = argv[4];
	string grand_file = argv[5];
	string s;
	string* s_part = new string[total_num_proc];

	int quotient = total_neurons/total_num_proc;
	assert(total_neurons%total_num_proc==0);
	ifstream* in = new ifstream[total_num_proc];
	ifstream* inter = new ifstream[total_num_proc];
	ofstream out;
	ofstream out2;
	ofstream out3;
	double time = 0;
	int count;

	char** file_to_open_read = new char*[total_num_proc];
	char** file_to_open_write = new char*[total_num_proc];
	char* grand_L23_file;
	char* grand_L5_file;
	char* grand_L6_file;

	double *s_value = new double[total_neurons];
	double *s_value_prev = new double[total_neurons];
	int *spiked = new int[quotient];
	int g_spiked;
	cout<<"Reading result files and processing."<<endl;
	for (int i=0; i<total_num_proc; i++){
		append_file_with_proc_number(result_file_to_read, i, file_to_open_read[i], "_proc_");
		append_file_with_proc_number(intermediate_file, i, file_to_open_write[i], "_proc_");
		in[i].open(file_to_open_read[i]);
		out.open(file_to_open_write[i]);
		cout<<"File is "<<file_to_open_read[i]<<endl;
		count = 0;
		for (int m=0; m<quotient; m++){
		s_value_prev[m] = 10;
		}
		time = 0;
		while((getline(in[i],s)) && !s.empty() && time<=50000){
			if (s != ""){
				size_t found = s.find("#");
				if (found>=string::npos){
					count++;
					stringstream ss(s);
					ss>>time;
					//cout<<"Time is "<<time<<endl;
					if (count%10!=0) continue; // Sample every 10 ms
					for (int j=0; j<quotient; j++){
						//if (j==0){cout<<"time is "<<time<<endl;}
						ss>>s_value[j];
						double diff = s_value[j] - s_value_prev[j];
						if (diff > 1e-15){
							//there is a spike
							spiked[j] = 1;
						}
						else{
							spiked[j] = 0;
						}
						out<<spiked[j]<<"\t";
						s_value_prev[j] = s_value[j];
					}
					out<<"\n";
					//count++;
					//cout<<"Count is "<<count<<endl;
				}
			}
		}
		out.close();
		in[i].close();
	}

	for (int i = 0; i<total_num_proc; i++){
		//char* file_to_open_write;
		//append_file_with_proc_number(intermediate_file, i, file_to_open_write);
		inter[i].open(file_to_open_write[i]);
	}
	count = -1;
	cout<<"Writing movie files."<<endl;
	do{
		count++;
		s="";
		for (int i= 0; i<total_num_proc; i++){
			getline(inter[i], s_part[i]);
			s+=s_part[i];
		}
		
		append_file_with_proc_number(grand_file, count, grand_L23_file, "_f_");
	//	append_file_with_proc_number(grand_file, count, grand_L5_file, "_l5_f_");
//append_file_with_proc_number(grand_file, count, grand_L6_file, "_l6_f_");
		out.open(grand_L23_file);
		//out2.open(grand_L5_file);
//out3.open(grand_L6_file);

		if (s != "" && !s_part[0].empty()){
			stringstream ss(s);
					for (int j=0; j<total_neurons; j++){
						int x, y, cn;
						reverse_indexing_neurons(j, x, y, cn);
						ss>>g_spiked;
						switch(cn){
						case(0):
							out<<2*x<<"\t"<<2*y<<"\t"<<g_spiked<<"\n";
							break;
						case(1):
							out<<2*x+1<<"\t"<<2*y<<"\t"<<g_spiked<<"\n";
							break;
						case(2):
							out<<2*x<<"\t"<<2*y+1<<"\t"<<g_spiked<<"\n";
						case(3):
							out<<2*x+1<<"\t"<<2*y+1<<"\t"<<g_spiked<<"\n";
							break;
						case(5):
							//out2<<2*x<<"\t"<<2*y<<"\t"<<g_spiked<<"\n";
							break;
						case(6):
							//out2<<2*x+1<<"\t"<<2*y<<"\t"<<g_spiked<<"\n";
							break;
						case(7):
							//out2<<2*x<<"\t"<<2*y+1<<"\t"<<g_spiked<<"\n";
						case(8):
							//out2<<2*x+1<<"\t"<<2*y+1<<"\t"<<g_spiked<<"\n";
							case(9):
							//out3<<2*x<<"\t"<<2*y<<"\t"<<g_spiked<<"\n";
							break;
						case(10):
							//out3<<2*x+1<<"\t"<<2*y<<"\t"<<g_spiked<<"\n";
							break;
						case(11):
							//out3<<2*x<<"\t"<<2*y+1<<"\t"<<g_spiked<<"\n";
						case(12):
							//out3<<2*x+1<<"\t"<<2*y+1<<"\t"<<g_spiked<<"\n";

							break;

						default:
							break;
						}
					}
		}
	out.close();
	//out2.close();
	//out3.close();
	}while (!(s_part[0]).empty());
			

		for (int i= 0; i<total_num_proc; i++){
			inter[i].close();
	
		}

		//out.close();
	
	delete [] in;
	delete [] inter;
	delete [] s_part;
	cout<<"Programme ends."<<endl;
		return 0;
}




					


