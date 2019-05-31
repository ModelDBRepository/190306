#include "Neuron.h"
#include "Neuron_read_write.h"
#include <omp.h>
#include <unordered_map>
#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "srng.h"
#include "gasdev.h"
#include <mpi.h>
#include "cortical_assignment.h"
#include "LFP.h"
#include "LFP_tools.h"
#include "Diffusion.h"
#include "Diffusion_tools.h"

// End change
using namespace std;

int main(int argc, char* argv[]){
  
  // MPI starts here
  int required = MPI_THREAD_SERIALIZED;
  int provided;
  //int ierr = MPI_Init(&argc, &argv);
  int ierr=MPI_Init_thread(&argc, &argv, required, &provided);
  int rank, n_workers;
  
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &n_workers);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    

  
  // Time counting
  time_t thistime, starttime, compustarttime;
  starttime=time(&thistime);
  //rout<<starttime<<endl;
  //
  
  // Variable declarations
  typedef unordered_map<string, double> Fmap;
  Neuron** neurons_row;  // array to store neuron information (1 of the 2 dimensions is for openmp padding)
	Neuron** neurons_row_transposed;
  double **syn_sorted; // synaptic connectivitiy array
	int **syn_location;
	int *syn_counter;
	int *syn_size;
  //int **valid_entries; // mapping of synaptic connectivity to valid syn entries
  //int *num_valid_entries; // mapping of valid_entries to syn
  double negative_eps=-1e-30;
	double *esyn_array; // array storing Esyn of the neurons in each worker
	double *esyn_grand_array; // array storing Esyn for ALL the neurons
	int num_neurons = 0; // Total number of all neurons
	double *array_to_receive; // array to store si for ALL the neurons
	double *array_to_send; // array to store si for neurons in each worker (padded for openmp)
	double *K_extracellular_array_to_receive; // array to store K_extracellular for ALL the neurons (DIFFUSION)
	double *K_extracellular_array_to_send; // array to store K_extracellular for each worker (DIFFUSION)
	
	double *potentials_to_send; // array to store membrane potential readings for each worker
	double *potentials_to_receive; // array to store membrane potential readings for all the neurons (after MPI gather)

	double *Ca_intracellular_array_to_send;
	double *Ca_intracellular_array_to_receive;

	double *KB_array_to_send;
	double *KB_array_to_receive;

	Diffusion *diffusion_array; // for DIFFUSION of K 


	//float *Ca_array_to_send; // array to store ICa readings for each worker
	//float *Ca_array_to_receive; // array to store ICa readings for all the neurons (after MPI gather)
	int use_vmd_flag; // flag to decide whether to use VmD or not (1 = use; 0 = use poisson)



	double **array_to_send_transposed;
	int thid; // Thread ID for openmp
	int network_spiked; // flag to tell whether the network has spiked

	// Integration time step
  double eps=1e-15; // to compensate the imprecise representation of double by C++ (e.g. 0.1 is different from 0.1000).  Last digits might affect the condition on for loop.
  //double dt=0.005-eps;
  double dt=0.01-eps;
	// Random seed initialization
  long **seed;  // For Poisson spike train  
   
  //Files definition:
	ifstream conn_file_stream; //ifstream handler for sorted connectivity file 
	ifstream connsize_file_stream; //ifstream handler for connecitivity matrix size file
  ifstream infile_cp; //ifstream handler for config file (seed and total number of neurons)
  ofstream outfile_cp; // ofstream handler (for outputting to checkppint file)
  ofstream rout; // handler for writing to data file
	ofstream sout; // handler for writing spike file
	ofstream lfp23out; //handler for writing LFP L23 file
	ofstream lfp4out; //handler for writing LFP L4 file
	ofstream lfp5out; //handler for writing LFP L5 file
	ofstream lfp6out; //handler for writing LFP L6 file
	ofstream lfpinout; //handler for writing LFP interneurons file
	ofstream K23out; //handler for writing K L23 file
	ofstream K4out; //handler for writing K L4 file
	ofstream K5out; //handler for writing K L5 file
	ofstream Kinout; //handler for writing K interneuron file
	ofstream Ca23out;
	ofstream Cainout;
	ofstream KB23out;
	ofstream KBinout;

  char rank_char[10]; //char for MPI rank
  sprintf(rank_char, "%d", rank);
	std::string rank_string(rank_char);
  rank_string = std::string("_proc_") + rank_string;

  
  cout<<"CP 1"<<endl;
  string conn_file_to_read_str = argv[1];  //sorted connectivity file (by processor)
	conn_file_to_read_str.insert(conn_file_to_read_str.length()-4, rank_string);
  conn_file_stream.open(conn_file_to_read_str.c_str());

	cout<<"CP 2"<<endl;
  string connsize_file_to_read_str = argv[2];  //sorted connectivity file (by processor)
	connsize_file_to_read_str.insert(connsize_file_to_read_str.length()-4, rank_string);
	connsize_file_stream.open(connsize_file_to_read_str.c_str());


	cout<<"CP 3"<<endl;
  string cp_file_to_read_str = argv[3];  //config file (neuron number and seed)
	cp_file_to_read_str.insert(cp_file_to_read_str.length()-4, rank_string);
	infile_cp.open(cp_file_to_read_str.c_str(), ios::binary);
	

  cout<<"CP 4"<<endl;
  string cp_file_to_write_str = argv[4];  //file for checkpoint to write
	cp_file_to_write_str.insert(cp_file_to_write_str.length()-4, rank_string);
	outfile_cp.open(cp_file_to_write_str.c_str(), ios::app);
  outfile_cp.precision(22);

  cout<<"CP 5"<<endl;
  string results_file_str = argv[5];  //file for results to write
	results_file_str.insert(results_file_str.length()-4, rank_string);
	rout.open(results_file_str.c_str()); //no append
	rout.precision(7);


	cout<<"CP 6"<<endl;
  string spike_file_str = argv[6];  //file for spike results to write
	spike_file_str.insert(spike_file_str.length()-4, rank_string);
  sout.open(spike_file_str.c_str()); //no append
  sout.precision(7);

	cout<<"CP 7"<<endl;
	if (rank==0){
		string lfp23_results_file_str = argv[7];  //file for LFP 2/3 results to write
		lfp23out.open(lfp23_results_file_str.c_str()); //no append
		lfp23out.precision(7);

		string lfp4_results_file_str = argv[8];  //file for LFP 4 results to write
		lfp4out.open(lfp4_results_file_str.c_str()); //no append
		lfp4out.precision(7);


		string lfp5_results_file_str = argv[9];  //file for LFP 5 results to write
		lfp5out.open(lfp5_results_file_str.c_str()); //no append
		lfp5out.precision(7);


		string lfp6_results_file_str = argv[10];  //file for LFP 6 results to write
		lfp6out.open(lfp6_results_file_str.c_str()); //no append
		lfp6out.precision(7);


		string lfpin_results_file_str = argv[11];  //file for LFP in results to write
		lfpinout.open(lfpin_results_file_str.c_str()); //no append
		lfpinout.precision(7);

		string K23_results_file_str = argv[12];  //file for K 2/3 results to write
		K23out.open(K23_results_file_str.c_str()); //no append
		K23out.precision(7);

		string K4_results_file_str = argv[13];  //file for K 4 results to write
		K4out.open(K4_results_file_str.c_str()); //no append
		K4out.precision(7);

		string K5_results_file_str = argv[14];  //file for K 5 results to write
		K5out.open(K5_results_file_str.c_str()); //no append
		K5out.precision(7);

		string Kin_results_file_str = argv[15];  //file for K in results to write
		Kinout.open(Kin_results_file_str.c_str()); //no append
		Kinout.precision(7);

		string Ca23_results_file_str = argv[16];  //file for Ca intracellular L23 results to write
		Ca23out.open(Ca23_results_file_str.c_str()); //no append
		Ca23out.precision(7);

		string Cain_results_file_str = argv[17];  //file for Ca intracellular in results to write
		Cainout.open(Cain_results_file_str.c_str()); //no append
		Cainout.precision(7);

		string KB23_results_file_str = argv[18];  //file for Ca intracellular L23 results to write
		KB23out.open(KB23_results_file_str.c_str()); //no append
		KB23out.precision(7);

		string KBin_results_file_str = argv[19];  //file for Ca intracellular in results to write
		KBinout.open(KBin_results_file_str.c_str()); //no append
		KBinout.precision(7);





	}
cout<<"Getting past all CPs"<<endl;
  int conti = 0;
  
  Fmap cell_vars, cell_pars;
  Fmap *var_vec;
  Fmap *par_vec;
	cout<<"After var par declaration"<<endl;

	//Iext parameters
	double rs_iext = atof(argv[20]);
	double ib_iext = atof(argv[21]);
	double fs_iext = atof(argv[22]);
	cout<<"After iext assignment"<<endl;
	//Use vmd flag
	use_vmd_flag = atoi(argv[23]);
	assert(use_vmd_flag==0 || use_vmd_flag==1);
	cout<<"After use vmd flag"<<endl;

	//Poisson rate
	double poss_rate = atof(argv[24]); // in events per ms
	double poss_increment = atof(argv[25]); //default 0.2
	cout<<"Poss rate is "<<poss_rate<<endl;
	cout<<"Poss increment is "<<poss_increment<<endl;

	double poss_rate_inhib = atof(argv[26]);
	double poss_increment_inhib = atof(argv[27]);


	//VmD parameters
	double ge_zero = atof(argv[28]); // to e cells
	double sigma_e = atof(argv[29]);
	double gi_zero = atof(argv[30]);
	double sigma_i = atof(argv[31]);
	double ge_zero_i = atof(argv[32]); // to i cells
	double sigma_e_i = atof(argv[33]);
	double gi_zero_i = atof(argv[34]);
	double sigma_i_i = atof(argv[35]);

	double maxtime = atof(argv[37]);
	cout<<"After poisson and vmd parameters"<<endl;

	//Delay parameters
	double syn_delay = atof(argv[38]);
	double syn_delay_var = atof(argv[39]);
	assert(syn_delay>(syn_delay_var/2.0));

	//Glial parameters
	double k_for = atof(argv[40]);
	double k_back = atof(argv[41]);
	double rise = atof(argv[42]);
	double Ko_eq_pump = atof(argv[43]);
	double Ko_eq_glia = atof(argv[44]);
	double max_pump_current = atof(argv[45]);
	cout<<"AFter first glial"<<endl;

	double k_for_in = atof(argv[46]);
	double k_back_in = atof(argv[47]);
	double rise_in = atof(argv[48]);
	double Ko_eq_pump_in = atof(argv[49]);
	double Ko_eq_glia_in = atof(argv[50]);
	double max_pump_current_in = atof(argv[51]);

	double k_for_in_2 = atof(argv[52]);
	double k_back_in_2 = atof(argv[53]);
	double rise_in_2 = atof(argv[54]);
	double Ko_eq_pump_in_2 = atof(argv[55]);
	double Ko_eq_glia_in_2 = atof(argv[56]);
	double max_pump_current_in_2 = atof(argv[57]);

	cout<<"After second glial"<<endl;
	cout<<"Max pum current in is "<<max_pump_current_in<<endl;


	//double gglia_in = atof(argv[38]);
	//double kzeroinf_in = atof(argv[39]);

	// Stimulation parameters
	double stim_start = atof(argv[58]);
	double stim_end = atof(argv[59]);
	double stim_strength = atof(argv[60]);

	cout<<"Stim end is"<<stim_end<<endl;

	double stim_start_in = atof(argv[61]);
	double stim_end_in = atof(argv[62]);
		double stim_strength_in = atof(argv[63]);
  // Inhibitory ramp factor
		double inhib_ramp_factor = atof(argv[64]);

  // LFP setup
	cout<<"Before LFP setup, rank = "<<rank<<endl;
	int LFP_n_per_dim = atoi(argv[36]);
	int LFP_n_per_dim_sq = LFP_n_per_dim*LFP_n_per_dim;
	LFP* LFParray = new LFP[LFP_n_per_dim_sq];
	float* LFP_res_L23 = new float[LFP_n_per_dim_sq];
	float* LFP_recv_L23 = new float[LFP_n_per_dim_sq];
	float* LFP_res_L4 = new float[LFP_n_per_dim_sq];
	float* LFP_recv_L4 = new float[LFP_n_per_dim_sq];
	float* LFP_res_L5 = new float[LFP_n_per_dim_sq];
	float* LFP_recv_L5 = new float[LFP_n_per_dim_sq];
	float* LFP_res_L6 = new float[LFP_n_per_dim_sq];
	float* LFP_recv_L6 = new float[LFP_n_per_dim_sq];
	float* LFP_res_in = new float[LFP_n_per_dim_sq];
	float* LFP_recv_in = new float[LFP_n_per_dim_sq];

	int* LFP_rcount = new int[n_workers];  //how many LFP objects are processed for each processor
	int* LFP_displ = new int[n_workers];
	int* LFP_pts_per_dim = new int[LFP_n_per_dim];
	int* LFP_start = new int[n_workers];
	int* LFP_end = new int[n_workers];
	double LFP_radius;
	get_all_LFP_points(LFP_pts_per_dim, LFP_radius, LFP_n_per_dim);
	initialize_LFP(LFParray,  LFP_pts_per_dim, LFP_radius, LFP_n_per_dim);
	LFP_jobs_partition(LFP_start, LFP_end, LFP_displ, LFP_rcount, n_workers, LFP_n_per_dim_sq);

	// K setup
	cout<<"Before K setup, rank = "<<rank<<endl;
	int K_array_length = XDIM*YDIM;
	int K_per_dim = XDIM;
	LFP* K_array = new LFP[K_array_length];
	float* K_res_L23 = new float[K_array_length];
	float* K_recv_L23 = new float[K_array_length];
	//float* K_res_L4 = new float[K_array_length];
	//float* K_recv_L4 = new float[K_array_length];
  //float* K_res_L5 = new float[K_array_length];
	//float* K_recv_L5 = new float[K_array_length];
	//float* K_res_L6 = new float[K_array_length];
	//float* K_recv_L6 = new float[K_array_length];
	float K_res_L4_dummy, K_res_L5_dummy, K_res_L6_dummy;
	float* K_res_in = new float[K_array_length]; // used
	float* K_recv_in = new float[K_array_length]; // used


	int* K_rcount = new int[n_workers];  //how many LFP objects are processed for each processor
	int* K_displ = new int[n_workers];
	int* K_pts_per_dim = new int[K_per_dim];
	int* K_start = new int[n_workers];
	int* K_end = new int[n_workers];
	initialize_Ca(K_array, XDIM*YDIM);
	LFP_jobs_partition(K_start, K_end, K_displ, K_rcount, n_workers, K_array_length);


// Ca setup
	cout<<"Before Ca setup, rank = "<<rank<<endl;
	int Ca_array_length = XDIM*YDIM;
	int Ca_per_dim = XDIM;
	LFP* Ca_array = new LFP[Ca_array_length];
	float* Ca_res_L23 = new float[Ca_array_length];
	float* Ca_recv_L23 = new float[Ca_array_length];
	//float* K_res_L4 = new float[K_array_length];
	//float* K_recv_L4 = new float[K_array_length];
  //float* K_res_L5 = new float[K_array_length];
	//float* K_recv_L5 = new float[K_array_length];
	//float* K_res_L6 = new float[K_array_length];
	//float* K_recv_L6 = new float[K_array_length];
	float Ca_res_L4_dummy, Ca_res_L5_dummy, Ca_res_L6_dummy;
	float* Ca_res_in = new float[Ca_array_length]; // used
	float* Ca_recv_in = new float[Ca_array_length]; // used

	int* Ca_rcount = new int[n_workers];  //how many LFP objects are processed for each processor
	int* Ca_displ = new int[n_workers];
	int* Ca_pts_per_dim = new int[Ca_per_dim];
	int* Ca_start = new int[n_workers];
	int* Ca_end = new int[n_workers];
	initialize_Ca(Ca_array, XDIM*YDIM);
	LFP_jobs_partition(Ca_start, Ca_end, Ca_displ, Ca_rcount, n_workers, Ca_array_length);

	// KB setup
	cout<<"Before KB setup, rank = "<<rank<<endl;
	int KB_array_length = XDIM*YDIM;
	int KB_per_dim = XDIM;
	LFP* KB_array = new LFP[KB_array_length];
	float* KB_res_L23 = new float[KB_array_length];
	float* KB_recv_L23 = new float[KB_array_length];
	//float* K_res_L4 = new float[K_array_length];
	//float* K_recv_L4 = new float[K_array_length];
  //float* K_res_L5 = new float[K_array_length];
	//float* K_recv_L5 = new float[K_array_length];
	//float* K_res_L6 = new float[K_array_length];
	//float* K_recv_L6 = new float[K_array_length];
	float KB_res_L4_dummy, KB_res_L5_dummy, KB_res_L6_dummy;
	float* KB_res_in = new float[KB_array_length]; // used
	float* KB_recv_in = new float[KB_array_length]; // used

	int* KB_rcount = new int[n_workers];  //how many LFP objects are processed for each processor
	int* KB_displ = new int[n_workers];
	int* KB_pts_per_dim = new int[Ca_per_dim];
	int* KB_start = new int[n_workers];
	int* KB_end = new int[n_workers];
	initialize_Ca(KB_array, XDIM*YDIM);
	LFP_jobs_partition(KB_start, KB_end, KB_displ, KB_rcount, n_workers, KB_array_length);





  // fill synaptic matrix; (Testing only)
  //syn[0][0]=0;
  //syn[1][1]=0;
  //syn[0][1]=0.2;
  //syn[1][0]=0.5;
  //gi<-j
  //ref: http://www.ibiblio.org/pub/languages/fortran/append-c.html
  //In this scenario, all conductances are e-e but can be to different groups of e-cells.
  //syn[0][0]=atof(argv[3]); // homo
  //
  //
  
  // Reading neuron info (neuron and synapse)
	seed = new long*[16]; // max 16 openmp threads
	for (int i=0; i<16;i++){
		*(seed+i)=new long[16]; // This is openmp padding, no real use other than allegedly speeding up the calculation
	}
  //read_other_info(infile_cp, &num_neurons, *(seed+0)+0, &start_int_time);
	cout<<"Reading number of neurons, rank ="<<rank<<endl;
	read_num_neurons(infile_cp, &num_neurons);
	read_seed(infile_cp, *(seed+0)+0);
	for (int i=0; i<16; i++){
		seed[i][0]=seed[0][0]-i;
		//seed[i][0]=seed[0][0]-0; // This is wrong!
		cout<<"Seed initialization: MPI process "<<rank<<" i = "<<i<<" seed = "<<seed[i][0]<<endl;
	}

  cout<<"num_neurons is "<<num_neurons<<endl;
  
  // Divide the job to different workers
  assert(num_neurons%(n_workers)==0);
  cout<<"n_workers is "<<n_workers<<endl;
  int quotient = num_neurons/(n_workers);
  int start_neuron = (rank)*quotient;
  int end_neuron = (rank+1)*quotient;
  
  assert(num_neurons>0);
  var_vec = new Fmap[quotient];
  par_vec = new Fmap[quotient];

#ifdef __HETERO__
		assign_neurons_cressman_2(var_vec, par_vec, quotient, num_neurons, n_workers, rank, rs_iext, ib_iext, fs_iext, syn_delay, syn_delay_var);
#else
		assign_neurons_cressman(var_vec, par_vec, quotient, num_neurons, n_workers, rank, rs_iext, ib_iext, fs_iext, syn_delay, syn_delay_var);
#endif

	// allocate memory for diffusion array
	diffusion_array = new Diffusion[quotient];
	// initialize diffusion array
	for (int i=0; i<quotient; i++){
		Diffusion temp(start_neuron + i);
		*(diffusion_array + i) = temp;
		(diffusion_array + i)->find_diffusion_neighbours(); // finding the cell neighbours indicies for DIFFUSION
		//cout<<"After finding diffusion neighbours"<<endl;
	}
    
	// allocate memory for synaptic matrix
	cout<<"Now getting to synapse matrix init, rank="<<rank<<endl;
	syn_size = new int[quotient];
	read_syn_size(syn_size, quotient, connsize_file_stream); // reading syn_size
	connsize_file_stream.close();
  syn_sorted=new double*[quotient];
	syn_location=new int*[quotient];
	syn_counter=new int[quotient];
  for (int i=0; i<quotient; i++){
		//cout<<"syn sorted i is "<<i<<endl;
		//cout<<"syn size i is"<<syn_size[i]<<endl;
    *(syn_sorted+i) = new double[syn_size[i]];
		*(syn_location+i) = new int[syn_size[i]];
		for (int j=0; j<syn_size[i]; j++){
			syn_sorted[i][j] = 0;
			syn_location[i][j] = 0;
			//syn_counter[i] = 0;
		}
		syn_counter[i]=0;
  }
	cout<<"Can we get past synapse matrix init? rank="<<rank<<endl;

	read_connectivity_matrix_custom(syn_sorted, syn_location, syn_counter, syn_size, conn_file_stream);
	cout<<"Getting past synapse matrix init. Rank="<<rank<<endl;
	//read_connectivity_matrix(syn, infile_cp, quotient, num_neurons, n_workers, rank);
	//read_connectivity_matrix_custom(syn, connfile_cp, quotient, num_neurons, n_workers, rank);
	//for (int i=0; i<quotient; i++){
	//	for (int j=0; j<num_neurons; j++){
	//		if (syn[i][j]>0) cout<<"Rank = "<<rank<<"; syn "<<i<<"<-"<<j<<" "<<syn[i][j]<<endl;
	//	}
  //}

  conn_file_stream.close();

  //Initialize neuron (synapses already initialized in both cases)
  neurons_row = new Neuron*[quotient];
	for (size_t i=0; i<quotient; i++){
		*(neurons_row+i) = new Neuron[16];//array padding for openmp
	}

  //Synapse prepping (for more efficient computation)
  //valid_entries = new int*[quotient];
  //num_valid_entries = new int[quotient];
	esyn_array = new double[quotient]; // This is the esyn array sent to gather
	esyn_grand_array = new double[num_neurons]; // After MPI gathering, this is the array to use for calculation
  cout<<"Before init"<<endl;
  for (int i=0; i<quotient; i++){
    string neuron_ei;
    string neuron_type;
		string syn_type; (use_vmd_flag==1)?syn_type="VmD2discont":syn_type="Poissondiscont";
		float ran_num = ran2(&(seed[0][0]));
		 ((*(par_vec+i))["ei"]>0)?neuron_ei="exc":neuron_ei="inh";
		if (use_vmd_flag==0){
				//syn_type = "Poissondiscont"; // 100% of neuronal population
				if (neuron_ei=="exc"){
				(*(par_vec+i))["poisson_rate"]=poss_rate; // 2 Hz (default 0.002)
				rout<<"# poss rate is "<<	(*(par_vec+i))["poisson_rate"]<<endl;

				(*(par_vec+i))["poisson_increment"]=poss_increment; // used to be 0.2
				}
				else{
				(*(par_vec+i))["poisson_rate"]=poss_rate_inhib; // 2 Hz (default 0.002)
				rout<<"# poss rate is "<<	(*(par_vec+i))["poisson_rate"]<<endl;

				(*(par_vec+i))["poisson_increment"]=poss_increment_inhib; // used to be 0.2
				}
				(*(par_vec+i))["tau_e"]=3.0;
				(*(par_vec+i))["ge_0"]=0;
			
		}
		else if(use_vmd_flag==1){
			if (neuron_ei=="exc"){
				(*(par_vec+i))["ge_0"] = ge_zero;  // VmD mean
				//(*(par_vec+i))["tau_e"] = 3.0;
				(*(par_vec+i))["tau_e"] = 3.0;
				(*(par_vec+i))["vmd_noise"] = sqrt(2.0)*sigma_e/sqrt((*(par_vec+i))["tau_e"]); 		
				(*(par_vec+i))["ge2_0"] = gi_zero;  // VmD mean (inhibitory)
				(*(par_vec+i))["tau_e2"] = 5.0;
				//(*(par_vec+i))["tau_e2"] = 8.0;
				(*(par_vec+i))["vmd_noise2"] = sqrt(2.0)*sigma_i/sqrt((*(par_vec+i))["tau_e2"]);
						}
			else{
				(*(par_vec+i))["ge_0"] = ge_zero_i;  // VmD mean
				//(*(par_vec+i))["tau_e"] = 3.0;
				(*(par_vec+i))["tau_e"] = 3.0;
				(*(par_vec+i))["vmd_noise"] = sqrt(2.0)*sigma_e_i/sqrt((*(par_vec+i))["tau_e"]); 		
				(*(par_vec+i))["ge2_0"] = gi_zero_i;  // VmD mean (inhibitory)
				(*(par_vec+i))["tau_e2"] = 5.0;
				//(*(par_vec+i))["tau_e2"] = 8.0;
				(*(par_vec+i))["vmd_noise2"] = sqrt(2.0)*sigma_i_i/sqrt((*(par_vec+i))["tau_e2"]);
							}

		}

    //((*(par_vec+i))["nt"]>0)?neuron_type="Anderson":neuron_type="Anderson_sans_Ca";
		//((*(par_vec+i))["nt"]>0)?neuron_type="Cressman":neuron_type="Cressman";
		((*(par_vec+i))["ei"]>0)?neuron_type="Froh":neuron_type="Froh";

    (*(par_vec+i))["time_step"] = dt;
		(*(par_vec+i))["inhib_ramp_factor"]=inhib_ramp_factor;

    //cout<<"Neuron type is "<<neuron_type<<endl;
		//
		//Glial parameters
		if (neuron_ei=="exc"){
		//	(*(par_vec+i))["cressman_G_glia"] = gglia;
		//(*(par_vec+i))["cressman_k_zero_inf"] = kzeroinf;
		(*(par_vec+i))["frohlich_k_forward"] = k_for;
		(*(par_vec+i))["frohlich_k_back"] = k_back;
		(*(par_vec+i))["frohlich_Keq_pump"] = Ko_eq_pump;
		(*(par_vec+i))["frohlich_Keq_glia"] = Ko_eq_glia;
		(*(par_vec+i))["frohlich_max_pump_current"] = max_pump_current;
		(*(par_vec+i))["frohlich_glia_rise_factor"] = rise;
		(*(par_vec+i))["frohlich_b_max"] = 500;



		//Stim parameters
		(*(par_vec+i))["stim_start"] = stim_start;
		(*(par_vec+i))["stim_end"] = stim_end;
		(*(par_vec+i))["stim_strength"] = stim_strength;
		
		}
		else{
//	(*(par_vec+i))["cressman_G_glia"] = gglia_in;
//		(*(par_vec+i))["cressman_k_zero_inf"] = kzeroinf_in;
			if ((*(par_vec+i))["inhibitory_type"]==1){

		(*(par_vec+i))["frohlich_k_forward"] = k_for_in;
		(*(par_vec+i))["frohlich_k_back"] = k_back_in;
		(*(par_vec+i))["frohlich_Keq_pump"] = Ko_eq_pump_in;
		(*(par_vec+i))["frohlich_Keq_glia"] = Ko_eq_glia_in;
		(*(par_vec+i))["frohlich_max_pump_current"] = max_pump_current_in;
		(*(par_vec+i))["frohlich_glia_rise_factor"] = rise_in;
			}
			else{

		(*(par_vec+i))["frohlich_k_forward"] = k_for_in_2;
		(*(par_vec+i))["frohlich_k_back"] = k_back_in_2;
		(*(par_vec+i))["frohlich_Keq_pump"] = Ko_eq_pump_in_2;
		(*(par_vec+i))["frohlich_Keq_glia"] = Ko_eq_glia_in_2;
		(*(par_vec+i))["frohlich_max_pump_current"] = max_pump_current_in_2;
		(*(par_vec+i))["frohlich_glia_rise_factor"] = rise_in_2;
			}

		(*(par_vec+i))["frohlich_b_max"] = 500;


		// Stim parameters
		(*(par_vec+i))["stim_start"] = stim_start_in;
		(*(par_vec+i))["stim_end"] = stim_end_in;
		(*(par_vec+i))["stim_strength"] = stim_strength_in;

		}



    Neuron this_neuron(neuron_type, syn_type, neuron_ei, (*(var_vec+i)), (*(par_vec+i)));
    //cout<<"Passed initialization"<<endl;
    neurons_row[i][0] = this_neuron;
    //*(valid_entries + i) = (*(neurons_row + i)+0)->find_valid_synapse_entries(*(syn+i), num_neurons, *(num_valid_entries+i));// returning an int pointer which points to an array of valid entries
    *(esyn_array + i) = (*(neurons_row + i)+0)->prep_synapse_reversal_potential(); //Have to rewrite this function; returning the value of Esyn
		//rout<<*(esyn_array+i)<<endl;
    //cout<<"After assignment"<<endl;
   }
  delete [] var_vec;
  delete [] par_vec;

	MPI_Gather(esyn_array, quotient, MPI_DOUBLE, esyn_grand_array, quotient, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(esyn_grand_array, num_neurons, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	//for (int m=0; m<=1; m++){
		//if (rank==m){
			//for (int i=m*quotient; i<m*quotient+16; i++){
			//	cout<<"i is "<<i<<" Esyn is "<<esyn_grand_array[i]<<endl;
			//}
		//}
	//}


  // End initialization
  
  // Do integration
  double t=0;
  double printint=1.0;
  int countint=int(printint/dt); 
  double sqroot_dt=sqrt(dt);
  double consum_inhib=0; //sum of si
  double consum_excit=0;
  double this_sum_i=0;
  int count=1;
  int nowprint=0;
  
  rout<<"# Method: Euler; dt="<<dt<<"ms."<<endl;
  // Change display (May 24th, 2010)
  rout<<"# Number of neurons="<<num_neurons<<" seed="<<seed<<endl;
  
  double this_start_time = dt;// + start_int_time;
  double this_end_time = maxtime;// + start_int_time;
  //cout<<"iset dump is "<<*iset_dump<<endl;
  //int array_to_send_buffer_size = end_neuron - start_neuron;
  
  array_to_receive = new double[num_neurons];
	K_extracellular_array_to_receive = new double[num_neurons];
	Ca_intracellular_array_to_receive = new double[num_neurons];
	potentials_to_receive = new double[num_neurons];
	KB_array_to_receive = new double[num_neurons];
	//Ca_array_to_receive = new float[num_neurons];
  for (int i=0; i<num_neurons; i++){
    array_to_receive[i]=0;
		potentials_to_receive[i]=0;
		K_extracellular_array_to_receive[i]=0;
		//Ca_array_to_receive[i]=0;
  }
  array_to_send = new double[quotient];
	K_extracellular_array_to_send = new double[quotient];
	potentials_to_send = new double[quotient];
	Ca_intracellular_array_to_send = new double[quotient];
	KB_array_to_send = new double[quotient];
	//Ca_array_to_send = new float[quotient];
	//for (int i=0; i<quotient; i++){
  //*(array_to_send+i) = new double[16]; //padding for open mp
	//}
  MPI_Barrier(MPI_COMM_WORLD);

  //cout<<"This is processor "<<rank<<endl;
  //cout<<"Size of array is "<<quotient<<endl;
  //omp_set_dynamic(0);     // Explicitly disable dynamic teams
  //omp_set_num_threads(16);

	neurons_row_transposed = new Neuron*[quotient];
	for (size_t i=0;i<quotient; i++){
		*(neurons_row_transposed+i) = *(neurons_row+i)+0;
		//cout<<*(neurons_row+i)+0<<endl;
	} 
	compustarttime = time(&thistime);


#pragma omp parallel firstprivate(t,thid)  num_threads(1)
  {
	thid = omp_get_thread_num(); // getting the thread ID for open mp
  cout<<"This starts parallel"<<endl;
  for (t=this_start_time; t<this_end_time; t+=dt){
#pragma omp for
	for (size_t i=0; i<quotient; i++){
		//for (size_t i=0; i<1; i++){

		//	neurons_row[i][0].get_all_derivatives(*(syn+i), num_neurons, *(neurons_row_transposed+0), array_to_receive, *(valid_entries+i), *(num_valid_entries+i));  
		//	neurons_row[i][0].get_all_derivatives(*(syn+i), num_neurons, esyn_grand_array, array_to_receive, *(valid_entries+i), *(num_valid_entries+i));
		//neurons_row[i][0].get_all_derivatives(*(syn+i), num_neurons, esyn_grand_array, array_to_receive, *(valid_entries+i), *(num_valid_entries+i));
	neurons_row[i][0].cellVariables::current_time=t;
	neurons_row[i][0].get_all_derivatives(*(syn_sorted+i), *(syn_location+i), *(syn_size+i), esyn_grand_array, array_to_receive);
		//if (i==15){
	//	for (int m=0; m<(*(syn_size+i)); m++){
	//		cout<<*(*(syn_sorted+i)+m)<<endl;
	//		cout<<*(*(syn_location+i)+m)<<endl;
	//	}
	//}
			  }
#pragma omp for
for (size_t i=0; i<quotient; i++){
	//for (size_t i=0; i<1; i++){
		//cout<<"Now doing i "<<i<<endl;

		neurons_row[i][0].cellVariables::K_extracellular_diffusion_term = lapacian(diffusion_array[i], K_extracellular_array_to_receive, start_neuron + i);
		//neurons_row[i][0].cellVariables::K_extracellular_diffusion_term = 0;
	 	neurons_row[i][0].integrate_neuron(dt, sqroot_dt, &(seed[thid][0]), thid);
		neurons_row[i][0].integrate_synapse(dt, sqroot_dt, &(seed[thid][0]), thid);   
	  neurons_row[i][0].update();
	//	*(esyn_array + i) = (*(neurons_row + i)+0)->prep_synapse_reversal_potential(); 

	  //cout<<"Before assembling array to send i is "<<i<<endl;
		//*(*(array_to_send+i)+0)=(*(neurons_row+i)+0)->cellVariables::synapse_i_delayed;		
		*(array_to_send+i)=(*(neurons_row+i)+0)->cellVariables::synapse_i_delayed;	
		*(potentials_to_send+i)=(*(neurons_row+i)+0)->cellVariables::membrane_potential;
		//if (rank==3)cout<<"Potential to send "<<i<<" is "<<*(potentials_to_send+i)<<endl;
		*(K_extracellular_array_to_send+i)=(*(neurons_row+i)+0)->cellVariables::K_extracellular;
		*(Ca_intracellular_array_to_send+i)=(*(neurons_row+i)+0)->cellVariables::Ca_intracellular;
		*(KB_array_to_send+i)=500.00-(*(neurons_row+i)+0)->cellVariables::frohlich_buffer;


	//	((*(neurons_row+i)+0)->display_type()=="Anderson")?(*(Ca_array_to_send+i))=(*(neurons_row+i)+0)->cellVariables::Ca_C6:(*(Ca_array_to_send+i))=0;
		//*(Ca_array_to_send+i) = 0;
		if ((*(neurons_row+i)+0)->spiked==1) network_spiked = 1;
		//if ((rank==0 && i==0)) rout<<t<<"\t"<<(*(neurons_row+i)+0)->cellVariables::membrane_potential<<"\t"<<(*(neurons_row+i)+0)->cellVariables::synapse_e<<"\n";
		//else if (rank==0) rout<<(*(neurons_row+i)+0)->cellVariables::membrane_potential<<"\t";
		//if (i==0) rout<<i<<"\t"<<t<<"\t"<<(*(neurons_row+i)+0)->cellVariables::membrane_potential<<"\t"<<(*(neurons_row+i)+0)->cellVariables::potassium_channel<<"\t"<<(*(neurons_row+i)+0)->cellVariables::sodium_channel<<"\t"<<(*(neurons_row+i)+0)->cellVariables::Ca_intracellular<<"\t"<<(*(neurons_row+i)+0)->cellVariables::K_extracellular<<"\t"<<(*(neurons_row+i)+0)->cellVariables::Na_intracellular<<"\t";
		//if (i==15) rout<<i<<"\t"<<t<<"\t"<<(*(neurons_row+i)+0)->cellVariables::membrane_potential<<"\t"<<(*(neurons_row+i)+0)->cellVariables::potassium_channel<<"\t"<<(*(neurons_row+i)+0)->cellVariables::sodium_channel<<"\t"<<(*(neurons_row+i)+0)->cellVariables::Ca_intracellular<<"\t"<<(*(neurons_row+i)+0)->cellVariables::K_extracellular<<"\t"<<(*(neurons_row+i)+0)->cellVariables::Na_intracellular<<"\t";

		//if (i==0) rout<<i<<"\t"<<t<<"\t"<<(*(neurons_row+i)+0)->cellVariables::membrane_potential<<"\t"<<(*(neurons_row+i)+0)->cellVariables::potassium_channel<<"\t"<<(*(neurons_row+i)+0)->cellVariables::sodium_channel<<"\t"<<(*(neurons_row+i)+0)->cellVariables::Ca_intracellular<<"\t"<<(*(neurons_row+i)+0)->cellVariables::K_extracellular<<"\t"<<(*(neurons_row+i)+0)->cellVariables::Na_intracellular<<"\t"<<(*(neurons_row+i)+0)->cellVariables::VK_var<<"\t"<<(*(neurons_row+i)+0)->cellVariables::VNa_var<<"\t";
//if (i==0) rout<<i<<"\t"<<t<<"\t"<<(*(neurons_row+i)+0)->cellVariables::membrane_potential<<"\t"<<(*(neurons_row+i)+0)->cellVariables::potassium_channel<<"\t"<<(*(neurons_row+i)+0)->cellVariables::sodium_channel<<"\t"<<(*(neurons_row+i)+0)->cellVariables::Ca_intracellular<<"\t"<<(*(neurons_row+i)+0)->cellVariables::K_extracellular<<"\t"<<(*(neurons_row+i)+0)->cellVariables::Na_intracellular<<"\t"<<(*(neurons_row+i)+0)->cellVariables::frohlich_buffer<<"\t";
//if (i==0) rout<<i<<"\t"<<t<<"\t"<<(*(neurons_row+i)+0)->cellVariables::membrane_potential<<"\t"<<(*(neurons_row+i)+0)->cellVariables::glia_current<<"\t"<<(*(neurons_row+i)+0)->cellVariables::pump_current<<"\t"<<(*(neurons_row+i)+0)->cellVariables::potassium_current<<"\t"<<(*(neurons_row+i)+0)->cellVariables::delta_K<<"\t";;

	}
	rout<<"\n";
	//if (rank==0)rout<<endl;
	//#pragma omp barrier
	//cout<<"Before mpi allgather"<<endl;
	//	MPI_Allgather(array_to_send, quotient, MPI_DOUBLE, array_to_receive, quotient, MPI_DOUBLE, MPI_COMM_WORLD);
#pragma omp master
			{	

			//MPI_Barrier(MPI_COMM_WORLD);
			//MPI_Gather(esyn_array, quotient, MPI_DOUBLE, esyn_grand_array, quotient, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			//MPI_Bcast(esyn_grand_array, num_neurons, MPI_DOUBLE, 0, MPI_COMM_WORLD);


			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Gather(array_to_send, quotient, MPI_DOUBLE, array_to_receive, quotient, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(array_to_receive, num_neurons, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			// Let's try putting the K gather inside the 1/10 loop
			MPI_Gather(K_extracellular_array_to_send, quotient, MPI_DOUBLE, K_extracellular_array_to_receive, quotient, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(K_extracellular_array_to_receive, num_neurons, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			// do Ca
			MPI_Gather(Ca_intracellular_array_to_send, quotient, MPI_DOUBLE, Ca_intracellular_array_to_receive, quotient, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(Ca_intracellular_array_to_receive, num_neurons, MPI_DOUBLE, 0, MPI_COMM_WORLD);

			// do KB
			MPI_Gather(KB_array_to_send, quotient, MPI_DOUBLE, KB_array_to_receive, quotient, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(KB_array_to_receive, num_neurons, MPI_DOUBLE, 0, MPI_COMM_WORLD);


			//Print spike statistics for every step;
			for (size_t i=0; i<quotient; i++){
				 if ((*(neurons_row+i)+0)->spiked==1) sout<<start_neuron+i<<"\t"<<t<<endl;
				}
			network_spiked=0;
	
			(count%countint==0)?nowprint=1:nowprint=0;
					if (nowprint==1){
				//cout<<"Now do LFP and Ca"<<endl;
				// do LFP
				MPI_Gather(potentials_to_send, quotient, MPI_DOUBLE, potentials_to_receive, quotient, MPI_DOUBLE, 0, MPI_COMM_WORLD);				
				MPI_Bcast(potentials_to_receive, num_neurons, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
							


							
				//  do Ca
				//	MPI_Gather(Ca_array_to_send, quotient, MPI_FLOAT, Ca_array_to_receive, quotient, MPI_FLOAT, 0, MPI_COMM_WORLD);				
				//	MPI_Bcast(Ca_array_to_receive, num_neurons, MPI_FLOAT, 0, MPI_COMM_WORLD);

				for (int m=LFP_start[rank]; m<LFP_end[rank]; m++){
				//sout<<"LFP_start["<<rank<<"]="<<LFP_start[rank]<<endl;
				//sout<<"LFP_end["<<rank<<"]="<<LFP_end[rank]<<endl;
					(LFParray+m)->calculate_LFP(potentials_to_receive, LFP_res_L23[m], LFP_res_L4[m], LFP_res_L5[m], LFP_res_L6[m], LFP_res_in[m]);
					//LFP_res_L23[m]=rank;
					//sout<<"L23 result is "<<LFP_res_L23[m]<<endl;
					//sout<<"LFP_rcount is "<<LFP_rcount[rank]<<endl;
					//sout<<"displ is "<<LFP_displ[rank]<<endl;
					//sout<<endl;
					//LFP_res[m] = m;
					//		lfp23out<<LFP_res<<"\t";
				}

				for (int m=K_start[rank]; m<K_end[rank]; m++){
					
					(K_array+m)->calculate_LFP(K_extracellular_array_to_receive, K_res_L23[m], K_res_L4_dummy, K_res_L5_dummy, K_res_L6_dummy, K_res_in[m]);

				}
				for (int m=Ca_start[rank]; m<Ca_end[rank]; m++){
										
					(Ca_array+m)->calculate_LFP(Ca_intracellular_array_to_receive, Ca_res_L23[m], Ca_res_L4_dummy, Ca_res_L5_dummy, Ca_res_L6_dummy, Ca_res_in[m]);

				}
				for (int m=KB_start[rank]; m<KB_end[rank]; m++){
										
					(KB_array+m)->calculate_LFP(KB_array_to_receive, KB_res_L23[m], KB_res_L4_dummy, KB_res_L5_dummy, KB_res_L6_dummy, KB_res_in[m]);

				}




				MPI_Barrier(MPI_COMM_WORLD);
				MPI_Gatherv(LFP_res_L23+LFP_start[rank], LFP_rcount[rank], MPI_FLOAT, LFP_recv_L23, LFP_rcount, LFP_displ, MPI_FLOAT, 0, MPI_COMM_WORLD);
			//	MPI_Gatherv(LFP_res_L4+LFP_start[rank], LFP_rcount[rank], MPI_FLOAT, LFP_recv_L4, LFP_rcount, LFP_displ, MPI_FLOAT, 0, MPI_COMM_WORLD);
			//	MPI_Gatherv(LFP_res_L5+LFP_start[rank], LFP_rcount[rank], MPI_FLOAT, LFP_recv_L5, LFP_rcount, LFP_displ, MPI_FLOAT, 0, MPI_COMM_WORLD);
			//	MPI_Gatherv(LFP_res_L6+LFP_start[rank], LFP_rcount[rank], MPI_FLOAT, LFP_recv_L6, LFP_rcount, LFP_displ, MPI_FLOAT, 0, MPI_COMM_WORLD);
				MPI_Gatherv(LFP_res_in+LFP_start[rank], LFP_rcount[rank], MPI_FLOAT, LFP_recv_in, LFP_rcount, LFP_displ, MPI_FLOAT, 0, MPI_COMM_WORLD);

				MPI_Gatherv(K_res_L23+K_start[rank], K_rcount[rank], MPI_FLOAT, K_recv_L23, K_rcount, K_displ, MPI_FLOAT, 0, MPI_COMM_WORLD);
				//MPI_Gatherv(Ca_res_L4+Ca_start[rank], Ca_rcount[rank], MPI_FLOAT, Ca_recv_L4, Ca_rcount, Ca_displ, MPI_FLOAT, 0, MPI_COMM_WORLD);
				//MPI_Gatherv(Ca_res_L5+Ca_start[rank], Ca_rcount[rank], MPI_FLOAT, Ca_recv_L5, Ca_rcount, Ca_displ, MPI_FLOAT, 0, MPI_COMM_WORLD);
				//MPI_Gatherv(Ca_res_L6+Ca_start[rank], Ca_rcount[rank], MPI_FLOAT, Ca_recv_L6, Ca_rcount, Ca_displ, MPI_FLOAT, 0, MPI_COMM_WORLD);
				MPI_Gatherv(K_res_in+K_start[rank], K_rcount[rank], MPI_FLOAT, K_recv_in, K_rcount, K_displ, MPI_FLOAT, 0, MPI_COMM_WORLD);

				MPI_Gatherv(Ca_res_L23+Ca_start[rank], Ca_rcount[rank], MPI_FLOAT, Ca_recv_L23, Ca_rcount, Ca_displ, MPI_FLOAT, 0, MPI_COMM_WORLD);
							MPI_Gatherv(Ca_res_in+Ca_start[rank], Ca_rcount[rank], MPI_FLOAT, Ca_recv_in, Ca_rcount, Ca_displ, MPI_FLOAT, 0, MPI_COMM_WORLD);


				MPI_Gatherv(KB_res_L23+KB_start[rank], KB_rcount[rank], MPI_FLOAT, KB_recv_L23, KB_rcount, KB_displ, MPI_FLOAT, 0, MPI_COMM_WORLD);
								MPI_Gatherv(KB_res_in+KB_start[rank], KB_rcount[rank], MPI_FLOAT, KB_recv_in, KB_rcount, KB_displ, MPI_FLOAT, 0, MPI_COMM_WORLD);


				if (rank==0){
					//lfp4out<<esyn_grand_array[4048]<<"\t"<<esyn_grand_array[4062];
					for (int i=0; i<LFP_n_per_dim_sq; i++){
						lfp23out<<LFP_recv_L23[i]<<"\t";
							lfpinout<<LFP_recv_in[i]<<"\t";


						//lfp4out<<LFP_recv_L4[i]<<"\t";
						//lfp5out<<LFP_recv_L5[i]<<"\t";
						//lfp6out<<LFP_recv_L6[i]<<"\t";
						
										}
					for (int i=0; i<K_array_length; i++){
						K23out<<K_recv_L23[i]<<"\t";
						Kinout<<K_recv_in[i]<<"\t";
					
					//	Ca4out<<Ca_recv_L4[i]<<"\t";
					//	Ca5out<<Ca_recv_L5[i]<<"\t";
					//	Ca6out<<Ca_recv_L6[i]<<"\t";
					}
					for (int i=0; i<Ca_array_length; i++){
						Ca23out<<Ca_recv_L23[i]<<"\t";
						Cainout<<Ca_recv_in[i]<<"\t";
					
					//	Ca4out<<Ca_recv_L4[i]<<"\t";
					//	Ca5out<<Ca_recv_L5[i]<<"\t";
					//	Ca6out<<Ca_recv_L6[i]<<"\t";
					}
					for (int i=0; i<KB_array_length; i++){
						KB23out<<KB_recv_L23[i]<<"\t";
						KBinout<<KB_recv_in[i]<<"\t";
					
					//	Ca4out<<Ca_recv_L4[i]<<"\t";
					//	Ca5out<<Ca_recv_L5[i]<<"\t";
					//	Ca6out<<Ca_recv_L6[i]<<"\t";
					}


				lfp23out<<"\n";
				lfp4out<<"\n";
				//lfp5out<<"\n";
				//lfp6out<<"\n";
				lfpinout<<"\n";

				K23out<<"\n";
				Kinout<<"\n";
				Ca23out<<"\n";
				Cainout<<"\n";
				KB23out<<"\n";
				KBinout<<"\n";

			//
			//	Ca4out<<"\n";
			//	Ca5out<<"\n";
			//	Ca6out<<"\n";

				//lfp23out<<*(LFP_res+LFP_start[rank])<<endl;
				}
				//	lfp23out<<endl;
				//}
			}
			count++;
			MPI_Barrier(MPI_COMM_WORLD);
      }
#pragma omp barrier
    }
  }

cout<<"This ends parallel"<<endl;
   
  outfile_cp.close();

  //for (int i=0; i<quotient; i++){
   // delete [] *(syn+i);
  //}
  //delete [] syn;
  //delete [] neurons_row;
  rout<<"# Seed is now:"<<seed[0][0]<<endl;
  rout<<"# Time taken:"<<time(&thistime)-starttime<<" second(s)."<<endl;
	rout<<"# Time taken for computation:"<<time(&thistime)-compustarttime<<" second(s)."<<endl;
  rout.close();
	sout.close();
	if (rank==0){
		lfp23out.close();
		lfp4out.close();
		lfp5out.close();
		lfp6out.close();
		lfpinout.close();
		K23out.close();
		K4out.close();
		K5out.close();
		//K6out.close();
		Kinout.close();


	}
	MPI_Finalize();
  return 0;
}

// test code
	//for (int i=0; i<quotient; i++){
	//	for (int j=0; j<num_neurons; j++){
		//	if (syn[i][j]==0)continue;
		//	else{
		//		rout<<"MPI task is "<<rank<<" i is "<<i<<" j is "<<j<<" Value is "<<syn[i][j]<<endl;
		//		int xf, yf, layerf, celltypef, cellnumberf, cnf;
		//		int xt, yt, layert, celltypet, cellnumbert, cnt;
		//		reverse_indexing_neurons(i+rank*quotient, xt, yt, cnt);
		//		reverse_indexing_neurons(j, xf, yf, cnf);
		//		custom_convert_cn_to_lcc(layerf, celltypef, cellnumberf, cnf);
		//		custom_convert_cn_to_lcc(layert, celltypet, cellnumbert, cnt);
		//		rout<<" i is "<<xt<<"\t"<<yt<<"\t"<<layert<<"\t"<<celltypet<<"\t"<<cellnumbert<<endl;
		//		rout<<" j is "<<xf<<"\t"<<yf<<"\t"<<layerf<<"\t"<<celltypef<<"\t"<<cellnumberf<<endl;
		//	}

		//}
	//}
	//return 0;
	//
	//
  //outfile_cp.open(cp_file_to_write);
  //done_dump_gasdev_flag = 0;
  //done_dump_ran_flag = 0;      
  //ran2(&seed_hetero, NULL, NULL, NULL, idum2_dump,iy_dump, iv_dump, done_init_ran_flag,  done_dump_ran_flag); //Here we are just dumping the internal state to the ..._dump variables.  We don't care what the actualy random number output is. We use &seed_hetero because we don't want to destroy the state of &seed
  //gasdev(&seed_hetero, NULL, NULL, iset_dump, gset_dump, NULL, NULL, NULL, NULL, NULL, NULL, done_init_gasdev_flag,  done_init_ran_flag, done_dump_gasdev_flag, done_dump_ran_flag);//same as above
  //write_checkpoint_data(outfile_cp, &num_neurons, &seed, &this_end_time,  iset_dump, gset_dump, idum2_dump, iy_dump, iv_dump, iv_size, syn, neurons_row);

