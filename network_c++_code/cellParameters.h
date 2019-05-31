#include <unordered_map>
#include <string>
using namespace std;

#ifndef CELLPARAMETERS_H
#define CELLPARAMETERS_H
class cellParameters{
 	public:
  	typedef unordered_map<string,double> Fmap;
  	cellParameters();
  	cellParameters(Fmap, string, string);
  	double get_values(string name);
  	Fmap get_hash();
	cellParameters& operator=(const cellParameters& right_side);
	private:
		void init_discont(Fmap&);
		void init_wb(Fmap&);
		void init_Iz(Fmap&);
		void init_abcsimple(Fmap&);
		void init_kc(Fmap&);
		void init_vmd(Fmap&);
		void init_poisson(Fmap&);
#if ANDERSON
		void init_anderson(Fmap&);
		void init_anderson_sans_Ca(Fmap&);
#endif
		void init_ahmed_modified(Fmap&);
		void init_cressman(Fmap&);
		void class_init(Fmap&, string, string);
		void init_frohlich(Fmap&);

 	public:
		string cell_type;
		string synapse_type;
 	 	Fmap cell_pars;
		double membrane_potential;
		double external_current;
		double reversal_potential_e;
		double reversal_potential_i;
		double synapse_potential;
		double tau_i;
		double synapse_delay;
		double current_noise;
		double time_step;
		double rate_depression;
		double tau_recovery;
		double inhib_ramp_factor;

		//wb
		double leak_conductance;
		double leak_potential;
		double sodium_conductance;
		double sodium_potential;
		double potassium_conductance;
		double potassium_potential;
		double sodium_half_activation;
		double sodium_tau;
		double alpha;

		//kc
		double tau_n;
	
		//abc simple
		double Vsn;
		double Vreset;
		double Vmax;
		double Rreset;
		double a;
		double bsn;
		double c;
		double d;
		double q;
		double l;

		//Ahmed modified
		double k;
		double k_after;
		double Vchangek;
		double Vtshift;
		double d_step;
		double b;

		//Iz
		double Vr;
#if ANDERSON
		//Anderson
		double aA;
		double vhalfA;
		double aX;	
		double vhalfX;
		double aB;
		double vhalfB;
		double aW;
		double vhalfW;
		double am;
		double vhalfm;
		double lambda;
		double tauX;
		double tauB;
		double VNa;
		double VK;
		double VL;
		double gNa;
		double gK;
		double gL;
		double gA;
		double ei;
		double nt;

		double VCa;
		double gKCa;
		double zCa;
		double Kc;
		double Kd;
		double PCa;
		double C0;
		double Rp;
		double Kp;
		double b;
		double f;
		double B0;
		double D;
		double length0;
		double length1;
		double length2;
		double length3;
		double length4;
		double length5;
		double length6;
#endif

		//cressman
		double cressman_beta;
		double cressman_rho;
		double cressman_D;
		double cressman_G_glia;
		double cressman_varepsilon;
		double cressman_k_zero_inf;
		double gL;
		double gNa;
		double gK;
		double gCa;
		double gAHP;
		double gKL;
		double gNaL;
		double gClL;
		double VCl;
		double VCa;
		double Cl_intracellular;
		double Cl_extracellular;

		// Frohlich
		double frohlich_k_forward;
		double frohlich_k_back;
		double frohlich_b_max;
		double frohlich_buffer;
		double frohlich_max_pump_current;
		double frohlich_Keq_pump;
		double frohlich_Keq_glia;
		double frohlich_glia_rise_factor;
		
		//Poisson
		double poisson_rate;
		double poisson_increment;
	
		//vmd
		double ge_0; //ge
		double tau_e;
		double vmd_noise;
		double ge2_0; //gi
		double tau_e2;
		double vmd_noise2;

		//stim
		double stim_strength;
		double stim_start;
		double stim_end;

		//discont
		double Vopens;
		double synapse_i_max;	
}
;
#endif
