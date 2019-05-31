#include <unordered_map>
#include <string>
//#include <utility>

#ifndef CELLVARIABLES_H
#define CELLVARIABLES_H
using namespace std;

class cellVariables{
 public:
  typedef unordered_map<string,double> Fmap;
  cellVariables();
  cellVariables(Fmap, string);
  double get_values(string);
  Fmap get_hash();
	cellVariables& operator=(const cellVariables& right_side);

	private:
#if ANDERSON
	void init_anderson_sans_Ca(Fmap&);
	void init_anderson(Fmap&);
#endif
	void init_generic(Fmap&);
	void init_cressman(Fmap&);
	void class_init(Fmap&, string);
	void init_frohlich(Fmap&);

 public:
	string cell_type;
  Fmap cell_vars;
	double membrane_potential;
	double synapse_e;
	double synapse_e2;
	double synapse_i;
	double synapse_i_delayed;
	double synapse_i_max_current;
	double current_time;

	//generic
	double potassium_channel;
	double sodium_channel;
#if ANDERSON
	//Anderson sans Ca
	double W;
	double X;
	double B;

	//Anderson
	double Ca_U0;
	double Ca_U1;
	double Ca_U2;
	double Ca_U3;
	double Ca_U4;
	double Ca_U5;
	double Ca_U6;
	double Ca_C0;
	double Ca_C1;
	double Ca_C2;
	double Ca_C3;
	double Ca_C4;
	double Ca_C5;
	double Ca_C6;
	double ICa;
	double ICa_next;
#endif

	//Cressman
	double Na_extracellular;
	double Na_intracellular;
	double K_extracellular;
	double K_intracellular;
	//double Cl_extracellular;
	//double Cl_intracellular;
	double Ca_intracellular;
	double VNa_var;
	double VK_var;
	//double VCl_var;
	double VL_var;
	double potassium_current;
	double sodium_current;
	double K_extracellular_diffusion_term;

	//Frohlich
	double frohlich_buffer;
	double glia_current;
	double pump_current;
	double delta_K;
	double tau_now_Ca;
}
;
#endif
