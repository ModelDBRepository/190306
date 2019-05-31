//#ifdef NEURON_H
//#define NEURON_H
//#include <omp.h>
#include <unordered_map>
#include <string>
#include "cellVariables.h"
#include "cellParameters.h"
using namespace std;
class Neuron : public cellVariables, public cellParameters
{
 public:
  typedef unordered_map<string, double> Fmap;
  typedef double (**FuncPtr)(cellVariables&, cellParameters&);
  typedef double (*DFuncPtr)(cellVariables&, cellParameters&);
  Neuron(); // default constructor
  Neuron(string type_neuron,string type_synapse, string i_or_e, Fmap cell_vars, Fmap cell_pars); // constructor
  Neuron(const Neuron& right_side); // copy costructor
  ~Neuron(); // destructor
  Neuron& operator=(const Neuron& right_side); // assignment operator
  void get_all_derivatives(double* synaptic_matrix_row, int number_of_neurons, Neuron* other_neurons); // getting the derivative numbers into the appropriate spaces
  void get_all_derivatives(double si, double syn_strength); // getting the derivative numbers into the appropriate spaces
	void get_all_derivatives(double* synaptic_matrix_row, int number_of_neurons, Neuron* other_neurons, double* sv); // externally supplied synaptic variables (for MPI)
	void get_all_derivatives(double* synaptic_matrix_row, int number_of_neurons, Neuron* other_neurons, double* sv, int* valid_entries, int); //ditto + prepped synapse for efficient looping
	void get_all_derivatives(double* synaptic_matrix_row, int number_of_neurons, double* prepped_esyn_array, double* sv, int* valid_entries, int); //ditto + prepped synapse AND esyn for efficient looping
	void get_all_derivatives(double* synaptic_sorted_row, int* synaptic_location, int synaptic_size, double* prepped_esyn_array, double* sv); //ditto + even more prepping
  void integrate_neuron(const double&, const double&, long*, int*, float*, int*, float*, long*, long*, long*, long*, long*, long*, int&, int&, int&, int&);
	void integrate_neuron(const double&, const double&, long*, int*, float*, int*, float*, long*, long*, long*, long*, long*, long*, int&, int&, int&, int&, const int&);// for openmp threads
  void integrate_neuron(const double&, const double&, long*); // integrates the neuron and put the results into appropriate temporary spaces
	void integrate_neuron(const double&, const double&, long*, const int&); // for open mp threads 
  void integrate_synapse(const double&, const double&, long*, int*, float*, int*, float*, long*, long*, long*, long*, long*, long*, int&, int&, int&, int&);
void integrate_synapse(const double&, const double&, long*, int*, float*, int*, float*, long*, long*, long*, long*, long*, long*, int&, int&, int&, int&, const int&);// for open mp

  void integrate_synapse(const double&, const double&, long*); // integrates the synapse and put the results into appropriate temporary spaces
	void integrate_synapse(const double&, const double&, long*, const int&);// for open mp threads
  void update_synapse(); // updates a particular neuron's associated synapse information into cellVariables
  void update_neuron(); // updates a particular neuron's information into cellVariables
	void update(); // wrapper function for update_synapse and update neuron
  double display_current_balance_derivative(int); // Added October 22nd, 2010.
	string display_type();
	int* find_valid_synapse_entries(double* synapse_matrix_row, int number_of_neurons, int&); //synapse prepping routine for efficient looping
	double prep_synapse_reversal_potential(); //same as above
	int spiked; // flag to tell whether a neuron has just spiked (applies only to discontinuous synapses)



 private:
  FuncPtr current_balance_equations;  //pointer to array of functions (current balance)
  FuncPtr synapse_equations; //pointer to array of functions (synapse)
  double* current_balance_readings; // store temp integrated results
  double* synapse_readings;
  double* current_balance_derivatives; // store derivatives for integration
  double* synapse_derivatives;
	double* synapse_array; // array to store delayed synapse variables
	double* synapse_array_current_position; // pointer that points to the current position of the synapse array

  string type_of_neuron;
  string type_of_synapse;
  string inh_or_exc;
  int number_of_current_balance_equations;
  int number_of_synapse_equations;
	int synapse_array_length; // size of synapse array

  //Flags
  int got_current_balance_equations;
  int got_synapse_equations;
  int got_space_current_balance_results;
  int got_space_synapse_results;
  int update_neuron_happened;
	int synapse_array_assigned;
	
  void get_current_balance_equations();
  void get_space_current_balance_results();
  void get_synapse_equations();
  void get_space_synapse_results();
  void get_current_balance_derivatives_no_psps();
  void check_and_get_spaces_all_derivatives();
  void assign_all_derivatives_to_neuron(double summed_psps);
	void assign_synapse_array();
  double get_summed_psps(double* synaptic_matrix_row, int number_of_neurons, Neuron* other_neurons);
	double get_summed_psps(double* synaptic_matrix_row, int number_of_neurons, Neuron* othre_neurons, double* delayed_synapse_array);
	double get_summed_psps(double* synaptic_matrix_row, int number_of_neurons, Neuron* other_neurons, double* delayed_synapse_array, int* valid_entries, int number_of_entries);
	double get_summed_psps(double* synaptic_matrix_row, int number_of_neurons, double* array_of_esyn, double* delayed_synapse_array, int* valid_entries, int number_of_entries);
	double get_summed_psps(double* synaptic_sorted_row, int* synaptic_location, int synaptic_size, double* array_of_esyn, double* delayed_synapse_array);
	double compute_one_psp(const double& si, Neuron* other_neuron);
	double compute_one_psp(const double& si, const double& esyn);
  double compute_one_psp(const Neuron& other_neuron);
  double compute_one_psp(double, string);
	void get_synapse_derivatives();
  void copy_current_synapse_wrapper(const Neuron& right_side);
  void copy_synapse_equations(const Neuron& right_side);
  void copy_current_equations(const Neuron& right_side);
  void copy_space_current_balance_results(const Neuron& right_side);
  void copy_space_synapse_results(const Neuron& right_side);
	void copy_synapse_array(const Neuron& right_side);
	void synapse_push();





};
//#endif
