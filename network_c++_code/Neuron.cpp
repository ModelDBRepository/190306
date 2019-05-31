#include "Neuron.h"
#include "cellVariables.h"
#include "cellParameters.h"

//#include <omp.h>
#include "currentEquations.h"
#include "synapseEquations.h"
#include "integration_routines.h"
#include <iostream>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "cellVariables.h"
//#include "cellParameters.h"
//#include <map>
//#include <string>
//using namespace std;

//Default Constructor
Neuron::Neuron():cellVariables(),cellParameters()
{
  //type of neuron
  type_of_neuron="0";
  type_of_synapse="0";
  inh_or_exc="0";
  
  
  //get current equations
  number_of_current_balance_equations=0;
  current_balance_equations=NULL;
  
  //get synapse equations
  number_of_synapse_equations=0;
  synapse_equations=NULL;

  //get space of integration and derivative results
  synapse_readings=NULL;
  synapse_derivatives=NULL;
  current_balance_readings=NULL;
  current_balance_derivatives=NULL;


  //set integration neuron flag
  update_neuron_happened=0;

	// delay
	synapse_array=NULL;

}


//Constructor
Neuron::Neuron(string type_neuron,string type_synapse, string i_or_e, Fmap vars, Fmap pars) : cellVariables(vars, type_neuron), cellParameters(pars, type_neuron, type_synapse){

  //Set all flags to zero
  got_current_balance_equations=0;
  got_synapse_equations=0;
  got_space_current_balance_results=0;
  got_space_synapse_results=0;
  update_neuron_happened=0;
	synapse_array_assigned=0;
	spiked = 0;

  //cout<<"INside constructor"<<endl;
  //cout<<" got_space_synapse_results is "<<got_space_synapse_results<<endl;
  //cout<<" get_space_current_balance_results is "<<got_space_current_balance_results<<endl;
  //type of neuron
  type_of_neuron=type_neuron;
  type_of_synapse=type_synapse;
  inh_or_exc=i_or_e;

  //get current equations
  get_current_balance_equations();
  get_space_current_balance_results();

  //get synapse equations
  get_synapse_equations();
  get_space_synapse_results();

	//assign synapse array
	assign_synapse_array();

}


// Copy constructor
Neuron::Neuron(const Neuron& right_side) : cellVariables(right_side.cell_vars, right_side.type_of_neuron), cellParameters(right_side.cell_pars, right_side.type_of_neuron, right_side.type_of_synapse){
  //cout<<"Inside copy constructor"<<endl;
   //Set all flags to zero
  got_current_balance_equations=0;
  got_synapse_equations=0;
  got_space_current_balance_results=0;
  got_space_synapse_results=0;
  update_neuron_happened=0;
	synapse_array_assigned=0;
	spiked = 0;

  type_of_neuron = right_side.type_of_neuron;
  type_of_synapse = right_side.type_of_synapse;
  inh_or_exc = right_side.inh_or_exc;
  copy_current_synapse_wrapper(right_side);

  get_space_current_balance_results();
  get_space_synapse_results();
  copy_space_current_balance_results(right_side);
  copy_space_synapse_results(right_side);
  assign_synapse_array();
}


//Assignment operator overloading
Neuron& Neuron::operator=(const Neuron& right_side){
  if (this == &right_side){
    //cout<<"Assignment operator same address"<<endl;
    return *this;
  }
  else{
    //cout<<"Assignment operator different address"<<endl;
    // Base class assignments
		//cout<<"Base class init"<<endl;
    cellVariables::operator=(right_side);
		//cout<<"After cellVariable init"<<endl;
    cellParameters::operator=(right_side);
    //Free all dynamic arrays on the left side
		//cout<<"Passed bass class init"<<endl;
    delete [] current_balance_equations;
    delete [] synapse_equations;
    delete [] current_balance_readings;
    delete [] synapse_readings;
    delete [] current_balance_derivatives;
    delete [] synapse_derivatives;
		delete [] synapse_array;
    // Copy functors from the right side
    type_of_neuron=right_side.type_of_neuron;
    type_of_synapse=right_side.type_of_synapse;
		spiked=right_side.spiked;
    inh_or_exc=right_side.inh_or_exc;
    got_current_balance_equations=0; // reset flags to zero
    got_synapse_equations=0;
    (*this).copy_current_synapse_wrapper(right_side);
    got_current_balance_equations=1; // Done copying, reset flags to one
    got_synapse_equations=1;
    // Copy reading and derivative information from the right side
    got_space_current_balance_results=0; //reset flags to zero
    got_space_synapse_results=0;
    (*this).get_space_current_balance_results();
    (*this).copy_space_current_balance_results(right_side);
    (*this).get_space_synapse_results();
    (*this).copy_space_synapse_results(right_side);
    //Set integration neuron flag
    update_neuron_happened=0;

		//Copy synapse array
		synapse_array_assigned = 0;
		(*this).assign_synapse_array();
		(*this).copy_synapse_array(right_side);
    return *this;
  }
}


//Destructor
Neuron::~Neuron()
{
  //cout<<"Destructor: begin freeing memories  current_balance_equations "<<endl;
  delete [] current_balance_equations;
  //cout<<"Destructor: begin freeing memories synapse_equations"<<endl;
  delete [] synapse_equations;
  //cout<<"Destructor: begin freeing memories current_balance_readings"<<endl;
  delete [] current_balance_readings;
  
  //cout<<"Destructorr: begin freeing memories  current_balance_derivatives"<<endl;
  delete [] current_balance_derivatives;
  //cout<<"Destructor: begin freeing memories synapse_derivatives"<<endl;
  delete [] synapse_derivatives;
  //cout<<"Destructor: begin freeing memories synapse_readings "<<endl;
  delete [] synapse_readings;
  //cout<<"End freeing memories"<<endl;
	delete [] synapse_array;
}
    


//-----------------------------------------------------------------------------------//
// copy current equations (from the right side to the left side) MUST BE PRIVATE
void Neuron::copy_current_equations(const Neuron& right_side){
  current_balance_equations = new DFuncPtr[number_of_current_balance_equations];
  assert(current_balance_equations!=NULL);
  for (int i=0; i<number_of_current_balance_equations; i++){
    current_balance_equations[i]=right_side.current_balance_equations[i];
  }
  return;
}

// copy synapse equations (from the right side to the left side) MUST BE PRIVATE
void Neuron::copy_synapse_equations(const Neuron& right_side){
  synapse_equations = new DFuncPtr[number_of_synapse_equations];
  assert(synapse_equations!=NULL);
  for (int i=0; i<number_of_synapse_equations; i++){
    synapse_equations[i]=right_side.synapse_equations[i];
    }
  return;
}

// Wrapper function for copying both current and synapse functions from the right to the left MUST BE PRIVATE
void Neuron::copy_current_synapse_wrapper(const Neuron& right_side){
  number_of_current_balance_equations = right_side.number_of_current_balance_equations;
  if (right_side.current_balance_equations!=NULL){
    copy_current_equations(right_side);
    
  }
  else{
    current_balance_equations=NULL;
  }
  number_of_synapse_equations = right_side.number_of_synapse_equations;
  if (right_side.synapse_equations!=NULL){
    copy_synapse_equations(right_side);
    
  }
  else{
    synapse_equations=NULL;
  }
  return;
}

void Neuron::copy_space_current_balance_results(const Neuron& right_side){
  assert(got_space_current_balance_results==1);
  for (int i=0; i<number_of_current_balance_equations; i++){
    *(current_balance_readings+i)=right_side.current_balance_readings[i];
    *(current_balance_derivatives+i)=right_side.current_balance_derivatives[i];
  }
  return;
}


void Neuron::copy_space_synapse_results(const Neuron& right_side){
 assert(got_space_synapse_results==1);
 for (int i=0; i<number_of_synapse_equations; i++){
   *(synapse_readings+i)=right_side.synapse_readings[i];
   *(synapse_derivatives+i)=right_side.synapse_derivatives[i];
 }
 return;
}

void Neuron::copy_synapse_array(const Neuron& right_side){
	assert(synapse_array_assigned==1);
	for (size_t i=0; i<synapse_array_length; i++){
		synapse_array[i]=right_side.synapse_array[i];
	}
	return;
}
	


//------------------------------------------------------------------------------------//


//------------------------------------------------------------------------------------//
// MUST BE PRIVATE
void Neuron::get_current_balance_equations(){
  // Do only once
  if (got_current_balance_equations==1){return;}
  got_current_balance_equations=1;
  
  // Assign appropriate equations
	if (type_of_neuron=="KC"){
    	number_of_current_balance_equations=2; 
    	current_balance_equations=new DFuncPtr[number_of_current_balance_equations];
    	current_balance_equations[0]=&KC_current_balance;
    	current_balance_equations[1]=&KC_potassium;
	}
	else if (type_of_neuron=="WB"){
    	number_of_current_balance_equations=6; 
    	current_balance_equations=new DFuncPtr[number_of_current_balance_equations];
    	current_balance_equations[0]=&wb_current_balance;
    	current_balance_equations[1]=&wb_potassium;
    	current_balance_equations[2]=&wb_sodium;
			current_balance_equations[3]=&cressman_Ca_total_derivative;
			current_balance_equations[4]=&cressman_K_total_derivative;
			current_balance_equations[5]=&cressman_Na_total_derivative;

	}
	else if (type_of_neuron=="Froh"){
			number_of_current_balance_equations=6; 
    	current_balance_equations=new DFuncPtr[number_of_current_balance_equations];
    	current_balance_equations[0]=&wb_current_balance;
    	current_balance_equations[1]=&wb_potassium;
    	current_balance_equations[2]=&wb_sodium;
			current_balance_equations[3]=&cressman_Ca_total_derivative;
			current_balance_equations[4]=&frohlich_K_total_derivative;
			current_balance_equations[5]=&frohlich_b_equation;
	}
#if ANDERSON
	else if (type_of_neuron=="Anderson"){
	 		number_of_current_balance_equations=18; 
			current_balance_equations=new DFuncPtr[number_of_current_balance_equations];
			current_balance_equations[0]=&anderson_current_balance;
			current_balance_equations[1]=&anderson_w_decay;
			current_balance_equations[2]=&anderson_b_decay;
			current_balance_equations[3]=&anderson_x_decay;
			current_balance_equations[4]=&anderson_Ca_U0_total_derivative;
			current_balance_equations[5]=&anderson_Ca_U1_total_derivative;
			current_balance_equations[6]=&anderson_Ca_U2_total_derivative;
			current_balance_equations[7]=&anderson_Ca_U3_total_derivative;
			current_balance_equations[8]=&anderson_Ca_U4_total_derivative;
			current_balance_equations[9]=&anderson_Ca_U5_total_derivative;
			current_balance_equations[10]=&anderson_Ca_U6_total_derivative;
			current_balance_equations[11]=&anderson_Ca_C0_total_derivative;
			current_balance_equations[12]=&anderson_Ca_C1_total_derivative;
			current_balance_equations[13]=&anderson_Ca_C2_total_derivative;
			current_balance_equations[14]=&anderson_Ca_C3_total_derivative;
			current_balance_equations[15]=&anderson_Ca_C4_total_derivative;
			current_balance_equations[16]=&anderson_Ca_C5_total_derivative;
			current_balance_equations[17]=&anderson_Ca_C6_total_derivative;	
	}
	else if (type_of_neuron=="Anderson_sans_Ca"){
	 		number_of_current_balance_equations=4; 
			current_balance_equations=new DFuncPtr[number_of_current_balance_equations];
			current_balance_equations[0]=&anderson_sans_Ca_current_balance;
			current_balance_equations[1]=&anderson_w_decay;
			current_balance_equations[2]=&anderson_b_decay;
			current_balance_equations[3]=&anderson_x_decay;
	}
#endif
	else if (type_of_neuron=="Cressman"){
			number_of_current_balance_equations=6;
			current_balance_equations=new DFuncPtr[number_of_current_balance_equations];
			current_balance_equations[0]=&cressman_current_balance;
			current_balance_equations[1]=&cressman_potassium;
			current_balance_equations[2]=&cressman_sodium;
			current_balance_equations[3]=&cressman_Ca_total_derivative;
			current_balance_equations[4]=&cressman_K_total_derivative;
			current_balance_equations[5]=&cressman_Na_total_derivative;
	}

	else if (type_of_neuron=="ABCsimp"){
    	number_of_current_balance_equations=2; 
    	current_balance_equations=new DFuncPtr[number_of_current_balance_equations];
    	current_balance_equations[0]=&abc_simple;
    	current_balance_equations[1]=&abc_simple_rc;
    	//double*(map<char,double>,map<char,double>) temp[]={&KC_crrent_balance,&KC_potassium_gating};
	}
	else if (type_of_neuron=="Iz"){
			number_of_current_balance_equations=2;
			current_balance_equations=new DFuncPtr[number_of_current_balance_equations];
			current_balance_equations[0]=&Iz_current_balance;
    	current_balance_equations[1]=&Iz_u;
	}
	else if (type_of_neuron=="Ahmed_modified"){
			number_of_current_balance_equations=2;
			current_balance_equations=new DFuncPtr[number_of_current_balance_equations];
			current_balance_equations[0]=&ahmed_modified_current_balance;
    	current_balance_equations[1]=&ahmed_modified_u;
	}
	else{		
    	number_of_current_balance_equations=0;
    	current_balance_equations=NULL;
	}
  // Put functors into appropriate spaces
  //current_balance_equations=new double*(map<char,double>,map<char,double>)[number_of_current_balance_equations];
  assert(current_balance_equations!=NULL);
  //for (size_t i=0; i<number_of_current_balance_equations; i++){
  //  *(current_balance_equations+i)=temp[i];
  //}
  return;
}

//PRIVATE
void Neuron::get_space_current_balance_results(){
  // Create space for derivatives and integrated results
  if (got_space_current_balance_results==1){return;}
  got_space_current_balance_results=1;
  current_balance_readings=new double[number_of_current_balance_equations]; // integrated results
  assert(current_balance_readings!=NULL);
  current_balance_derivatives=new double[number_of_current_balance_equations]; // derivatives
  assert(current_balance_derivatives!=NULL);
}


// MUST BE PRIVATE
void Neuron::get_synapse_equations(){
  // Do only once
  if (got_synapse_equations==1){return;}
  got_synapse_equations=1;
  
  if (type_of_synapse=="Normal"){
    number_of_synapse_equations=1;
    synapse_equations=new DFuncPtr[number_of_synapse_equations];
    synapse_equations[0]=&synapse_gating;
    //double*(map<char,double>,map<char,double>) temp[]={&synaptic_gating};
  }
  else if (type_of_synapse=="Normaldiscont"){
    number_of_synapse_equations=1;
    synapse_equations=new DFuncPtr[number_of_synapse_equations];
    synapse_equations[0]=&synapse_gating_decay_part;
  }
  else if (type_of_synapse=="VmD"){
    number_of_synapse_equations=2;
    synapse_equations=new DFuncPtr[number_of_synapse_equations];
    synapse_equations[0]=&synapse_gating;
    synapse_equations[1]=&vmd_gating;
    //double*(map<char,double>,map<char,double>) temp[]={&synaptic_gating, &vmd_gating};
  }
   else if (type_of_synapse=="VmDdiscont"){
    number_of_synapse_equations=2;
    synapse_equations=new DFuncPtr[number_of_synapse_equations];
    synapse_equations[0]=&synapse_gating_decay_part;
    synapse_equations[1]=&vmd_gating;
  }
	 else if (type_of_synapse=="VmD2discont"){
    number_of_synapse_equations=3;
    synapse_equations=new DFuncPtr[number_of_synapse_equations];
    synapse_equations[0]=&synapse_gating_decay_part;
    synapse_equations[1]=&vmd_gating; //ge
		 synapse_equations[2]=&vmd_gating2; //gi
  }
	else if (type_of_synapse=="Poissondiscont"){
		number_of_synapse_equations=2;
		synapse_equations=new DFuncPtr[number_of_synapse_equations];
    synapse_equations[0]=&synapse_gating_decay_part;
    synapse_equations[1]=&vmd_gating;
	}
  else{
    number_of_synapse_equations=0;
    synapse_equations=NULL;
    return;
  }
  //synapse_equations=new double*[number_of_synapse_equations](map<char,double>,map<char,double>); 
  assert(synapse_equations!=NULL);
  // for (size_t i=0; i<number_of_synapse_equations; i++){
  // *(synapse_equations+i)=temp[i];
  //}
  return;
}
void Neuron::synapse_push(){
	assert(synapse_array_assigned==1);
	double out = *(synapse_array+0);
	for (size_t i=0; i<(synapse_array_length-1); i++){
		*(synapse_array+i) = *(synapse_array+i+1);
	}
	*(synapse_array+synapse_array_length-1) = cellVariables::synapse_i;
	cellVariables::synapse_i_delayed = out;

	//double out = *(synapse_array_current_position);
	//cellVariables::["synapse_i_delayed"] = out;
	//*(synapse_array_current_position) = cellVariables::["synapse_i"];
	//int mem_diff = synapse_array_current_position - (synapse_array + 0);
	//cout<<"mem_diff is "<<mem_diff<<endl;
	//if (mem_diff < (synapse_array_length-1)*sizeof(double)){
	//if (mem_diff < (synapse_array_length-1)){
	//	synapse_array_current_position++;
	//}
	//else{
	//	synapse_array_current_position = synapse_array + 0;
	//}
	return;
}
void Neuron::assign_synapse_array(){
	if (synapse_array_assigned==1){return;}
	double delay = cellParameters::synapse_delay;
	double time_step = cellParameters::time_step;
	//cout<<"time step is "<<time_step<<endl;
	//cout<<"synapse delay is "<<delay<<endl;
	int num_slots = int(delay/time_step);
	if (num_slots<1){num_slots=1;}
	//cout<<"num slots is "<<num_slots<<endl;
	synapse_array = new double[num_slots];
	synapse_array_length = num_slots;
	for (size_t i=0; i<synapse_array_length; i++){
		synapse_array[i]=0;
	}
	synapse_array_assigned=1;
	synapse_array_current_position = synapse_array + 0;
	return;
}
//PRIVATE
void Neuron::get_space_synapse_results(){
  // Create space for derivatives and integrated results
  if (got_space_synapse_results==1){return;}
  got_space_synapse_results=1;
  //cout<<"In :get_space_synapse_results(), before new [number_of_synapse_equations is "<<number_of_synapse_equations<<endl;
  synapse_readings=new double[number_of_synapse_equations]; // integrated results
  assert(synapse_readings!=NULL);
  synapse_derivatives=new double[number_of_synapse_equations]; // derivatives
  assert(synapse_derivatives!=NULL);
  //cout<<"In :get_space_synapse_results(), before return"<<endl;
}

//----------------------------------------------------------------------------------------------//

void Neuron::get_current_balance_derivatives_no_psps(){
  for (size_t i=0; i<number_of_current_balance_equations; i++){
    //#pragma omp critical
    *(current_balance_derivatives+i)=(**(current_balance_equations+i))(*this, *this);
    //cout<<"In get_current_balance_derivative_no_psps: i is "<<i<<" value is "<< *(current_balance_derivatives+i)<<endl;
  }
  return;
}

double Neuron::get_summed_psps(double* synaptic_matrix_row, int number_of_neurons, Neuron* other_neurons){
  if (synaptic_matrix_row==NULL){return 0;}
  double summed_psps = 0;
  size_t i;
  //double one_psp;

   for (i=0; i<number_of_neurons; i++){
     if (synaptic_matrix_row[i]<=0) continue;
     //one_psp=compute_one_psp(*(other_neurons+i));
     //if (synaptic_matrix_row[i]<0) continue;
     //summed_psps+=(synaptic_matrix_row[i])*one_psp;
			summed_psps+=(synaptic_matrix_row[i])*compute_one_psp(*(other_neurons+i));
		//cout<<"summed psps is "<<summed_psps<<endl;
   }
 return summed_psps;
}

double Neuron::get_summed_psps(double* synaptic_matrix_row, int number_of_neurons, Neuron* other_neurons, double* delayed_synapse_array){
  if (synaptic_matrix_row==NULL){return 0;}
  double summed_psps = 0;
  size_t i;
  //double one_psp;

   for (i=0; i<number_of_neurons; i++){
     if (synaptic_matrix_row[i]<=0) continue;
     //one_psp=compute_one_psp(*(other_neurons+i));
     //if (synaptic_matrix_row[i]<0) continue;
     //summed_psps+=(synaptic_matrix_row[i])*one_psp;
			summed_psps+=(synaptic_matrix_row[i])*compute_one_psp(*(delayed_synapse_array+i), (other_neurons+i));
		//cout<<"summed psps is "<<summed_psps<<endl;
   }
 return summed_psps;
}

double Neuron::get_summed_psps(double* synaptic_matrix_row, int number_of_neurons, Neuron* other_neurons, double* delayed_synapse_array, int* valid_entries, int number_of_entries){//more efficient for looping
  if (synaptic_matrix_row==NULL){return 0;}
  double summed_psps = 0;
  size_t i;
  //double one_psp;

   for (i=0; i<number_of_entries; i++){
		int j = valid_entries[i];
		//if (j==-1)break;
		summed_psps+=synaptic_matrix_row[j]*compute_one_psp(*(delayed_synapse_array+j), other_neurons+j);
	
	}
   return summed_psps;
}
double Neuron::get_summed_psps(double* synaptic_matrix_row, int number_of_neurons, double* array_of_esyn, double* delayed_synapse_array, int* valid_entries, int number_of_entries){//more efficient for looping (both esyn and syn entries are prepped)
  if (synaptic_matrix_row==NULL){return 0;}
  double summed_psps = 0;
  size_t i;
  //double one_psp;

   for (i=0; i<number_of_entries; i++){
		int j = valid_entries[i];
		//if (j==-1)break;
		summed_psps+=synaptic_matrix_row[j]*compute_one_psp(*(delayed_synapse_array+j), array_of_esyn[j]);
	
	}
   return summed_psps;
}


double Neuron::get_summed_psps(double* synaptic_sorted_row, int* synaptic_location, int synaptic_size, double* array_of_esyn, double* delayed_synapse_array){//EVEN more efficient for looping (both esyn and syn entries are prepped)
  if (synaptic_sorted_row==NULL){return 0;}
  double summed_psps = 0;
  size_t i;
	assert(synaptic_size>=0);
  //double one_psp;

   for (i=0; i<synaptic_size; i++){
		int j = synaptic_location[i];
		//if (j==-1)break;
		summed_psps+=synaptic_sorted_row[i]*compute_one_psp(*(delayed_synapse_array+j), array_of_esyn[j]);	
	}
   return summed_psps;
}

int* Neuron::find_valid_synapse_entries(double* synapse_matrix_row, int number_of_neurons, int& number_of_entries){
		int* valid = new int[number_of_neurons];
		for (int i=0; i<number_of_neurons; i++){
			valid[i] = -1;
		}
		int count = 0;
		for (int i=0; i<number_of_neurons; i++){
			if (*(synapse_matrix_row+i)>1e-50){
				valid[count] = i;
				count++;
			}
		}
	number_of_entries = count;
		return valid;
}

double Neuron::prep_synapse_reversal_potential(){
	if(this->inh_or_exc=="exc"){cellParameters::synapse_potential=cellParameters::reversal_potential_e;}
  else if(this->inh_or_exc=="inh"){cellParameters::synapse_potential=cellParameters::reversal_potential_i;}
  else{cout<<"In Neuron::prep_synapse_reversal_potential: Unknown inh or exc"<<endl; exit (1);}
	return cellParameters::synapse_potential;
}


// This is for externally supplied synaptic variables
double Neuron::compute_one_psp(const double& si, Neuron* other_neuron){
	//double rpi =other_neuron.cellParameters::synapse_potential;
	//if(other_neuron->inh_or_exc=="exc"){rpi=cellParameters::reversal_potential_e;}
  //else if(other_neuron->inh_or_exc=="inh"){rpi=cellParameters::reversal_potential_i;}
  //else{cout<<"Unknown inh or exc"<<endl; exit;}
//	double mp=cellVariables::membrane_potential;
  //double rpi=cellParameters::["reversal_potential_i"];
  return (other_neuron->cellParameters::synapse_potential-cellVariables::membrane_potential)*si;
}
// This is for externally supplied synaptic variables
double Neuron::compute_one_psp(const double& si, const double& esyn){
	  return (esyn-cellVariables::membrane_potential)*si;
}



double Neuron::compute_one_psp(double si, string i_or_e){
  double rpi;
  if(i_or_e=="exc"){rpi=cellParameters::reversal_potential_e;}
  else if(i_or_e=="inh"){rpi=cellParameters::reversal_potential_i;}
  else{cout<<"Unknown inh or exc"<<endl; exit (1);}
  double mp=cellVariables::membrane_potential;
  //double rpi=cellParameters::["reversal_potential_i"];
  return (rpi-mp)*si;
}

double Neuron::compute_one_psp(const Neuron& other_neuron){
  double mp=cellVariables::membrane_potential;
  double rpi;
  //cellVariables other_neuron_cellVariables::=&(other_neuron.cellVariables);
  double si=other_neuron.cellVariables::synapse_i_delayed;
  if(other_neuron.inh_or_exc=="exc"){rpi=cellParameters::reversal_potential_e;}
  else if(other_neuron.inh_or_exc=="inh"){rpi=cellParameters::reversal_potential_i;}
  else{cout<<"Unknown inh or exc"<<endl; exit;}
  //cout<<"mp is "<<mp<<endl;
  //cout<<"rpi is "<<rpi<<endl;
  //cout<<"si is "<<si<<endl;
  return (rpi-mp)*si;
}

void Neuron::get_synapse_derivatives(){
  for (size_t i=0; i<number_of_synapse_equations; i++){
    *(synapse_derivatives+i)=(**(synapse_equations+i))(*this,*this);
    //cout<<"(**(synapse_equations+i))(cellVariables::,cellParameters::) is "<<(**(synapse_equations+i))(cellVariables::,cellParameters::)<<endl;
  }
  return;
}

//PUBLIC BELOW
void Neuron::get_all_derivatives(double* synaptic_matrix_row, int number_of_neurons, Neuron* other_neurons){
  check_and_get_spaces_all_derivatives();
  double summed_psps=get_summed_psps(synaptic_matrix_row, number_of_neurons, other_neurons);
  assign_all_derivatives_to_neuron(summed_psps);
}

//This version of get_all_derivatives assumes all inhibitory connection
void Neuron::get_all_derivatives(double si,double syn_strength){
  check_and_get_spaces_all_derivatives();
  double summed_psps=syn_strength*compute_one_psp(si, "inh");
  assign_all_derivatives_to_neuron(summed_psps);
}

//Externally supplied synapse variables
void Neuron::get_all_derivatives(double* synaptic_matrix_row, int number_of_neurons, Neuron* other_neurons, double* sv){
  check_and_get_spaces_all_derivatives();
  double summed_psps=get_summed_psps(synaptic_matrix_row, number_of_neurons, other_neurons, sv);
  assign_all_derivatives_to_neuron(summed_psps);
}

//Externally supplied synapse variables + synapse prepped for efficient looping
void Neuron::get_all_derivatives(double* synaptic_matrix_row, int number_of_neurons, Neuron* other_neurons, double* sv, int* valid_entries, int number_of_entries){
  check_and_get_spaces_all_derivatives();
  double summed_psps=get_summed_psps(synaptic_matrix_row, number_of_neurons, other_neurons, sv, valid_entries, number_of_entries);
  assign_all_derivatives_to_neuron(summed_psps);
}
//Externally supplied synapse variables + synapse AND esyn prepped for efficient looping
void Neuron::get_all_derivatives(double* synaptic_matrix_row, int number_of_neurons, double* prepped_esyn_array, double* sv, int* valid_entries, int number_of_entries){
  check_and_get_spaces_all_derivatives();
  double summed_psps=get_summed_psps(synaptic_matrix_row, number_of_neurons, prepped_esyn_array, sv, valid_entries, number_of_entries);
  assign_all_derivatives_to_neuron(summed_psps);
}
//Externally supplied synapse variable + synapse AND esyn prepped (even more) for efficient MPI looping
void Neuron::get_all_derivatives(double* synaptic_sorted_row, int* synaptic_location, int synaptic_size, double* prepped_esyn_array, double* sv){
  check_and_get_spaces_all_derivatives();
	double summed_psps = 0;
	if (synaptic_size>0) summed_psps=get_summed_psps(synaptic_sorted_row, synaptic_location, synaptic_size, prepped_esyn_array, sv);
  assign_all_derivatives_to_neuron(summed_psps);
	return;
}

void Neuron::check_and_get_spaces_all_derivatives(){
  assert(current_balance_derivatives!=NULL);
  assert(got_current_balance_equations==1);
  assert(got_synapse_equations==1);
  get_space_synapse_results();
  get_space_current_balance_results();
  return;
}

void Neuron::assign_all_derivatives_to_neuron(double summed_psps){
		get_current_balance_derivatives_no_psps();
		*(current_balance_derivatives+0)+=summed_psps;
		//cellParameters::D=summed_psps;
		if (type_of_synapse=="VmD" || type_of_synapse=="VmDdiscont" || type_of_synapse=="Poissondiscont" || type_of_synapse=="VmD2discont"){
			 *(current_balance_derivatives+0)+=cellVariables::synapse_e*(cellParameters::reversal_potential_e-cellVariables::membrane_potential); //ge
			 //cellParameters::aA=cellVariables::synapse_e*(cellParameters::reversal_potential_e-cellVariables::membrane_potential); 
			 //cout<<"synapse e is "<<cellVariables::synapse_e<<endl;
		}
		if (type_of_synapse=="VmD2discont"){
			*(current_balance_derivatives+0)+=cellVariables::synapse_e2*(cellParameters::reversal_potential_i-cellVariables::membrane_potential); //gi
			//cout<<"synapse e2 is "<<cellVariables::synapse_e2<<endl;
		}

	get_synapse_derivatives();
	return;
}



void Neuron::update_neuron(){
  assert(current_balance_readings!=NULL);
  assert(got_current_balance_equations==1);
  assert(got_space_current_balance_results==1);

	if (type_of_neuron=="KC"){
    cellVariables::membrane_potential=*(current_balance_readings+0);
    cellVariables::potassium_channel=*(current_balance_readings+1);
    cellVariables::sodium_channel=KC_sodium(cellVariables::membrane_potential);
  	}
	else if (type_of_neuron=="WB"){
    cellVariables::membrane_potential=*(current_balance_readings+0);
    cellVariables::potassium_channel=*(current_balance_readings+1);
    cellVariables::sodium_channel=*(current_balance_readings+2);
		cellVariables::Ca_intracellular=*(current_balance_readings+3);
		cellVariables::K_extracellular=*(current_balance_readings+4);
		cellVariables::Na_intracellular=*(current_balance_readings+5);

  	}
		else if (type_of_neuron=="Froh"){
    cellVariables::membrane_potential=*(current_balance_readings+0);
    cellVariables::potassium_channel=*(current_balance_readings+1);
    cellVariables::sodium_channel=*(current_balance_readings+2);
		cellVariables::Ca_intracellular=*(current_balance_readings+3);
		cellVariables::K_extracellular=*(current_balance_readings+4);
		cellVariables::frohlich_buffer=*(current_balance_readings+5);

  	}

#if ANDERSON	
	else if (type_of_neuron=="Anderson"){
		cellVariables::membrane_potential=*(current_balance_readings+0);
		cellVariables::W=*(current_balance_readings+1);
		cellVariables::B=*(current_balance_readings+2);
		cellVariables::X=*(current_balance_readings+3);
		cellVariables::Ca_U0=*(current_balance_readings+4);
		cellVariables::Ca_U1=*(current_balance_readings+5);
		cellVariables::Ca_U2=*(current_balance_readings+6);
		cellVariables::Ca_U3=*(current_balance_readings+7);
		cellVariables::Ca_U4=*(current_balance_readings+8);
		cellVariables::Ca_U5=*(current_balance_readings+9);
		cellVariables::Ca_U6=*(current_balance_readings+10);
		cellVariables::Ca_C0=*(current_balance_readings+11);
		cellVariables::Ca_C1=*(current_balance_readings+12);
		cellVariables::Ca_C2=*(current_balance_readings+13);
		cellVariables::Ca_C3=*(current_balance_readings+14);
		cellVariables::Ca_C4=*(current_balance_readings+15);
		cellVariables::Ca_C5=*(current_balance_readings+16);
		cellVariables::Ca_C6=*(current_balance_readings+17);
		cellVariables::ICa=cellVariables::ICa_next;
		}
	else if (type_of_neuron=="Anderson_sans_Ca"){
		cellVariables::membrane_potential=*(current_balance_readings+0);
		cellVariables::W=*(current_balance_readings+1);
		cellVariables::B=*(current_balance_readings+2);
		cellVariables::X=*(current_balance_readings+3);
		}
#endif
	else if (type_of_neuron=="Cressman"){
		cellVariables::membrane_potential=*(current_balance_readings+0);
		cellVariables::potassium_channel=*(current_balance_readings+1);
		cellVariables::sodium_channel=*(current_balance_readings+2);
		cellVariables::Ca_intracellular=*(current_balance_readings+3);
		cellVariables::K_extracellular=*(current_balance_readings+4);
		cellVariables::Na_intracellular=*(current_balance_readings+5);
		//HARD CODE the change of Egaba for the moment.  When it becomes more mature we will find better places to code it.
		//if (this->inh_or_exc=="inh"){
			//double tau1=1000*(5+395/(1+exp(-3.0*(*(current_balance_readings+4)-5.0))));// in ms #8.5 #5.5
			//double tau2=1000*(2+395/(1+exp(4.0*(*(current_balance_readings+4)-5.0))));// in ms #8.5 #5.5
			// cellParameters::reversal_potential_i+=-time_step*(cellParameters::reversal_potential_i+72.0)/tau1-time_step*(cellParameters::reversal_potential_i+66.5)/tau2;
			// cellParameters::reversal_potential_i+=-time_step*(cellParameters::reversal_potential_i+72.0)/tau1-time_step*(cellParameters::reversal_potential_i+58.0)/tau2; //para 3
			// cellParameters::reversal_potential_i+=-time_step*(cellParameters::reversal_potential_i+72.0)/tau1-time_step*(cellParameters::reversal_potential_i+20.0)/tau2; // para 4
			//cellParameters::reversal_potential_i+=-time_step*(cellParameters::reversal_potential_i+72.0)/tau1-time_step*(cellParameters::reversal_potential_i+50.0)/tau2; // para 5


		//	}
	}
		else if (type_of_neuron=="ABCsimp"){
			if (*(current_balance_readings+0) > cellParameters::Vmax){
				cellVariables::membrane_potential=cellParameters::Vreset;
      	cellVariables::potassium_channel=cellParameters::Rreset;
    	}
    	else{
      	cellVariables::membrane_potential=*(current_balance_readings+0);
      	cellVariables::potassium_channel=*(current_balance_readings+1);
    	}
		}
	else if (type_of_neuron=="Iz" || type_of_neuron=="Ahmed_modified"){
			if (*(current_balance_readings+0) > cellParameters::Vmax){
      	cellVariables::membrane_potential=cellParameters::Vreset;
      	cellVariables::potassium_channel=cellParameters::d_step*cellVariables::potassium_channel+cellParameters::Rreset;
				//cout<<"Rreset is "<<cellParameters::Rreset<<endl;
    	}
    	else{
      	cellVariables::membrane_potential=*(current_balance_readings+0);
      	cellVariables::potassium_channel=*(current_balance_readings+1);
    	}
	}
	return;
}


void Neuron::update_synapse(){
  assert(synapse_readings!=NULL);
  assert(got_synapse_equations==1);
  assert(got_space_synapse_results==1);
  //cout<<"*(synapse_readings+0) is "<<*(synapse_readings+0)<<endl;
	if (spiked==1){spiked = 0;}
  cellVariables::synapse_i=*(synapse_readings+0);
	double d_factor=1.0;
#ifdef __FIGTWO__
	// Fig 2 gamma/theta (linearly ramping up inhibition right after stimulation)

	if (inh_or_exc=="inh" && cellParameters::stim_end>cellVariables::current_time && cellVariables::current_time>cellParameters::stim_start){ 
					double delta_t = (cellParameters::stim_end-cellParameters::stim_start);
					double tc = delta_t/40.0;
					double terminal_g = cellParameters::inhib_ramp_factor*cellParameters::synapse_i_max;
					double factor_g = -(cellVariables::synapse_i_max_current - terminal_g)*cellParameters::time_step/(tc*cellVariables::synapse_i_max_current);
					d_factor=1.0+factor_g;
					cellVariables::synapse_i_max_current=cellVariables::synapse_i_max_current*d_factor;
//cout<<"d_factor is "<<d_factor<<endl;
	//		cout<<"synapse i is "<<cellVariables::synapse_i_max_current<<endl;

				}
#endif


	if (type_of_synapse=="Normaldiscont" || type_of_synapse=="VmDdiscont"  || type_of_synapse=="Poissondiscont" || type_of_synapse=="VmD2discont"){
		 	if (cellVariables::membrane_potential<cellParameters::Vopens && *(current_balance_readings+0)>=cellParameters::Vopens){
#ifdef __DEPLETE__				
			//	cout<<"Vopens is "<<cellParameters::Vopens<<endl;
      // hard reset to synapse_i_max after spike; if use this type of synapse, please update synapse first then neuron
			// HARD CODE SYNAPTIC DEPRESSION VARIABLE HERE...
			//	double d_factor = 1.0;
			if (inh_or_exc=="inh" && cellVariables::current_time>49500){
				//d_factor = 0.4;
				d_factor = 0.1;
				cellVariables::synapse_i_max_current=cellVariables::synapse_i_max_current*d_factor;
				}
#endif
			//if (cellVariables::synapse_i_max_current>=0.50){
			//	d_factor = cellParameters::rate_depression;
			//}
			//else if (cellVariables::synapse_i_max_current<0.50 && cellVariables::synapse_i_max_current>0.0001){ // original 0.005
				//	d_factor = 0.2;
			//}
			//else{
			//		d_factor = 1;
			//}
			//cellVariables::synapse_i_max_current*=cellParameters::rate_depression;
			//cellVariables::synapse_i_max_current=cellVariables::synapse_i_max_current*d_factor;
			cellVariables::synapse_i+=cellVariables::synapse_i_max_current;
			//cout<<"d_factor is "<<d_factor<<endl;
			//cout<<"synapse i is "<<cellVariables::synapse_i<<endl;

      //cellVariables::synapse_i+=cellParameters::synapse_i_max;
			spiked = 1;
			//cout<<"synapse i max is "<<cellParameters::synapse_i_max<<endl;
			//
    	}
			// HARD CODE SYNAPTIC RECOVERY VARIABLE HERE
			//cellVariables::synapse_i_max_current+=(cellParameters::synapse_i_max-cellVariables::synapse_i_max_current)*time_step/cellParameters::tau_recovery;
	}
	if (type_of_synapse=="VmD" || type_of_synapse=="VmDdiscont" || type_of_synapse=="Poissondiscont"  || type_of_synapse=="VmD2discont"){
			cellVariables::synapse_e=*(synapse_readings+1);
	} 
	if (type_of_synapse=="VmD2discont"){
			cellVariables::synapse_e2=*(synapse_readings+2);
	}
	synapse_push();  //for delay
	return;
}

void Neuron::update(){
	update_synapse();
	update_neuron();
	return;
}


void Neuron::integrate_synapse(const double& dt, const double& sqroot_dt, long *seed, int* iset_in, float* gset_in, int* iset_out, float* gset_out, long* idum2_in, long* iy_in, long* iv_in, long* idum2_out, long* iy_out, long* iv_out, int& done_init_gasdev_flag, int& done_init_ran2_flag, int& done_dump_gasdev_flag, int& done_dump_ran2_flag, const int& thread){
  assert(got_synapse_equations==1);
  assert(got_space_synapse_results==1);
  *(synapse_readings+0)=euler_next(cellVariables::synapse_i, *(synapse_derivatives+0), dt); // inhib
	if (type_of_synapse=="VmD" || type_of_synapse=="VmDdiscont" || type_of_synapse=="VmD2discont"){
			*(synapse_readings+1)=euler_next(cellVariables::synapse_e, *(synapse_derivatives+1), dt, sqroot_dt, cellParameters::vmd_noise, seed, iset_in, gset_in, iset_out, gset_out, idum2_in, iy_in, iv_in, idum2_out, iy_out, iv_out, done_init_gasdev_flag, done_init_ran2_flag, done_dump_gasdev_flag, done_dump_ran2_flag, thread); // ge
	}
	else if (type_of_synapse=="Poissondiscont"){
	*(synapse_readings+1)=poisson_next(cellVariables::synapse_e, *(synapse_derivatives+1), cellParameters::poisson_rate, cellParameters::poisson_increment, dt, sqroot_dt, seed, thread);
	}
	if (type_of_synapse=="VmD2discont"){
		*(synapse_readings+2)=euler_next(cellVariables::synapse_e2, *(synapse_derivatives+2), dt, sqroot_dt, cellParameters::vmd_noise2, seed, iset_in, gset_in, iset_out, gset_out, idum2_in, iy_in, iv_in, idum2_out, iy_out, iv_out, done_init_gasdev_flag, done_init_ran2_flag, done_dump_gasdev_flag, done_dump_ran2_flag, thread); // gi
	}
  return;
}
void Neuron::integrate_synapse(const double& dt, const double& sqroot_dt, long *seed, int* iset_in, float* gset_in, int* iset_out, float* gset_out, long* idum2_in, long* iy_in, long* iv_in, long* idum2_out, long* iy_out, long* iv_out, int& done_init_gasdev_flag, int& done_init_ran2_flag, int& done_dump_gasdev_flag, int& done_dump_ran2_flag){
	int thread = 0;
	integrate_synapse(dt, sqroot_dt, seed, iset_in, gset_in, iset_out, gset_out, idum2_in, iy_in, iv_in, idum2_out, iy_out, iv_out, done_init_gasdev_flag, done_init_ran2_flag, done_dump_gasdev_flag, done_dump_ran2_flag,thread);
	return;
}

void Neuron::integrate_synapse(const double& dt, const double& sqroot_dt, long *seed, const int& thread){
  int dummy = 1;
  integrate_synapse(dt, sqroot_dt, seed, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, dummy, dummy, dummy, dummy, thread);
  return;
} 
void Neuron::integrate_synapse(const double& dt, const double& sqroot_dt, long *seed){
  int dummy = 1;
	int thread = 0;
  integrate_synapse(dt, sqroot_dt, seed, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, dummy, dummy, dummy, dummy, thread);
  return;
} 


void Neuron::integrate_neuron(const double& dt, const double& sqroot_dt, long *seed, int* iset_in, float* gset_in, int* iset_out, float* gset_out, long* idum2_in, long* iy_in, long* iv_in, long* idum2_out, long* iy_out, long* iv_out, int& done_init_gasdev_flag, int& done_init_ran2_flag, int& done_dump_gasdev_flag, int& done_dump_ran2_flag, const int& thread){
  assert(got_current_balance_equations==1);
  assert(got_space_current_balance_results==1);
  //cout<<"Integrating neuron"<<endl;
  //*(current_balance_readings+0)=euler_next(cellVariables::["membrane_potential"], *(current_balance_derivatives+0), dt, sqroot_dt, cellParameters::["current_noise"], seed);
  //*(current_balance_readings+0)=euler_next(cellVariables::["membrane_potential"], *(current_balance_derivatives+0), dt); // test without noise first
	if (type_of_neuron=="KC" || type_of_neuron=="Iz" || type_of_neuron=="Ahmed_modified"){
		*(current_balance_readings+0)=euler_next(cellVariables::membrane_potential, *(current_balance_derivatives+0), dt, sqroot_dt, cellParameters::current_noise, seed, iset_in, gset_in, iset_out, gset_out, idum2_in, iy_in, iv_in, idum2_out, iy_out, iv_out, done_init_gasdev_flag, done_init_ran2_flag, done_dump_gasdev_flag, done_dump_ran2_flag, thread);
    *(current_balance_readings+1)=euler_next(cellVariables::potassium_channel, *(current_balance_derivatives+1), dt);
		//cout<<"Derivative 0 is "<<*(current_balance_derivatives+0)<<endl;
		//cout<<"Derivative 1 is "<<*(current_balance_derivatives+1)<<endl;
		//cout<<endl;

	}
	else if (type_of_neuron=="WB"){
    *(current_balance_readings+0)=euler_next(cellVariables::membrane_potential, *(current_balance_derivatives+0), dt, sqroot_dt, cellParameters::current_noise, seed, iset_in, gset_in, iset_out, gset_out, idum2_in, iy_in, iv_in, idum2_out, iy_out, iv_out, done_init_gasdev_flag, done_init_ran2_flag, done_dump_gasdev_flag, done_dump_ran2_flag, thread);
    *(current_balance_readings+1)=euler_next(cellVariables::potassium_channel, *(current_balance_derivatives+1), dt);
    *(current_balance_readings+2)=euler_next(cellVariables::sodium_channel, *(current_balance_derivatives+2), dt);
		*(current_balance_readings+3)=euler_next(cellVariables::Ca_intracellular, *(current_balance_derivatives+3), dt);
	*(current_balance_readings+4)=euler_next(cellVariables::K_extracellular, *(current_balance_derivatives+4), dt);
	*(current_balance_readings+5)=euler_next(cellVariables::Na_intracellular, *(current_balance_derivatives+5), dt);
	}
	else if (type_of_neuron=="Froh"){
*(current_balance_readings+0)=euler_next(cellVariables::membrane_potential, *(current_balance_derivatives+0), dt, sqroot_dt, cellParameters::current_noise, seed, iset_in, gset_in, iset_out, gset_out, idum2_in, iy_in, iv_in, idum2_out, iy_out, iv_out, done_init_gasdev_flag, done_init_ran2_flag, done_dump_gasdev_flag, done_dump_ran2_flag, thread, cellVariables::current_time, cellParameters::stim_start, cellParameters::stim_end, cellParameters::stim_strength);
		 //*(current_balance_readings+0)=euler_next(cellVariables::membrane_potential, *(current_balance_derivatives+0), dt, sqroot_dt, cellParameters::current_noise, seed, iset_in, gset_in, iset_out, gset_out, idum2_in, iy_in, iv_in, idum2_out, iy_out, iv_out, done_init_gasdev_flag, done_init_ran2_flag, done_dump_gasdev_flag, done_dump_ran2_flag, thread);

    	*(current_balance_readings+1)=euler_next(cellVariables::potassium_channel, *(current_balance_derivatives+1), dt);
    	*(current_balance_readings+2)=euler_next(cellVariables::sodium_channel, *(current_balance_derivatives+2), dt);
	*(current_balance_readings+3)=euler_next(cellVariables::Ca_intracellular, *(current_balance_derivatives+3), dt);
	*(current_balance_readings+4)=euler_next(cellVariables::K_extracellular, *(current_balance_derivatives+4), dt);
	*(current_balance_readings+5)=euler_next(cellVariables::frohlich_buffer, *(current_balance_derivatives+5), dt);
	}
	else if (type_of_neuron=="ABCsimp"){
    *(current_balance_readings+0)=euler_next(cellVariables::membrane_potential, *(current_balance_derivatives+0), dt, sqroot_dt, cellParameters::c*cellParameters::current_noise, seed,iset_in, gset_in, iset_out, gset_out, idum2_in, iy_in, iv_in, idum2_out, iy_out, iv_out, done_init_gasdev_flag, done_init_ran2_flag, done_dump_gasdev_flag, done_dump_ran2_flag, thread);
    *(current_balance_readings+1)=euler_next(cellVariables::potassium_channel, *(current_balance_derivatives+1), dt);
	}
#if ANDERSON
	else if (type_of_neuron=="Anderson"){
 		*(current_balance_readings+0)=euler_next(cellVariables::membrane_potential, *(current_balance_derivatives+0), dt, sqroot_dt, cellParameters::current_noise, seed, iset_in, gset_in, iset_out, gset_out, idum2_in, iy_in, iv_in, idum2_out, iy_out, iv_out, done_init_gasdev_flag, done_init_ran2_flag, done_dump_gasdev_flag, done_dump_ran2_flag, thread);
    *(current_balance_readings+1)=euler_next(cellVariables::W, *(current_balance_derivatives+1), dt);
		*(current_balance_readings+2)=euler_next(cellVariables::B, *(current_balance_derivatives+2), dt);
		*(current_balance_readings+3)=euler_next(cellVariables::X, *(current_balance_derivatives+3), dt);
		*(current_balance_readings+4)=euler_next(cellVariables::Ca_U0, *(current_balance_derivatives+4), dt);
		*(current_balance_readings+5)=euler_next(cellVariables::Ca_U1, *(current_balance_derivatives+5), dt);
		*(current_balance_readings+6)=euler_next(cellVariables::Ca_U2, *(current_balance_derivatives+6), dt);
		*(current_balance_readings+7)=euler_next(cellVariables::Ca_U3, *(current_balance_derivatives+7), dt);
		*(current_balance_readings+8)=euler_next(cellVariables::Ca_U4, *(current_balance_derivatives+8), dt);
		*(current_balance_readings+9)=euler_next(cellVariables::Ca_U5, *(current_balance_derivatives+9), dt);
		*(current_balance_readings+10)=euler_next(cellVariables::Ca_U6, *(current_balance_derivatives+10), dt);
		*(current_balance_readings+11)=euler_next(cellVariables::Ca_C0, *(current_balance_derivatives+11), dt);
		*(current_balance_readings+12)=euler_next(cellVariables::Ca_C1, *(current_balance_derivatives+12), dt);
		*(current_balance_readings+13)=euler_next(cellVariables::Ca_C2, *(current_balance_derivatives+13), dt);
		*(current_balance_readings+14)=euler_next(cellVariables::Ca_C3, *(current_balance_derivatives+14), dt);
		*(current_balance_readings+15)=euler_next(cellVariables::Ca_C4, *(current_balance_derivatives+15), dt);
		*(current_balance_readings+16)=euler_next(cellVariables::Ca_C5, *(current_balance_derivatives+16), dt);
		*(current_balance_readings+17)=euler_next(cellVariables::Ca_C6, *(current_balance_derivatives+17), dt);
	}
	else if (type_of_neuron=="Anderson_sans_Ca"){
		*(current_balance_readings+0)=euler_next(cellVariables::membrane_potential, *(current_balance_derivatives+0), dt, sqroot_dt, cellParameters::current_noise, seed, iset_in, gset_in, iset_out, gset_out, idum2_in, iy_in, iv_in, idum2_out, iy_out, iv_out, done_init_gasdev_flag, done_init_ran2_flag, done_dump_gasdev_flag, done_dump_ran2_flag, thread);
    *(current_balance_readings+1)=euler_next(cellVariables::W, *(current_balance_derivatives+1), dt);
		*(current_balance_readings+2)=euler_next(cellVariables::B, *(current_balance_derivatives+2), dt);
		*(current_balance_readings+3)=euler_next(cellVariables::X, *(current_balance_derivatives+3), dt);
	}
#endif
	else if (type_of_neuron=="Cressman"){
*(current_balance_readings+0)=euler_next(cellVariables::membrane_potential, *(current_balance_derivatives+0), dt, sqroot_dt, cellParameters::current_noise, seed, iset_in, gset_in, iset_out, gset_out, idum2_in, iy_in, iv_in, idum2_out, iy_out, iv_out, done_init_gasdev_flag, done_init_ran2_flag, done_dump_gasdev_flag, done_dump_ran2_flag, thread);
    *(current_balance_readings+1)=euler_next(cellVariables::potassium_channel, *(current_balance_derivatives+1), dt);
		*(current_balance_readings+2)=euler_next(cellVariables::sodium_channel, *(current_balance_derivatives+2), dt);
		*(current_balance_readings+3)=euler_next(cellVariables::Ca_intracellular, *(current_balance_derivatives+3), dt);
	*(current_balance_readings+4)=euler_next(cellVariables::K_extracellular, *(current_balance_derivatives+4), dt);
	*(current_balance_readings+5)=euler_next(cellVariables::Na_intracellular, *(current_balance_derivatives+5), dt);
	}
	//else if (type_of_neuron=="Ahmed_modified"){
	//*(current_balance_readings+0)=euler_next(cellVariables::membrane_potential, *(current_balance_derivatives+0), dt, sqroot_dt, cellParameters::current_noise, seed, iset_in, gset_in, iset_out, gset_out, idum2_in, iy_in, iv_in, idum2_out, iy_out, iv_out, done_init_gasdev_flag, done_init_ran2_flag, done_dump_gasdev_flag, done_dump_ran2_flag, thread);
  //*(current_balance_readings+1)=euler_next(cellVariables::potassium_channel, *(current_balance_derivatives+1), dt);		
	//}
  return;
}
void Neuron::integrate_neuron(const double& dt, const double& sqroot_dt, long *seed, int* iset_in, float* gset_in, int* iset_out, float* gset_out, long* idum2_in, long* iy_in, long* iv_in, long* idum2_out, long* iy_out, long* iv_out, int& done_init_gasdev_flag, int& done_init_ran2_flag, int& done_dump_gasdev_flag, int& done_dump_ran2_flag){
	int thread = 0;
	integrate_neuron(dt, sqroot_dt, seed, iset_in, gset_in, iset_out, gset_out, idum2_in, iy_in, iv_in, idum2_out, iy_out, iv_out, done_init_gasdev_flag, done_init_ran2_flag, done_dump_gasdev_flag, done_dump_ran2_flag, thread);
	return;
}

void Neuron::integrate_neuron(const double& dt, const double& sqroot_dt, long *seed, const int& thread){
  int dummy = 1;
  integrate_neuron(dt, sqroot_dt, seed, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, dummy, dummy, dummy, dummy, thread);
  return;
} 
void Neuron::integrate_neuron(const double& dt, const double& sqroot_dt, long *seed){
  int dummy = 1;
	int thread = 0;
  integrate_neuron(dt, sqroot_dt, seed, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, dummy, dummy, dummy, dummy, thread);
  return;
} 

//INPUT AND OUTPUT FUNCTIONS

//Added October 22nd, 2010
double Neuron::display_current_balance_derivative(int i){
  assert(current_balance_derivatives!=NULL);
  return (*(current_balance_derivatives+i));
}

string Neuron::display_type(){
	return type_of_neuron;
}
