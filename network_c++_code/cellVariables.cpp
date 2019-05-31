#include "cellVariables.h"
#include "mapFunctions.h"
#include <iostream>

//#include <map>
//#include <string>
//using namespace std;

// Default constructor

cellVariables::cellVariables(){
//initialize cell variables
	return;
 }

// Constructor
cellVariables::cellVariables(Fmap var, string cell_type){
 //initialize cell variables
	this->class_init(var, cell_type);
	return;
}

void cellVariables::class_init(Fmap& var, string cell_type){
	//cout<<"Class init cellVariables"<<endl;
 for (auto it = var.begin(); it != var.end(); ++it){
		this->cell_vars[ it->first ] = it->second;
	}
	this->cell_type = cell_type;
	cout<<"cell type is "<<cell_type<<endl;
	//cout<<"jfkasjkfsd"<<endl;
	membrane_potential = cell_vars["membrane_potential"];
	//&membrane_potential = temp;
	synapse_e = cell_vars["synapse_e"];
	synapse_e2 = cell_vars["synapse_e2"];
	synapse_i = cell_vars["synapse_i"];
	synapse_i_delayed = cell_vars["synapse_i_delayed"];
	synapse_i_max_current = cell_vars["synapse_i_max_current"];

	if (cell_type=="Cressman" || cell_type=="WB"){
		init_cressman(var);
	}
#if ANDERSON	
	else if (cell_type=="Anderson"){
		init_anderson(var);
	}
	else if (cell_type=="Anderson_sans_Ca"){
		init_anderson_sans_Ca(var);
	}
#endif
	else if (cell_type=="Froh"){
		//cout<<"Just before init_frohlich"<<endl;
		init_frohlich(var);
	}

	else{
		init_generic(var);
	}
	return;
}

//Assignment operator overloading
cellVariables& cellVariables::operator=(const cellVariables& right_side){
	if (this==&right_side){
		return *this;
	}
	else{
		//cout<<"cellV init"<<endl;
		Fmap v = right_side.cell_vars;
		string t = right_side.cell_type;
		this->class_init(v,t);
		//cout<<"before return; this is "<<this->membrane_potential<<endl;
		return *this;
	}
}

double cellVariables::get_values(string name){
  return cell_vars[name];
}

void cellVariables::init_generic(Fmap& v){
	potassium_channel = v["potassium_channel"];
	sodium_channel = v["sodium_channel"];
	return;
}

#if ANDERSON
void cellVariables::init_anderson_sans_Ca(Fmap& v){
	W = v["W"];
	X = v["X"];
	B = v["B"];
	}

void cellVariables::init_anderson(Fmap& v){
	init_anderson_sans_Ca(v);
	Ca_U0 = v["Ca_U0"];
	Ca_U1 = v["Ca_U1"];
	Ca_U2 = v["Ca_U2"];
	Ca_U3 = v["Ca_U3"];
	Ca_U4 = v["Ca_U4"];
	Ca_U5 = v["Ca_U5"];
	Ca_U6 = v["Ca_U6"];
	Ca_C0 = v["Ca_C0"];
	Ca_C1 = v["Ca_C1"];
	Ca_C2 = v["Ca_C2"];
	Ca_C3 = v["Ca_C3"];
	Ca_C4 = v["Ca_C4"];
	Ca_C5 = v["Ca_C5"];
	Ca_C6 = v["Ca_C6"];
	ICa = v["ICa"];
	ICa_next = v["ICa_next"];
	return;
}
#endif

void cellVariables::init_cressman(Fmap& v){
	potassium_channel = v["potassium_channel"];
	sodium_channel = v["sodium_channel"];
	Na_extracellular = v["Na_extracellular"];
	Na_intracellular = v["Na_intracellular"];
	K_extracellular = v["K_extracellular"];
	K_intracellular = v["K_intracellular"];
	//Cl_extracellular = v["Cl_extracellular"];
	//Cl_intracellular = v["Cl_intracellular"];
	Ca_intracellular = v["Ca_intracellular"];
	VNa_var = v["VNa_var"];
	VK_var = v["VK_var"];
	//VCl_var = v["VCl_var"];
	VL_var = v["VL_var"];
	potassium_current = v["potassium_current"];
	sodium_current = v["sodium_current"];
	K_extracellular_diffusion_term = v["K_extracellular_diffusion_term"];
	return;
}

void cellVariables::init_frohlich(Fmap &v){
	init_cressman(v);  // Most of the variables are not used in Frohlich
	//cout<<"K_extracellular "<<K_extracellular<<endl;
	frohlich_buffer=v["frohlich_buffer"];
	tau_now_Ca=v["tau_now_Ca"];
	return;
}








// Better to have a template for that, hopefully fix it later.
Fmap cellVariables::get_hash(){
  return cell_vars;
}
