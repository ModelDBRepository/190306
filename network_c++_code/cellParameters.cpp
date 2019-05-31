#include "cellParameters.h"
#include "mapFunctions.h"
//#include <map>
#include <string>
#include <iostream>
//using namespace std

// Default constructor

cellParameters::cellParameters(){
	return;	
}

// Constructor
cellParameters::cellParameters(Fmap p, string cell_type, string synapse_type){
	this->class_init(p, cell_type, synapse_type);
	return;
}

void cellParameters::class_init(Fmap& p, string cell_type, string synapse_type){
	for (auto it = p.begin(); it != p.end(); ++it){
		this->cell_pars[ it->first ] = it->second;
	}
	this->cell_type = cell_type;
	this->synapse_type = synapse_type;
	external_current = cell_pars["external_current"];
	reversal_potential_e = cell_pars["reversal_potential_e"];
	reversal_potential_i = cell_pars["reversal_potential_i"];
	tau_i = cell_pars["tau_i"];
	synapse_delay = cell_pars["synapse_delay"];
	current_noise = cell_pars["current_noise"];
	time_step = cell_pars["time_step"];
	rate_depression = cell_pars["rate_depression"];
	tau_recovery = cell_pars["tau_recovery"];
	inhib_ramp_factor = cell_pars["inhib_ramp_factor"];

	if (cell_type=="WB"){
		init_wb(p);
		init_cressman(p); // only some of the parameters are useful
	}
	else if (cell_type=="KC"){
		init_kc(p);
	}
	else if (cell_type=="ABCsimple"){
		init_abcsimple(p);
	}
	else if (cell_type=="Iz"){
		init_Iz(p);
	}
#if ANDERSON 
	else if (cell_type=="Anderson_sans_Ca"){
		init_anderson_sans_Ca(p);
	}
	else if (cell_type=="Anderson"){
		init_anderson(p);
	}
#endif
	else if (cell_type=="Cressman"){
		init_cressman(p);
	}
	else if (cell_type=="Ahmed_modified"){
		init_ahmed_modified(p);
	}
	else if (cell_type=="Froh"){
		init_frohlich(p);
	}
	else{
		cout<<"cellParameters construction: no matching cell type"<<endl;
		exit (1);
	}

	if (synapse_type=="VmD" || synapse_type=="VmDdiscont" || synapse_type=="VmD2discont"){
		init_vmd(p);
	}
	if (synapse_type=="VmDdiscont" || synapse_type=="VmD2discont" || synapse_type=="Normaldiscont" || synapse_type=="Poissondiscont"){
		init_discont(p);
	}
	if (synapse_type=="Poissondiscont"){
		init_poisson(p);
	}
	return;

}

//Assignment operator overloading
cellParameters& cellParameters::operator=(const cellParameters& right_side){
	if (this==&right_side){
		return *this;
	}
	else{
		//cout<<"cellP init; cell type is"<<right_side.cell_type<<endl;
		Fmap p = right_side.cell_pars;
		string t = right_side.cell_type;
		string s = right_side.synapse_type;
		//cout<<"Before P class init"<<endl;
		this->class_init(p, t, s);
		//cout<<"cellpara init before return this"<<endl;
		return *this;
	}
}

void cellParameters::init_wb(Fmap& cell_paras){
	leak_conductance=cell_paras["leak_conductance"];
  leak_potential=cell_paras["leak_potential"];
  sodium_conductance=cell_paras["sodium_conductance"];
  sodium_potential=cell_paras["sodium_potential"];
  potassium_conductance=cell_paras["potassium_conductance"];
  potassium_potential=cell_paras["potassium_potential"];
  sodium_half_activation=cell_paras["sodium_half_activation"];
  sodium_tau=cell_paras["sodium_tau"];
	alpha=cell_paras["alpha"]; // synaptic opening rate
	return;
}

void cellParameters::init_kc(Fmap& p){
	//tau_n for KC model (default=0ms)
	init_wb(p);
  tau_n=p["tau_n"];
	return;
}

void cellParameters::init_abcsimple(Fmap& p){
  //ABC simple model
  Vsn=p["Vsn"];
  Vreset=p["Vreset"];
  Vmax=p["Vmax"];
  Rreset=p["Rreset"];
  a=p["a"];
  bsn=p["bsn"];
  c=p["c"];
  d=p["d"];
  q=p["q"];
	l=p["l"];
	return;
  }

void cellParameters::init_vmd(Fmap& p){
	vmd_noise=p["vmd_noise"]; // ge
	vmd_noise2=p["vmd_noise2"]; //gi

  tau_e=p["tau_e"];
	tau_e2=p["tau_e2"];

  ge_0=p["ge_0"]; //ge
	ge2_0=p["ge2_0"]; // gi

	return;
}

void cellParameters::init_discont(Fmap& p){
	synapse_i_max=p["synapse_i_max"];
	//Vopens=p["Vopens"];
	Vopens=0;
	return;
}

void cellParameters::init_Iz(Fmap& p){
	init_abcsimple(p);
	Vr=p["Vr"];
	d_step=0;
	return;
}
#if ANDERSON
void cellParameters::init_anderson_sans_Ca(Fmap& p){
	aA=p["aA"];
	vhalfA=p["vhalfA"];
	aX=p["aX"];
	vhalfX=p["vhalfX"];
	aB=p["aB"];
	vhalfB=p["vhalfB"];
	aW=p["aW"];
	vhalfW=p["vhalfW"];
	am=p["am"];
	vhalfm=p["vhalfm"];
	lambda=p["lambda"];
	tauX=p["tauX"];
	tauB=p["tauB"];
	VNa=p["VNa"];
	VK=p["VK"];
	VL=p["VL"];
	gNa=p["gNa"];
	gK=p["gK"];
	gL=p["gL"];
	gA=p["gA"];
	ei=p["ei"];
	nt=p["nt"];	
	return;
}

void cellParameters::init_anderson(Fmap& p){
		init_anderson_sans_Ca(p);
		VCa=p["VCa"];
		gKCa=p["gKCa"];
		zCa=p["zCa"];
		Kc=p["Kc"];
		Kd=p["Kd"];
		PCa=p["PCa"];
		C0=p["C0"];
		Rp=p["Rp"];
		Kp=p["Kp"];
		b=p["b"];
		f=p["f"];
		B0=p["B0"];
		D=p["D"];
		length0=17e-4;
		length1=17.5e-4;
		length2=18e-4;
		length3=18.5e-4;
		length4=19e-4;
		length5=19.5e-4;
		length6=20e-4;
		return;
}
#endif

void cellParameters::init_poisson(Fmap& p){
	poisson_rate=p["poisson_rate"];
	poisson_increment=p["poisson_increment"];
	tau_e=p["tau_e"];
	ge_0=0;
	return;
}

void cellParameters::init_ahmed_modified(Fmap& p){
a = p["a"];
b = p["b"];
Vr = p["Vr"];
k = p["k"];
d_step = p["d_step"];
k_after = p["k_after"];
Vchangek = p["Vchangek"];
Vtshift = p["Vtshift"];
Vreset = p["Vreset"];
Rreset = p["Rreset"];
Vmax = p["Vmax"];
return;
}

void cellParameters::init_frohlich(Fmap& p){
init_wb(p);
cressman_D = p["cressman_D"];
frohlich_k_forward = p["frohlich_k_forward"];
frohlich_k_back = p["frohlich_k_back"];
frohlich_b_max = p["frohlich_b_max"];
frohlich_buffer = p["frohlich_buffer"];
frohlich_max_pump_current = p["frohlich_max_pump_current"];
frohlich_Keq_pump = p["frohlich_Keq_pump"];
frohlich_Keq_glia = p["frohlich_Keq_glia"];
frohlich_glia_rise_factor = p["frohlich_glia_rise_factor"];
stim_start = p["stim_start"];
stim_end = p["stim_end"];
stim_strength = p["stim_strength"]; 
gCa = p["gCa"];
gAHP = p["gAHP"];
Cl_intracellular = p["Cl_intracellular"];
Cl_extracellular = p["Cl_extracellular"];
VCa = p["VCa"];
return;
}


void cellParameters::init_cressman(Fmap& p){
cressman_beta = p["cressman_beta"];
cressman_rho = p["cressman_rho"];
cressman_D = p["cressman_D"]; // diffusion constant
cressman_G_glia = p["cressman_G_glia"];
cressman_varepsilon = p["cressman_varepsilon"];
cressman_k_zero_inf = p["cressman_k_zero_inf"];
gL = p["gL"];
gNa = p["gNa"];
gK = p["gK"];
gCa = p["gCa"];
gAHP = p["gAHP"];
gKL = p["gKL"];
gNaL = p["gNaL"];
gClL = p["gClL"];
VCl = p["VCl"];
VCa = p["VCa"];
Cl_intracellular = p["Cl_intracellular"];
Cl_extracellular = p["Cl_extracellular"];
return;
}

double cellParameters::get_values(string name){
  return cell_pars[name];
}


Fmap cellParameters::get_hash(){
  return cell_pars;
}
