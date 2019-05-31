#include "synapseEquations.h"
#include <string>
#include <unordered_map>
#include <math.h>
//#include <iostream.h>


double synapse_gating(cellVariables& cell_vars, cellParameters& cell_pars){
  double mp=cell_vars.membrane_potential;
  double gv=cell_vars.synapse_i;
  double tau_i=cell_pars.tau_i;
  double alpha=cell_pars.alpha; // initially set to 12
  //cout<<"In synapse gating; mp is "<<mp<<" gv is "<<gv<<" tau_i is "<<tau_i<<endl;
  return (alpha*(1.0-gv)/(1+exp(-mp/2))-gv/tau_i);
}

double synapse_gating_decay_part(cellVariables& cell_vars, cellParameters& cell_pars){
  //double mp=cell_vars["membrane_potential"];
  double gv=cell_vars.synapse_i;
  double tau_i=cell_pars.tau_i;
  //double alpha=cell_pars["alpha"]; // initially set to 12
  //cout<<"In synapse gating; mp is "<<mp<<" gv is "<<gv<<" tau_i is "<<tau_i<<endl;
  return -gv/tau_i;


}

double vmd_gating(cellVariables& cell_vars, cellParameters& cell_pars){
  double ge_0=cell_pars.ge_0;
  double tau_e=cell_pars.tau_e;
  double res=-(cell_vars.synapse_e-ge_0)/tau_e;
  return res;
}

double vmd_gating2(cellVariables& cell_vars, cellParameters& cell_pars){
  double ge_0=cell_pars.ge2_0;
  double tau_e=cell_pars.tau_e2;
  double res=-(cell_vars.synapse_e2-ge_0)/tau_e;
  return res;
}
