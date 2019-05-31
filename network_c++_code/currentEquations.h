#include <string>
//#include <unordered_map>
#include "cellVariables.h"
#include "cellParameters.h"
using namespace std;

//typedef unordered_map<string,double> Fmap;
double KC_current_balance(cellVariables&, cellParameters&);
double KC_potassium(cellVariables&, cellParameters&);
double KC_sodium(const double&);

double wb_current_balance(cellVariables&, cellParameters&);
double wb_potassium(cellVariables&, cellParameters&);
double wb_sodium(cellVariables&, cellParameters&);
double wb_minf(const double& x);
double wb_alphan(const double& x);
double wb_betan(const double& x);
double wb_ninf(const double& x);
double wb_taun(const double& x);

double abc_simple(cellVariables&, cellParameters&);
double abc_simple_rc(cellVariables&, cellParameters&);
double Iz_current_balance(cellVariables&, cellParameters&);
double Iz_u(cellVariables&, cellParameters&);

#if ANDERSON
double anderson_Ca_U_total_derivative(cellVariables&, cellParameters&, int);
double anderson_Ca_U0_total_derivative(cellVariables&, cellParameters&);
double anderson_Ca_U1_total_derivative(cellVariables&, cellParameters&);
double anderson_Ca_U2_total_derivative(cellVariables&, cellParameters&);
double anderson_Ca_U3_total_derivative(cellVariables&, cellParameters&);
double anderson_Ca_U4_total_derivative(cellVariables&, cellParameters&);
double anderson_Ca_U5_total_derivative(cellVariables&, cellParameters&);
double anderson_Ca_U6_total_derivative(cellVariables&, cellParameters&);
double anderson_Ca_C_total_derivative(cellVariables&, cellParameters&, int i);
double anderson_Ca_C6_total_derivative(cellVariables&, cellParameters&);
double anderson_Ca_C5_total_derivative(cellVariables&, cellParameters&);
double anderson_Ca_C4_total_derivative(cellVariables&, cellParameters&);
double anderson_Ca_C3_total_derivative(cellVariables&, cellParameters&);
double anderson_Ca_C2_total_derivative(cellVariables&, cellParameters&);
double anderson_Ca_C1_total_derivative(cellVariables&, cellParameters&);
double anderson_Ca_C0_total_derivative(cellVariables&, cellParameters&);
double anderson_Ca_C_partial_derivative(cellVariables&, cellParameters&);
double anderson_Ca_C_partial_derivative(cellVariables&, cellParameters&, int i);
double anderson_Ca_C6_partial_derivative(cellVariables&, cellParameters&);
void anderson_get_Ca_variables(cellVariables& cell_vars, cellParameters& cell_pars, int i, double& Ca_U, double& Ca_C, double& length);

double anderson_get_channel_decay_variable(cellVariables& cell_vars, string label);
double anderson_channel_infinity(cellVariables&, cellParameters&, string label);
double anderson_channel_tau(cellVariables&, cellParameters&, string label);
double anderson_channel_decay(cellVariables&, cellParameters&, string label);
double anderson_get_a(cellParameters&, string label);
double anderson_get_vhalf(cellParameters&, string label);
double anderson_w_decay(cellVariables&, cellParameters&);
double anderson_b_decay(cellVariables&, cellParameters&);
double anderson_x_decay(cellVariables&, cellParameters&);
double anderson_KCa_current(cellVariables&, cellParameters&);
double anderson_L_current(cellVariables&, cellParameters&);
double anderson_A_current(cellVariables&, cellParameters&);
double anderson_sodium_current(cellVariables&, cellParameters&);
double anderson_calcium_current(cellVariables&, cellParameters&);
double anderson_return_poisson_rate(cellVariables&, cellParameters&);
double anderson_current_balance(cellVariables&, cellParameters&);
double anderson_sans_Ca_current_balance(cellVariables&, cellParameters&);
double anderson_K_current(cellVariables&, cellParameters&);
#endif


double ahmed_modified_current_balance(cellVariables&, cellParameters&);
double ahmed_modified_u(cellVariables&, cellParameters&);

// Cressman
double cressman_current_balance(cellVariables& cell_vars, cellParameters& cell_pars);
void cressman_ion_reversal_potential(cellVariables& cell_vars, cellParameters& cell_pars);
double cressman_sodium_current(cellVariables& cell_vars, cellParameters& cell_pars);
double cressman_leak_current(cellVariables& cell_vars, cellParameters& cell_pars);
double cressman_potassium_current(cellVariables& cell_vars, cellParameters& cell_pars);
double cressman_ahp_current(cellVariables& cell_vars, cellParameters& cell_pars);
double cressman_sodium(cellVariables& cell_vars, cellParameters& cell_pars);
double cressman_minf(const double& x);
double cressman_potassium(cellVariables& cell_vars, cellParameters& cell_pars);
double cressman_Ca_total_derivative(cellVariables& cell_vars, cellParameters& cell_pars);
double cressman_K_total_derivative(cellVariables& cell_vars, cellParameters& cell_pars);
double cressman_Na_total_derivative(cellVariables& cell_vars, cellParameters& cell_pars);
double cressman_pump_current(cellVariables& cell_vars, cellParameters& cell_pars);
double cressman_glia_current(cellVariables& cell_vars, cellParameters& cell_pars);
double cressman_d_current(cellVariables& cell_vars, cellParameters& cell_pars);
void cressman_K_intracellular(cellVariables& cell_vars, cellParameters& cell_pars);
void cressman_Na_extracellular(cellVariables& cell_vars, cellParameters& cell_pars);

double frohlich_K_total_derivative(cellVariables& cell_vars, cellParameters& cell_pars);
double frohlich_K_pump(cellVariables& cell_vars, cellParameters& cell_pars);
double frohlich_g_equation(cellVariables& cell_vars, cellParameters& cell_pars);
double frohlich_k_back(cellVariables& cell_vars, cellParameters& cell_pars);
double frohlich_b_equation(cellVariables& cell_vars, cellParameters& cell_pars);
double change_Ca_removal_time_constant(cellVariables& cell_vars, cellParameters& cell_pars);
