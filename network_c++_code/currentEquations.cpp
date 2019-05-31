#include "currentEquations.h"
//#include "cellVariables.h"
//#include "cellParameters.h"
#include <string>
//#include <unordered_map>
#include <math.h>
#include <iostream>
//using namespace std;

double KC_current_balance(cellVariables& cell_vars, cellParameters& cell_pars){

  double membrane_p=cell_vars.membrane_potential;
  double potassium_ch=cell_vars.potassium_channel;
  double leak_c=cell_pars.leak_conductance;
  double leak_p=cell_pars.leak_potential;
  double sodium_c=cell_pars.sodium_conductance;
  double sodium_p=cell_pars.sodium_potential;
  double potassium_c=cell_pars.potassium_conductance;
  double potassium_p=cell_pars.potassium_potential;
  double iext=cell_pars.external_current;
  double sodium_ch=KC_sodium(membrane_p);
  // double sodium_ch=1/(1+exp((-20-membrane_p)/15));// separate equation


  double leak_df=leak_c*(membrane_p-leak_p);
  double sodium_df=sodium_c*(membrane_p-sodium_p)*sodium_ch;
  double potassium_df=potassium_c*(membrane_p-potassium_p)*potassium_ch;
  return -(leak_df+sodium_df+potassium_df)+iext;
}

double KC_potassium(cellVariables& cell_vars, cellParameters& cell_pars){
  double membrane_p=cell_vars.membrane_potential;
  double potassium_ch=cell_vars.potassium_channel;
  double tau_n=cell_pars.tau_n; // initially set to 1
  return (1/(1+exp((-25-membrane_p)/5))-potassium_ch)/tau_n;
}

double KC_sodium(const double& membrane_p){
 double result=1/(1+exp((-20-membrane_p)/15));
 return result;
}


double wb_current_balance(cellVariables& cell_vars,cellParameters& cell_pars){
	cressman_K_intracellular(cell_vars, cell_pars);
  double membrane_p=cell_vars.membrane_potential;
  double potassium_ch=cell_vars.potassium_channel;
  double wb_sodium_ch=cell_vars.sodium_channel;
  double leak_c=cell_pars.leak_conductance;
  //double leak_p=cell_pars.leak_potential;
  double sodium_c=cell_pars.sodium_conductance;
	//cout<<"MP is "<<membrane_p<<" Potassium ch is "<<potassium_ch<<endl;
 // double sodium_p=cell_pars.sodium_potential;
  double potassium_c=cell_pars.potassium_conductance;
	double gAHP=cell_pars.gAHP;
	//cout<<"gAHP is "<<gAHP;
	double Ca_i=cell_vars.Ca_intracellular;
	//cout<<"Ca_i "<<Ca_i<<endl;
  //double potassium_p=cell_pars.potassium_potential;
	double potassium_p=26.64*log(cell_vars.K_extracellular/cell_vars.K_intracellular);
	//cout<<"Potassium p is "<<potassium_p<<endl;
	//double sodium_p=26.64*log(cell_vars.Na_extracellular/cell_vars.Na_intracellular);
	//double potassium_p=-90;
	double sodium_p=55.0;
  double iext=cell_pars.external_current;
  double wb_sodium_ch_inf=wb_minf(membrane_p);
	//double leak_p=potassium_p;
	double leak_p=26.64*log((cell_vars.K_extracellular+0.085*130.0+0.1*8.0)/(cell_vars.K_intracellular+0.085*17.0+0.1*130.0));
  double leak_df=leak_c*(membrane_p-leak_p);
  double sodium_df=sodium_c*(membrane_p-sodium_p)*pow(wb_sodium_ch_inf,3)*wb_sodium_ch;
  double potassium_df=potassium_c*(membrane_p-potassium_p)*pow(potassium_ch,4);
	double ahp_df=gAHP*(Ca_i/(1.0+Ca_i))*(membrane_p-potassium_p);
	//cell_vars.potassium_current=-potassium_df-ahp_df-leak_df; // for cressman
//Dec 18, 2015
	double leak_factor=0.06;
	//if (cell_vars.K_extracellular>3.5){
	//	leak_factor=0.15;
	//}
	double leak_due_to_K_df = leak_factor*leak_c*(membrane_p-potassium_p); // (estimate 65% of total leak due to K)
	//cell_vars.potassium_current=-potassium_df-ahp_df; // for cressman (not including leak)
cell_vars.potassium_current=-potassium_df-ahp_df-leak_due_to_K_df; // for cressman (including leak due to K)
// End
	//cout<<"Potassium current is "<<	cell_vars.potassium_current<<endl;
	//cout<<"IN WB iext is "<<iext<<endl;
	//cout<<"IN WB potassium p is "<<potassium_p<<endl;
  return -(leak_df+sodium_df+potassium_df+ahp_df)+iext;
  return 0;
}


double wb_potassium(cellVariables& cell_vars, cellParameters& cell_pars){
  double membrane_potential=cell_vars.membrane_potential;
  double potassium_channel=cell_vars.potassium_channel;
  double ninf = wb_ninf(membrane_potential);
  double tau_n = wb_taun(membrane_potential);
  return (5*(ninf-potassium_channel)/tau_n);
}

double wb_sodium(cellVariables& cell_vars, cellParameters& cell_pars){
  double x=cell_vars.membrane_potential;
  double z=cell_vars.sodium_channel;
  double alphah=0.07*exp(-(x+58)/20);
  double betah=1/(exp(-0.1*(x+28))+1);

  return (5*(alphah*(1-z)-betah*z));
}


// x is membrane potential
double wb_minf(const double& x){
  double wb_alpham=-0.1*(x+35)/(exp(-0.1*(x+35))-1);
  double wb_betam=4*exp(-(x+60)/18);
  return (wb_alpham/(wb_alpham+wb_betam));
}

double wb_alphan(const double& x){
  return (-0.01*(x+34)/(exp(-0.1*(x+34))-1));
}

double wb_betan(const double& x){
  return 0.125*exp(-(x+44)/80);
}

double wb_ninf(const double& x){
  return wb_alphan(x)/(wb_alphan(x)+wb_betan(x));
}

double wb_taun(const double& x){
  return 1/(wb_alphan(x)+wb_betan(x));
}

double abc_simple_rc(cellVariables& cell_vars, cellParameters& cell_pars){
  double decay_rate = cell_pars.d;
  double rc = cell_vars.potassium_channel;
  return -rc*decay_rate;
}

double abc_simple(cellVariables& cell_vars, cellParameters& cell_pars){
  double x = cell_vars.membrane_potential;
  double vsn = cell_pars.Vsn;
  double rc = cell_vars.potassium_channel;
  //double vreset = cell_pars["Vreset"];
  //double vmax = cell_pars["Vmax"];
  double a = cell_pars.a;
  double bsn = cell_pars.bsn;
  double c = cell_pars.c;
  double q = cell_pars.q;
  // Added February 27th, 2010
  double l = cell_pars.l;
  //
  double iext = cell_pars.external_current;
  double diff = x-vsn;
  double diffsq = diff*diff;
  double diffquar = diffsq*diffsq;
  // Changed February 27th, 2010
  //return (c*(iext-bsn)+a*diffsq+q*diffquar)-rc;
  return (c*(iext-bsn)+a*diffsq+q*diffquar+l*fabs(diff))-rc;
  //

}
double Iz_current_balance(cellVariables& cell_vars, cellParameters& cell_pars){
	double x = cell_vars.membrane_potential;
	double u = cell_vars.potassium_channel;	
	double k = cell_pars.q;
	double vr = cell_pars.Vr;
	double vt = cell_pars.Vsn;
	double iext = cell_pars.external_current;
	double term1 = x-vr;
	double term2 = x-vt;
	if (x>vt){k=0.02;}

	return k*term1*term2-u+iext;
}

double Iz_u(cellVariables& cell_vars, cellParameters& cell_pars){
	double a = cell_pars.a;
	double b = cell_pars.d;
	double vr = cell_pars.Vr;
	double u = cell_vars.potassium_channel;
	double x = cell_vars.membrane_potential;
	double term1 = x-vr;

	return a*(b*term1-u);
}
//--------------------------------------------Anderson stuff--------------------------------------------------------------------//

// List of all Ca diffusion functions
//
#if ANDERSON
void anderson_get_Ca_variables(cellVariables& cell_vars, cellParameters& cell_pars, int i, double& Ca_U, double& Ca_C, double& length){
	switch (i){
		case 0:
			Ca_C=cell_vars.Ca_C0;
			Ca_U=cell_vars.Ca_U0;
			length=cell_pars.length0;
			break;
		case 1:
			Ca_C=cell_vars.Ca_C1;
			Ca_U=cell_vars.Ca_U1;
			length=cell_pars.length1;
			break;
		case 2:
			Ca_C=cell_vars.Ca_C2;
			Ca_U=cell_vars.Ca_U2;
			length=cell_pars.length2;
			break;
		case 3:
			Ca_C=cell_vars.Ca_C3;
			Ca_U=cell_vars.Ca_U3;
			length=cell_pars.length3;
			break;
		case 4:
			Ca_C=cell_vars.Ca_C4;
			Ca_U=cell_vars.Ca_U4;
			length=cell_pars.length4;
			break;
		case 5:
			Ca_C=cell_vars.Ca_C5;
			Ca_U=cell_vars.Ca_U5;
			length=cell_pars.length5;
			break;
		case 6:
			Ca_C=cell_vars.Ca_C6;
			Ca_U=cell_vars.Ca_U6;
			length=cell_pars.length6;
			break;
		default:
			cout<<"In anderson_return_Ca_variables:  Invalid Ca parameters"<<endl;
			exit (1);
	}
	return;
}

double anderson_Ca_U_total_derivative(cellVariables& cell_vars, cellParameters& cell_pars, int i){
	double b = cell_pars.b;
	double B0 = cell_pars.B0;
	double f = cell_pars.f;
	double Kp = cell_pars.Kp;
	double Ca_U, Ca_C, length;
	anderson_get_Ca_variables(cell_vars, cell_pars, i, Ca_U, Ca_C, length);
	return (-b*Ca_U + f*Ca_C*(1.0-Ca_U));
}

double anderson_Ca_U0_total_derivative(cellVariables& cell_vars, cellParameters& cell_pars){
	return anderson_Ca_U_total_derivative(cell_vars, cell_pars, 0);
}
double anderson_Ca_U1_total_derivative(cellVariables& cell_vars, cellParameters& cell_pars){
	return anderson_Ca_U_total_derivative(cell_vars, cell_pars, 1);
}
double anderson_Ca_U2_total_derivative(cellVariables& cell_vars, cellParameters& cell_pars){
	return anderson_Ca_U_total_derivative(cell_vars, cell_pars, 2);
}
double anderson_Ca_U3_total_derivative(cellVariables& cell_vars, cellParameters& cell_pars){
	return anderson_Ca_U_total_derivative(cell_vars, cell_pars, 3);
}
double anderson_Ca_U4_total_derivative(cellVariables& cell_vars, cellParameters& cell_pars){
	return anderson_Ca_U_total_derivative(cell_vars, cell_pars, 4);
}
double anderson_Ca_U5_total_derivative(cellVariables& cell_vars, cellParameters& cell_pars){
	return anderson_Ca_U_total_derivative(cell_vars, cell_pars, 5);
}
double anderson_Ca_U6_total_derivative(cellVariables& cell_vars, cellParameters& cell_pars){
	return anderson_Ca_U_total_derivative(cell_vars, cell_pars, 6);
}


double anderson_Ca_C_total_derivative(cellVariables& cell_vars, cellParameters& cell_pars, int i){
	double Ca_U_here, Ca_C_here, length_here;
	double b = cell_pars.b;
	double f = cell_pars.f;
	double B0 = cell_pars.B0;
	anderson_get_Ca_variables(cell_vars, cell_pars, i, Ca_U_here, Ca_C_here, length_here);
	double J_diff_here = anderson_Ca_C_partial_derivative(cell_vars, cell_pars, i);

	double term2 = b*B0*Ca_U_here-f*B0*Ca_C_here*(1-Ca_U_here);
	return J_diff_here + term2;
}
double anderson_Ca_C6_total_derivative(cellVariables& cell_vars, cellParameters& cell_pars){
	double Ca_U_here = cell_vars.Ca_U6;
	double Ca_C_here = cell_vars.Ca_C6;
	double ICa = cell_vars.ICa;
	double b = cell_pars.b;
	double f = cell_pars.f;
	double Rp = cell_pars.Rp;
	double Kp = cell_pars.Kp;
	double B0 = cell_pars.B0;

	//double term1 = -ICa*0.02;
	double term1 = -ICa*0.1063;
	double term2 = b*B0*Ca_U_here-f*B0*Ca_C_here*(1-Ca_U_here);
	double term3 = -Rp*(Ca_C_here/(Ca_C_here+Kp));
	double J_diff_here = anderson_Ca_C6_partial_derivative(cell_vars, cell_pars);
	return (term1 + term2 + term3 + J_diff_here);
	
}

double anderson_Ca_C5_total_derivative(cellVariables& cell_vars, cellParameters& cell_pars){
	return anderson_Ca_C_total_derivative(cell_vars, cell_pars, 5);
}
double anderson_Ca_C4_total_derivative(cellVariables& cell_vars, cellParameters& cell_pars){
	return anderson_Ca_C_total_derivative(cell_vars, cell_pars, 4);
}
double anderson_Ca_C3_total_derivative(cellVariables& cell_vars, cellParameters& cell_pars){
	return anderson_Ca_C_total_derivative(cell_vars, cell_pars, 3);
}
double anderson_Ca_C2_total_derivative(cellVariables& cell_vars, cellParameters& cell_pars){
	return anderson_Ca_C_total_derivative(cell_vars, cell_pars, 2);
}
double anderson_Ca_C1_total_derivative(cellVariables& cell_vars, cellParameters& cell_pars){
	return anderson_Ca_C_total_derivative(cell_vars, cell_pars, 1);
}
double anderson_Ca_C0_total_derivative(cellVariables& cell_vars, cellParameters& cell_pars){
	double Ca_U0 = cell_vars.Ca_U0;
	double Ca_C0 = cell_vars.Ca_C0;
	double b = cell_pars.b;
	double f = cell_pars.f;
	double B0 = cell_pars.B0;
	double J_diff0 = anderson_Ca_C_partial_derivative(cell_vars, cell_pars);
		
	double term2 = b*B0*Ca_U0-f*B0*Ca_C0*(1.0-Ca_U0);
	return J_diff0 + term2;
}


double anderson_Ca_C_partial_derivative(cellVariables& cell_vars, cellParameters& cell_pars){
	double diff_const = cell_pars.D;
	double Ca_C0 = cell_vars.Ca_C0;
	double Ca_C1 = cell_vars.Ca_C1;
	double length0 = cell_pars.length0;
	double length_diff0 = length0;
	double length_diffsq0 = length_diff0*length_diff0;
	double rCa_C0 = length0*Ca_C0;
	double term1 = 3.0*length0*(Ca_C1-Ca_C0)/length_diffsq0;
	return diff_const*term1/length0;
}

double anderson_Ca_C_partial_derivative(cellVariables& cell_vars, cellParameters& cell_pars, int i){
	double diff_const = cell_pars.D;
	double Ca_U_prev, Ca_C_prev, length_prev;
	double Ca_U_here, Ca_C_here, length_here;
	double Ca_U_next, Ca_C_next, length_next;
	anderson_get_Ca_variables(cell_vars, cell_pars, i-1, Ca_U_prev, Ca_C_prev, length_prev);
	anderson_get_Ca_variables(cell_vars, cell_pars, i, Ca_U_here, Ca_C_here, length_here);
	anderson_get_Ca_variables(cell_vars, cell_pars, i+1, Ca_U_next, Ca_C_next, length_next);

 	double length_diff = length_here - length_prev;
	double length_diffsq = length_diff*length_diff;
	double rCa_C_here = Ca_C_here*length_here;
	double rCa_C_prev = Ca_C_prev*length_here;
	double rCa_C_next = Ca_C_next*length_here;
	double term1 = (rCa_C_prev + rCa_C_next - 2*rCa_C_here)/length_diffsq;

	return diff_const*term1/length_here;
}

double anderson_Ca_C6_partial_derivative(cellVariables& cell_vars, cellParameters& cell_pars){
	double diff_const = cell_pars.D;
	double Ca_C_prev = cell_vars.Ca_C5;
	double Ca_C_here = cell_vars.Ca_C6;
	double length_prev = cell_pars.length5;
	double length_here = cell_pars.length6;
	double length_diff = length_here - length_prev;
	double length_diffsq = length_diff*length_diff;
	double rCa_C_here = Ca_C_here*length_prev;
	double rCa_C_prev = Ca_C_prev*length_prev;

	double term1 = -(rCa_C_here - rCa_C_prev)/length_diffsq;
	return diff_const*term1/length_here;
}


// List of all channel subsidary functions (Anderson)

double anderson_channel_infinity(cellVariables& cell_vars, cellParameters& cell_pars, string label){
	double V = cell_vars.membrane_potential;
	double a = anderson_get_a(cell_pars, label);
	double v_half = anderson_get_vhalf(cell_pars, label);
	//cout<<"In Anderson channel infinity"<<endl;
	//cout<<"Label is "<<label<<endl;
	//cout<<"V is "<<V<<endl;	
	//cout<<"a is "<<a<<endl;
	//cout<<"v half is "<<v_half<<endl;
	double term1 = 2.0*a*(V-v_half);
	double term2 = exp(-term1);
	double term3 = 1 + term2;
	//cout<<"result is "<<1/term3<<endl;
	return 1/term3;
}

double anderson_channel_tau(cellVariables& cell_vars, cellParameters& cell_pars, string label){
	if (label=="B"){return cell_pars.tauB;}
	else if (label=="X"){return cell_pars.tauX;}
		double V = cell_vars.membrane_potential;
		double lambda = cell_pars.lambda;
		double a = anderson_get_a(cell_pars, label);
		double vhalf = anderson_get_vhalf(cell_pars, label);
		double term1 = a*(V-vhalf);
		double term2 = 0.5/cosh(term1);
		return term2/lambda;
}

double anderson_channel_decay(cellVariables& cell_vars, cellParameters& cell_pars, string label){
	// To fix:
	double x = anderson_get_channel_decay_variable(cell_vars, label);
	//cout<<"label is "<<label<<"x is "<<x<<endl;
	double x_infinity = anderson_channel_infinity(cell_vars, cell_pars, label);
	double tau_infinity = anderson_channel_tau(cell_vars, cell_pars, label);
	return (x_infinity-x)/tau_infinity;
}

double anderson_get_channel_decay_variable(cellVariables& cell_vars, string label){
	if (label=="X"){
		return cell_vars.X;
	}
	else if (label=="B"){
		return cell_vars.B;
	}
	else if (label=="W"){
		return cell_vars.W;
	}
}

double anderson_get_a(cellParameters& cell_pars, string label){
	if (label=="W"){
		//cout<<"aW is "<<cell_pars["aW"]<<endl;
		return cell_pars.aW;
		}
	else if (label=="m"){
		//cout<<"am is "<<cell_pars["am"]<<endl;
		return cell_pars.am;
		}
	else if (label=="A"){
		//cout<<"aA is "<<cell_pars["aA"]<<endl;
		return cell_pars.aA;
		}
	else if (label=="B"){
		//cout<<"aB is "<<cell_pars["aB"]<<endl;
		return cell_pars.aB;
		}
	else if (label=="X"){
		//cout<<"aX is "<<cell_pars["aX"]<<endl;
		return cell_pars.aX;
	}
	return 0;
	}
	
double anderson_get_vhalf(cellParameters& cell_pars, string label){
	if (label=="W"){return cell_pars.vhalfW;}
	else if (label=="m"){return cell_pars.vhalfm;}
	else if (label=="A"){return cell_pars.vhalfA;}
	else if (label=="B"){return cell_pars.vhalfB;}
	else if (label=="X"){return cell_pars.vhalfX;}
	return 0;
}

// List of all Anderson decay functions
double anderson_w_decay(cellVariables& cell_vars, cellParameters& cell_pars){
	//cout<<"andersonn w decay"<<endl;
	return anderson_channel_decay(cell_vars, cell_pars, "W");
}

double anderson_b_decay(cellVariables& cell_vars, cellParameters& cell_pars){
	return anderson_channel_decay(cell_vars, cell_pars, "B");
}

double anderson_x_decay(cellVariables& cell_vars, cellParameters& cell_pars){
	return anderson_channel_decay(cell_vars, cell_pars, "X");
}


// List of all ionic currents (Anderson)
double anderson_KCa_current(cellVariables& cell_vars, cellParameters& cell_pars){
	double V = cell_vars.membrane_potential;
	double Ca_C6 = cell_vars.Ca_C6;
	double gKCa = cell_pars.gKCa;
	double VK = cell_pars.VK;
	double Kd = cell_pars.Kd;
	return gKCa*(Ca_C6/(Kd+Ca_C6))*(V-VK);

}

double anderson_K_current(cellVariables& cell_vars, cellParameters& cell_pars){
	double V = cell_vars.membrane_potential;
	double W = cell_vars.W;
	double gK = cell_pars.gK;
	double VK = cell_pars.VK;
	double W2 = W*W;
	double W4 = W2*W2;
	//cout<<"Getting into K current "<<endl;
	//cout<<"V is "<<V<<endl;
	//cout<<"W is "<<W<<endl;
	//cout<<"gK is "<<gK<<endl;
	//cout<<"W4 is "<<W4<<endl;
	double out = gK*W4*(V-VK);
	//cout<<"result is "<<out<<endl;
	return gK*W4*(V-VK);
}
double anderson_L_current(cellVariables& cell_vars, cellParameters& cell_pars){
	double V = cell_vars.membrane_potential;
	double gL = cell_pars.gL; //cout<<"gL is "<<gL<<endl;
	double VL = cell_pars.VL; //cout<<"VL is "<<VL<<endl;
	
	//cout<<"Getting into L current "<<endl;
	//cout<<"V is "<<V<<endl;
	//cout<<"gL is "<<gL<<endl;
	//cout<<"VL is "<<VL<<endl;
	double out = gL*(V-VL);

	//cout<<"Result is "<<out<<endl;
	return gL*(V-VL);

}

double anderson_A_current(cellVariables& cell_vars, cellParameters& cell_pars){
	double V = cell_vars.membrane_potential;
	double A_infinity = anderson_channel_infinity(cell_vars, cell_pars, "A");
	//cout<<"A infinity is "<<A_infinity<<endl;
	double B = cell_vars.B;
	//cout<<"B is "<<B<<endl;
	double VK = cell_pars.VK;
	double gA = cell_pars.gA;
	
	//cout<<"Getting into a current"<<endl;
	//cout<<"A_inf is "<<A_infinity<<endl;
	//cout<<"B is "<<B<<endl;
	//cout<<"VK is "<<VK<<endl;
	//cout<<"gA is "<<gA<<endl;

	double out =  gA*A_infinity*B*(V-VK);
	//cout<<"Result is "<<out<<endl;
	return gA*A_infinity*B*(V-VK);
}	

double anderson_sodium_current(cellVariables& cell_vars, cellParameters& cell_pars){
	double V = cell_vars.membrane_potential;
	double W = cell_vars.W;
	double gNa = cell_pars.gNa;
	double VNa = cell_pars.VNa;
	double m_infinity = anderson_channel_infinity(cell_vars, cell_pars, "m");
	double term1 = m_infinity*m_infinity*m_infinity;
	//cout<<"term 1 is "<<term1<<endl;
	//cout<<"W is "<<W<<endl;
	//cout<<"Getting into sodium current "<<endl;
	//cout<<"V is "<<V<<endl;
	//cout<<"W is "<<W<<endl;
	//cout<<"gNa is "<<gNa<<endl;
	//cout<<"VNa is "<<VNa<<endl;
	//cout<<"minf is "<<m_infinity<<endl;
	double out = gNa*term1*(1.0-W)*(V-VNa);
	//cout<<"Result is "<<out<<endl;
	//cout<<"jkajfksdlj "<<out<<endl;
	return out;
}

double anderson_calcium_current(cellVariables& cell_vars, cellParameters& cell_pars){
	double V = cell_vars.membrane_potential;
	//cout<<" mp IS "<<V<<endl;
	double X = cell_vars.X;
	double Ca_C6 = cell_vars.Ca_C6;
	double Kc = cell_pars.Kc;
	double zCa = cell_pars.zCa;
	double PCa = cell_pars.PCa;
	double C0 = cell_pars.C0;
	
	double term1 = Kc/(Kc + Ca_C6);
	double term2 = 15.05155*V;
	double term3 = exp(-0.078*V);
	
	double result = PCa*X*X*term1*term2*(Ca_C6-C0*term3)/(1-term3);
	cell_vars.ICa_next = result;
	return result;
}


double anderson_return_poisson_rate(cellVariables& cell_vars, cellParameters& cell_pars){
	return cell_pars.poisson_rate;
} 


// Anderson current balance
double anderson_current_balance(cellVariables& cell_vars, cellParameters& cell_pars){
	double calcium_df = anderson_calcium_current(cell_vars, cell_pars);
	//cout<<"ca is "<<calcium_df<<endl;
	double sodium_df = anderson_sodium_current(cell_vars, cell_pars);
	//cout<<"sodium is "<<sodium_df<<endl;
	double a_df = anderson_A_current(cell_vars, cell_pars);
	//cout<<"A is "<<a_df<<endl;
	double l_df = anderson_L_current(cell_vars, cell_pars);
	//cout<<"L is "<<l_df<<endl;
	double kca_df = anderson_KCa_current(cell_vars, cell_pars);
	//cout<<"Kca is "<<kca_df<<endl;
	double iext = cell_pars.external_current;
	//cout<<"Iext is "<<iext<<endl;
	double k_df = anderson_K_current(cell_vars, cell_pars);
	double out = -(calcium_df + sodium_df + a_df + l_df + kca_df + k_df) + iext;
	//cout<<"Out is "<<out<<endl;
	return out;

}

double anderson_sans_Ca_current_balance(cellVariables& cell_vars, cellParameters& cell_pars){
	double sodium_df = anderson_sodium_current(cell_vars, cell_pars);
	double a_df = anderson_A_current(cell_vars, cell_pars);
	double l_df = anderson_L_current(cell_vars, cell_pars);
	double k_df = anderson_K_current(cell_vars, cell_pars);
	//double kca_df = anderson_KCa_current(cell_vars, cell_pars);
	double iext = cell_pars.external_current;
	double total = -(sodium_df + a_df + l_df + k_df);
	//cout<<"Total current is "<<total<<endl;
	return -(sodium_df + a_df + l_df + k_df) + iext;

}
#endif
//-----Ahmed modified---------//
double ahmed_modified_current_balance(cellVariables& cell_vars, cellParameters& cell_pars){
	double iext = cell_pars.external_current;
	double k = cell_pars.k;
	double b = cell_pars.b;
	double k_after = cell_pars.k_after;
	double V_changek = cell_pars.Vchangek;
	double Vt_shift = cell_pars.Vtshift;
	double rc = cell_vars.potassium_channel;
	double V = cell_vars.membrane_potential;
	double Vr = cell_pars.Vr;
	double Vt = Vt_shift;
	//double Vt = Vr + Vt_shift - b/k;

	//cout<<"iext is "<<iext<<endl;
	//cout<<"k is "<<k<<endl;
	//cout<<"b is "<<b<<endl;
	//cout<<"k_after is "<<k_after<<endl;
	//cout<<"V_changek is "<<V_changek<<endl;
	//cout<<"Vt_shift is "<<Vt_shift<<endl;
	//cout<<"rc is "<<rc<<endl;
	//cout<<"V is "<<V<<endl;
	//cout<<"Vr is "<<Vr<<endl;
	//cout<<"Vt is "<<Vt<<endl;
	//cout<<endl;

	if (V>V_changek) k = k_after;
	//cell_pars.aX=k*(V-Vr)*(V-Vt)-rc+iext;
	return  k*(V-Vr)*(V-Vt)-rc+iext;
}

double ahmed_modified_u(cellVariables& cell_vars, cellParameters& cell_pars){
	double a = cell_pars.a;
	double b = cell_pars.b;
	double Vr = cell_pars.Vr;
	double rc = cell_vars.potassium_channel;
	double V = cell_vars.membrane_potential;
	return a*(b*(V-Vr)-rc);
}

//---------------Cressman----------------------//

double cressman_current_balance(cellVariables& cell_vars, cellParameters& cell_pars){
	double i_ext = cell_pars.external_current;
	cressman_ion_reversal_potential(cell_vars, cell_pars); // obtain the values of ionic reversal potentials
	cressman_K_intracellular(cell_vars, cell_pars);
	cressman_Na_extracellular(cell_vars, cell_pars);

	double sodium_current = cressman_sodium_current(cell_vars, cell_pars);
	double potassium_current = cressman_potassium_current(cell_vars, cell_pars);
	double ahp_current = cressman_ahp_current(cell_vars, cell_pars);
	cell_vars.potassium_current = potassium_current + ahp_current;

	double leak_current = cressman_leak_current(cell_vars, cell_pars);
	cell_vars.sodium_current = sodium_current;
	//cout<<"Na reversal is "<<cell_vars.VNa_var<<"\n";
	//cout<<"K reversal is "<<cell_vars.VK_var<<"\n";
	//cout<<"Ca reversal is "<<cell_pars.VCa<<"\n";
	//cout<<"Na current is "<<sodium_current<<"\n";
	//cout<<"K current is "<<potassium_current<<"\n";
	//cout<<"AHP current is "<<ahp_current<<"\n";
	//cout<<"Leak current is "<<leak_current<<"\n";

	return (sodium_current + potassium_current + ahp_current + leak_current) + i_ext;

}

void cressman_ion_reversal_potential(cellVariables& cell_vars, cellParameters& cell_pars){
	double Na_extracellular = cell_vars.Na_extracellular;
	double Na_intracellular = cell_vars.Na_intracellular;
	double K_extracellular = cell_vars.K_extracellular;
	double K_intracellular = cell_vars.K_intracellular;
	double Cl_extracellular = cell_pars.Cl_extracellular;
	double Cl_intracellular = cell_pars.Cl_intracellular;

	//cout<<"IN cressman ion reversal potential Na_e "<<Na_extracellular<<endl;
	//cout<<"IN cressman ion reversal potential Na_i "<<Na_intracellular<<endl;
	//cout<<"IN cressman ion reversal potential K_e "<<K_extracellular<<endl;
	//cout<<"IN cressman ion reversal potential K_i "<<K_intracellular<<endl;
	//cout<<"IN cressman ion reversal potential Cl_e "<<Cl_extracellular<<endl;
	//cout<<"IN cressman ion reversal potential Cl_i "<<Cl_intracellular<<endl;


	cell_vars.VNa_var = 26.64*log(Na_extracellular/Na_intracellular);
	cell_vars.VK_var = 26.64*log(K_extracellular/K_intracellular);
	//cell_vars.VCl_var = 26.64*log(Cl_extracellular/Cl_intracellular);
	cell_vars.VL_var = 26.64*log((K_extracellular + 0.065*Na_extracellular + 0.6*Cl_extracellular)/(K_intracellular + 0.065*Na_intracellular + 0.6*Cl_intracellular));

	return;
}

double cressman_sodium_current(cellVariables& cell_vars, cellParameters& cell_pars){
	double V = cell_vars.membrane_potential;
	double gNa = cell_pars.gNa;
	double VNa = cell_vars.VNa_var;
	double m_infinity = cressman_minf(V);
	double term1 = m_infinity*m_infinity*m_infinity;
	double sodium_channel = cell_vars.sodium_channel;

	//cout<<"IN Na current, m_infinity is "<<m_infinity<<endl;
	//cout<<"IN Na current, sodium_channel is "<<sodium_channel<<endl;
	//cout<<"IN Na current, gNa is"<<gNa<<endl;
	//cout<<"In Na current, term1 is"<<term1<<endl;
	//cout<<"In Na current, VNa is"<<VNa<<endl;
	//cout<<"IN Na current, V is "<<V<<endl;

	double out = -gNa*term1*sodium_channel*(V-VNa);
	return out;
}

double cressman_leak_current(cellVariables& cell_vars, cellParameters& cell_pars){
	double V = cell_vars.membrane_potential;
	double gNaL = cell_pars.gNaL;
	double gKL = cell_pars.gKL;
	double VNaL = cell_vars.VNa_var;
	double VKL = cell_vars.VK_var;
	double VClL = cell_pars.VCl;
	double gClL = cell_pars.gClL;
	//double VL = cell_vars.VL_var; // Notice VL is a variable depending on ion concentrations
	double leak_Na = -gNaL*(V - VNaL);
	double leak_K =  -gKL*(V - VKL);
	double leak_Cl = -gClL*(V - VClL);
	cell_vars.potassium_current+=leak_K;
	return (leak_Na + leak_K + leak_Cl);
}

double cressman_potassium_current(cellVariables& cell_vars, cellParameters& cell_pars){
	double V = cell_vars.membrane_potential;
	double VK = cell_vars.VK_var;
	double gK = cell_pars.gK;
	double potassium_channel = cell_vars.potassium_channel;
	double potassium_channel_4 = pow(potassium_channel, 4);
	//cout<<"In potassium current, V is "<<V<<endl;
	//cout<<"In potassium current, VK is "<<VK<<endl;
	//cout<<"In potassium current, gK is "<<gK<<endl;
	//cout<<"In potassium current, channel is "<<potassium_channel<<endl;
	return -gK*(V - VK)*potassium_channel_4;
}

double cressman_ahp_current(cellVariables& cell_vars, cellParameters& cell_pars){
	double VK = cell_vars.VK_var;
	double gAHP = cell_pars.gAHP;
	double Ca_i = cell_vars.Ca_intracellular;
	double V = cell_vars.membrane_potential;
	double ahp_current = -(V - VK)*(gAHP*Ca_i)/(1+Ca_i);
	//cout<<"gAHP is "<<gAHP<<endl;
	//cout<<"ahp current is "<<ahp_current<<endl;
	return ahp_current;
}

double cressman_sodium(cellVariables& cell_vars, cellParameters& cell_pars){
  double x=cell_vars.membrane_potential;
  double z=cell_vars.sodium_channel;
	//cout<<"In cressman sodium, x is "<<x<<endl;
	//cout<<"In cressman sodium, z is "<<z<<endl;
  double alphah=0.07*exp(-(x+44)/20);
  double betah=1/(exp(-0.1*(x+4))+1);
  double out=(3*(alphah*(1-z)-betah*z));
	//cout<<"In cressman sodium, out is "<<out<<endl;
	return out;
	// differential equation for h
}


// x is membrane potential
double cressman_minf(const double& x){
  double cressman_alpham=-0.1*(x+30)/(exp(-0.1*(x+30))-1);
  double cressman_betam=4*exp(-(x+55)/18);
  return (cressman_alpham/(cressman_alpham+cressman_betam));
}

// x is membrane potential
double cressman_potassium(cellVariables& cell_vars, cellParameters& cell_pars){
	double x=cell_vars.membrane_potential;
  double z=cell_vars.potassium_channel;

  double alphan=-0.01*(x+34)/(exp(-0.1*(x+34))-1);
  double betan=0.125*exp(-(x+44)/80);
  return  (3*(alphan*(1-z)-betan*z));
}


double cressman_Ca_total_derivative(cellVariables& cell_vars, cellParameters& cell_pars){//intracellular Ca2+
	double Ca_i = cell_vars.Ca_intracellular;
	double V = cell_vars.membrane_potential;
	double gCa = cell_pars.gCa;
	double VCa = cell_pars.VCa;
	//double tau_now_Ca = cell_vars.tau_now_Ca;
	//cout<<"VCA is "<<VCa<<endl;
	//cell_vars.tau_now_Ca+=change_Ca_removal_time_constant(cell_vars, cell_pars)*cell_pars.time_step;
	return -0.002*gCa*(V-VCa)/(1+exp(-(V+25)/2.5))-Ca_i/cell_vars.tau_now_Ca;
	}

double change_Ca_removal_time_constant(cellVariables& cell_vars, cellParameters& cell_pars){
	double tau_init_Ca = 80.0;
	double tau_final_Ca = 2000.0;
	double K_extracellular = cell_vars.K_extracellular;
	double tau_now_Ca = cell_vars.tau_now_Ca;
	if (cell_vars.current_time>12200){
		if (K_extracellular>3.5){
			return -(tau_now_Ca - tau_final_Ca)/50000.00;
		}
		else{
			return -(tau_now_Ca - tau_init_Ca)/20000.00;
		}
	}
	else{
		return 0;
	}

}



double cressman_K_total_derivative(cellVariables& cell_vars, cellParameters& cell_pars){//extracellular K+
	double beta = cell_pars.cressman_beta;
	double K_extracellular = cell_vars.K_extracellular;
	double K_current = cell_vars.potassium_current;
	double pump_current = cressman_pump_current(cell_vars, cell_pars);
	double glia_current = cressman_glia_current(cell_vars, cell_pars);
	double d_current = cressman_d_current(cell_vars, cell_pars);
	double diffusion_current = cell_pars.cressman_D*cell_vars.K_extracellular_diffusion_term;
	//cout<<"IN K total derivative, K_extracellular_diffusion_term is "<<diffusion_current<<endl;
	//cout<<"IN K total derivative, glia current is "<<glia_current<<endl;
	//cout<<"In K total derivative, pump_current is "<<pump_current<<endl;
	//cout<<"In K total derivative, beta is "<<beta<<endl;
	//cout<<"In K total derivative, K current is "<<K_current<<endl;

	
	return 0.001*(-0.33*K_current - 2*beta*pump_current - glia_current - d_current + diffusion_current);//0.001 factor to convert to (ms)^{-1}

}

double cressman_Na_total_derivative(cellVariables& cell_vars, cellParameters& cell_pars){//intracellular Na+
	double Na_current = cell_vars.sodium_current;
	double beta = cell_pars.cressman_beta;
	double pump_current = cressman_pump_current(cell_vars, cell_pars);

	//cout<<"In Na_total_derivative, Na current is "<<Na_current<<endl;
	//cout<<"In Na_total_derivative, beta is "<<beta<<endl;
	//cout<<"In Na_total_derivative, pump_current is "<<pump_current<<endl;
	return 0.001*(0.33*Na_current/beta - 3*pump_current);//0.001 factor to convert to (ms)^{-1}

}

double cressman_pump_current(cellVariables& cell_vars, cellParameters& cell_pars){
	double Na_intracellular = cell_vars.Na_intracellular;
	double rho = cell_pars.cressman_rho;
	double K_extracellular = cell_vars.K_extracellular;
	//cout<<"In Ipump, Na_i is "<<Na_intracellular<<endl;
	//cout<<"In Ipump, rho is "<<rho<<endl;
	//cout<<"In Ipump, K_e is "<<K_extracellular<<endl;

	return (rho/(1.0+exp((25-Na_intracellular)/3.0)))*(1.0/(1+exp(5.5-K_extracellular))); 
}

double cressman_glia_current(cellVariables& cell_vars, cellParameters& cell_pars){
	double G_glia = cell_pars.cressman_G_glia;
	double K_extracellular = cell_vars.K_extracellular;
	//cout<<"In Iglia, G_glia is "<<G_glia<<endl;
	//cout<<"In Iglia, K_e is "<<K_extracellular<<endl;
	return G_glia/(1+exp((18-K_extracellular)/2.5));
}

double cressman_d_current(cellVariables& cell_vars, cellParameters& cell_pars){
	double varepsilon = cell_pars.cressman_varepsilon;
	double K_extracellular = cell_vars.K_extracellular;
	double k_zero_inf = cell_pars.cressman_k_zero_inf;
	//cout<<"In d_current, varepsilon is "<<varepsilon<<endl;
	//cout<<"In d_current, K_e is "<<K_extracellular<<endl;
	//cout<<"In d_current, k_zero_inf is "<<k_zero_inf<<endl;
	return varepsilon*(K_extracellular-k_zero_inf);
}

void cressman_K_intracellular(cellVariables& cell_vars, cellParameters& cell_pars){
	double Na_intracellular = cell_vars.Na_intracellular;
	//cout<<"In K_intracellular, Na_i is "<<Na_intracellular<<endl;
	//cell_vars.K_intracellular = (140 + (18 - Na_intracellular));
	//cell_vars.K_intracellular = 130;
cell_vars.K_intracellular = 133.0;

	return;
}


void cressman_Na_extracellular(cellVariables& cell_vars, cellParameters& cell_pars){
	double Na_intracellular = cell_vars.Na_intracellular;
	double beta = cell_pars.cressman_beta;
	//cout<<"In Na_extracellular, Na_i is "<<Na_intracellular<<endl;
	//cout<<"In Na_extracellular, beta is "<<beta<<endl;
	cell_vars.Na_extracellular = (144 - beta*(Na_intracellular - 18.0));
	return;

}


//----Frohlich's glia and potassium pump----//

double frohlich_b_equation(cellVariables& cell_vars, cellParameters& cell_pars){

	double k_forward=cell_pars.frohlich_k_forward;
	double k_back=cell_pars.frohlich_k_back;
	double B_max=cell_pars.frohlich_b_max;
	double buffer=cell_vars.frohlich_buffer;
	double K_extracellular=cell_vars.K_extracellular;
	//cout<<"k_forward is "<<k_forward<<endl;
	//cout<<"k_back is "<<k_back<<endl;
	//cout<<"B max is "<<B_max<<endl;
	//cout<<"Buffer is "<<buffer<<endl;
	//cout<<"K extracellular is "<<K_extracellular<<endl;

return k_forward*(B_max-buffer)-frohlich_k_back(cell_vars, cell_pars)*K_extracellular*buffer;

}

double frohlich_k_back(cellVariables& cell_vars, cellParameters& cell_pars){

	double k_back = cell_pars.frohlich_k_back;
double K_oth = cell_pars.frohlich_Keq_glia;
double rise_factor = cell_pars.frohlich_glia_rise_factor;
double K_extracellular = cell_vars.K_extracellular;
//cout<<"rise factor is "<<rise_factor<<endl;
//cout<<"K_oth is "<<K_oth<<endl;

return k_back/(1+exp((K_extracellular-K_oth)/rise_factor));

}

double frohlich_g_equation(cellVariables& cell_vars, cellParameters& cell_pars){

double k_forward=cell_pars.frohlich_k_forward;
double B_max=cell_pars.frohlich_b_max;
double buffer=cell_vars.frohlich_buffer;
double K_extracellular=cell_vars.K_extracellular;
//cout<<"Frohlich g equation buffer is"<<buffer<<endl;
//cout<<"Frohlich g equation k forward is "<<k_forward<<endl;;
//cout<<"Frohlich g equation B max is "<<B_max<<endl;
//cout<<"Frohlich g equation K_extracellular is "<<K_extracellular<<endl;

return k_forward*(B_max-buffer)/1.1-frohlich_k_back(cell_vars,cell_pars)*K_extracellular*buffer;

}

double frohlich_K_pump(cellVariables& cell_vars, cellParameters& cell_pars){

	double I_max = cell_pars.frohlich_max_pump_current;
	double K_oeq = cell_pars.frohlich_Keq_pump;
	double K_extracellular=cell_vars.K_extracellular;
	//cout<<"Frolich K pump "<<K_extracellular<<endl;
	//cout<<"Frolich K pump "<<K_oeq<<endl;
	//cout<<"Frolich K pump "<<I_max<<endl;

	return I_max/pow((1+(K_oeq/K_extracellular)),2);
	//testing
	//return I_max;
}

double frohlich_K_total_derivative(cellVariables& cell_vars, cellParameters& cell_pars){//extracellular K+

	double K_extracellular = cell_vars.K_extracellular;
	double K_current = cell_vars.potassium_current;
	double pump_current = frohlich_K_pump(cell_vars, cell_pars);
	double glia_current = frohlich_g_equation(cell_vars, cell_pars);
	double diffusion_current = cell_pars.cressman_D*cell_vars.K_extracellular_diffusion_term;
	//cout<<"IN K total derivative, K_extracellular_diffusion_term is "<<diffusion_current<<endl;
	//cout<<"IN K total derivative, glia current is "<<glia_current<<endl;
	//cout<<"In K total derivative, pump_current is "<<pump_current<<endl;
	///cout<<"In K total derivative, beta is "<<beta<<endl;
	//cout<<"In K total derivative, K current is "<<K_current<<endl;

	//cout<<"IN frohlich K pump current is "<<pump_current<<endl;
	//cout<<"IN frohlich K glia current is "<<glia_current<<endl;
	//cout<<"IN frohlich K diffusion current is "<<diffusion_current<<endl;
	//cout<<"IN frohlich K current current is "<<K_current<<endl;
	cell_vars.glia_current=glia_current;
	cell_vars.pump_current=pump_current;

	
	double ans=-(50.0/96489.0)*(K_current+pump_current)+glia_current+0.001*diffusion_current;//unit mM/ms
	//double ans=-(20.0/96489.0)*(K_current+pump_current)+glia_current+0.001*diffusion_current;//unit mM/ms //testing

	//cout<<"IN K total derivative, ans is "<<ans<<endl;
	cell_vars.delta_K=ans;
	return ans;

	//return -(K_extracellular-4.0)/3.0; // pretend constant K

}

