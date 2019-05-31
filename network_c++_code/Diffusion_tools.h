#include "Diffusion.h"
#include "cortical_assignment.h"

#ifndef DIFFUSION_TOOLS_H
#define DIFFUSION_TOOLS_H
using namespace std;
// This is the class to find the *cell* neighbours in a 3D grid for diffusion equation calculation

double lapacian(Diffusion& diff_obj, double* variable_array, int cell_index);
double calculate_second_derivative(double x1, double x2, double v_prev, double v_central, double v_next);
double assign_variable_value(double* variable_array, int cell_index_next, int cell_index);
#endif	




