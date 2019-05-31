#include "Diffusion.h"
#include "cortical_assignment.h"
#include "Diffusion_tools.h"
#include <stdlib.h>
// This is the class to find the *cell* neighbours in a 3D grid for diffusion equation calculation

double lapacian(Diffusion& diff_obj, double* variable_array, int cell_index){

	int cell_plus_x = diff_obj.cell_neighbour_plus_x_dir;
	int cell_minus_x = diff_obj.cell_neighbour_minus_x_dir;

	int cell_plus_y = diff_obj.cell_neighbour_plus_y_dir;
	int cell_minus_y = diff_obj.cell_neighbour_minus_y_dir;

	int cell_plus_z = diff_obj.cell_neighbour_plus_z_dir;
	int cell_minus_z = diff_obj.cell_neighbour_minus_z_dir;

	float dist_plus_x = diff_obj.dist_plus_x;
	float dist_minus_x = diff_obj.dist_minus_x;

	float dist_plus_y = diff_obj.dist_plus_y;
	float dist_minus_y = diff_obj.dist_minus_y;

	float dist_plus_z = diff_obj.dist_plus_z;
	float dist_minus_z = diff_obj.dist_minus_z;

	double variable_plus_x = assign_variable_value(variable_array, cell_plus_x, cell_index);
	double variable_minus_x = assign_variable_value(variable_array, cell_minus_x, cell_index);

	double variable_plus_y = assign_variable_value(variable_array, cell_plus_y, cell_index);
	double variable_minus_y = assign_variable_value(variable_array, cell_minus_y, cell_index);
	
	double variable_plus_z =  assign_variable_value(variable_array, cell_plus_z, cell_index);
	double variable_minus_z =  assign_variable_value(variable_array, cell_minus_z, cell_index);

	double variable_here = variable_array[cell_index];

	double second_derivative_x = calculate_second_derivative(dist_minus_x, dist_plus_x, variable_minus_x, variable_here, variable_plus_x);
	double second_derivative_y = calculate_second_derivative(dist_minus_y, dist_plus_y, variable_minus_y, variable_here, variable_plus_y);
	double second_derivative_z = calculate_second_derivative(dist_minus_z, dist_plus_z, variable_minus_z, variable_here, variable_plus_z);

	return (second_derivative_x + second_derivative_y + second_derivative_z);
}

double calculate_second_derivative(double x1, double x2, double v_prev, double v_central, double v_next){
	x1 = 1e-4*x1;
	x2 = 1e-4*x2; // 1um = 10^-4cm
	double first_derivative_1 = (v_central - v_prev)/x1;
	double first_derivative_2 = (v_next - v_central)/x2;
	return 2*(first_derivative_2 - first_derivative_1)/(x1+x2);	
}
	
double assign_variable_value(double* variable_array, int cell_index_next, int cell_index){
	if (cell_index_next < 0 || cell_index_next>=(XDIM*YDIM*CPC)){
		return variable_array[cell_index]; // No flux boundary condition (more realistic)
		//return 0; // absorbing boundary condition
	}
	else{
		return variable_array[cell_index_next];
	}
}





