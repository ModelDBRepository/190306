#include "LFP.h"
#include "cortical_assignment.h"
#include <stdlib.h>
#include "Diffusion.h"

// This is the class to find the *cell* neighbours in a 3D grid for diffusion equation calculation (assume each column has 16 cells 0-15)

Diffusion::Diffusion()
{
	//cell_neighbour_plus_x_dir = NULL;
	//cell_neighbour_minux_x_dir = NULL;
	//cell_neighbour_plus_y_dir = NULL;
	//cell_neighbour_minus_y_dir = NULL;
	//cell_neighbour_plus_z_dir = NULL;
	//cell_neighbour_minus_z_dir = NULL;
	return;
}

Diffusion::Diffusion(int cell_index)
{
	int x,y,cn;
	reverse_indexing_neurons(cell_index,x,y,cn);
	Diffusion::x = x;
	Diffusion::y = y;
	Diffusion::cn = cn;
	Diffusion::cell_index = cell_index;
	Diffusion::cpc = CPC;
	Diffusion::ecpc = ECPC;

	cell_neighbour_plus_x_dir = -1;
	cell_neighbour_minus_x_dir = -1;
	cell_neighbour_plus_y_dir = -1;
	cell_neighbour_minus_y_dir = -1;
	cell_neighbour_plus_z_dir = -1;
	cell_neighbour_minus_z_dir = -1;

	return;
}



void Diffusion::find_diffusion_neighbours()
{
		find_cell_plus_x();
		find_cell_minus_x();
		find_cell_plus_y();
		find_cell_minus_y();
		find_cell_plus_z();
		find_cell_minus_z();
		
	return;
}

void Diffusion::find_cell_plus_x(){
	int res;	
	double dist;
	int next_col_x = Diffusion::x + 1; // To add: exception at boundary
	int next_col_y = Diffusion::y;
	int cell_number = Diffusion::cell_index;

	switch(Diffusion::cn){
		case 0:
			res = cell_number+1;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 1:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 0);
			dist = WIDTH_BTW_COL;
			break;
		case 2:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 3);
			dist = WIDTH_BTW_COL;
			break;
		case 3:
			res = cell_number-1;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 4:
			res = cell_number+1;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 5:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 4);
			dist = WIDTH_BTW_COL;
			break;
		case 6:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 7);
			dist = WIDTH_BTW_COL;
			break;
		case 7:
			res = cell_number-1;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 8:
			res = cell_number+1;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 9:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 8);
			dist = WIDTH_BTW_COL;
			break;
		case 10:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 11);
			dist = WIDTH_BTW_COL;
			break;
		case 11:
			res = cell_number-1;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 12:
			res = cell_number+1;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 13:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 12);
			dist = WIDTH_BTW_COL;
			break;
		case 14:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 15);
			dist = WIDTH_BTW_COL;
			break;
		case 15:
			res = cell_number-1;
			dist = WIDTH_IN_COL_XY;
			break;
		default:
			res = -1;
			break;
	}
	cell_neighbour_plus_x_dir=res;
	dist_plus_x = dist;

	return;
}

void Diffusion::find_cell_minus_x(){
	int res;
	double dist;
	int next_col_x = Diffusion::x - 1; // To add: exception at boundary
	int next_col_y = Diffusion::y;
	int cell_number = Diffusion::cell_index;

	switch(cn){
		case 0:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 1);
			dist = WIDTH_BTW_COL;
			break;
		case 1:
			res = cell_number-1;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 2:
			res = cell_number+1;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 3:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 2);
			dist = WIDTH_BTW_COL;
			break;
		case 4:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 5);
			dist = WIDTH_BTW_COL;
			break;
		case 5:
			res = cell_number-1;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 6:
			res = cell_number+1;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 7:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 6);
			dist = WIDTH_BTW_COL;
			break;
		case 8:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 9);
			dist = WIDTH_BTW_COL;
			break;
		case 9:
			res = cell_number-1;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 10:
			res = cell_number+1;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 11:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 10);
			dist = WIDTH_BTW_COL;
			break;
		case 12:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 13);
			dist = WIDTH_BTW_COL;
			break;
		case 13:
			res = cell_number-1;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 14:
			res = cell_number+1;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 15:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 14);
			dist = WIDTH_BTW_COL;
			break;
		default:
			res = -1;
			break;
	}
	cell_neighbour_minus_x_dir = res;
	dist_minus_x = dist;
	return;
}

				
void Diffusion::find_cell_plus_y(){
	int res;
	double dist;
	int next_col_x = Diffusion::x; // To add: exception at boundary
	int next_col_y = Diffusion::y + 1;
	int cell_number = Diffusion::cell_index;


	switch(cn){
		case 0:
			res = cell_number+3;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 1:
			res = cell_number+1;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 2:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 1);
			dist = WIDTH_BTW_COL;
			break;
		case 3:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 0);
			dist = WIDTH_BTW_COL;
			break;
		case 4:
			res = cell_number+3;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 5:
			res = cell_number+1;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 6:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 5);
			dist = WIDTH_BTW_COL;
			break;
		case 7:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 4);
			dist = WIDTH_BTW_COL;
			break;
		case 8:
			res = cell_number+3;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 9:
			res = cell_number+1;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 10:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 9);
			dist = WIDTH_BTW_COL;
			break;
		case 11:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 8);
			dist = WIDTH_BTW_COL;
			break;
		case 12:
			res = cell_number+3;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 13:
			res = cell_number+1;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 14:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 13);
			dist = WIDTH_BTW_COL;
			break;
		case 15:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 12);
			dist = WIDTH_BTW_COL;
			break;
		default:
			res = -1;
			dist = WIDTH_BTW_COL;
			break;
	}
	cell_neighbour_plus_y_dir = res;
	dist_plus_y = dist;
	return;
}



void Diffusion::find_cell_minus_y(){
	int res;
	double dist;
	int next_col_x = Diffusion::x; // To add: exception at boundary
	int next_col_y = Diffusion::y - 1;
	int cell_number = Diffusion::cell_index;

	switch(cn){
		case 0:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 3);
			dist = WIDTH_BTW_COL;
			break;
		case 1:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 2);
			dist = WIDTH_BTW_COL;
			break;
		case 2:
			res = cell_number - 1;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 3:
			res = cell_number - 3;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 4:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 7);
			dist = WIDTH_BTW_COL;
			break;
		case 5:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 6);
			dist = WIDTH_BTW_COL;
			break;
		case 6:
			res = cell_number - 1;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 7:
			res = cell_number + 3;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 8:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 11);
			dist = WIDTH_BTW_COL;
			break;
		case 9:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 10);
			dist = WIDTH_BTW_COL;
			break;
		case 10:
			res = cell_number - 1;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 11:
			res = cell_number - 3;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 12:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 15);
			dist = WIDTH_BTW_COL;
			break;
		case 13:
			forward_indexing_neurons_diffuse(res, next_col_x, next_col_y, 14);
			dist = WIDTH_BTW_COL;
			break;
		case 14:
			res = cell_number-1;
			dist =  WIDTH_IN_COL_XY;
			break;
		case 15:
			res = cell_number+3;
			dist =  WIDTH_IN_COL_XY;
			break;
		default:
			res = -1;
			break;
	}
	cell_neighbour_minus_y_dir = res;
	dist_minus_y = dist;
	return;
}

void Diffusion::find_cell_plus_z(){

	int res;
	double dist;
	int next_col_x = Diffusion::x; // To add: exception at boundary
	int next_col_y = Diffusion::y;
	int cell_number = Diffusion::cell_index;

	switch(cn){
		case 0:
			res = -1;
			dist = WIDTH_IN_COL_Z; 
			break;
		case 1:
			res = -1;
			dist = WIDTH_IN_COL_Z; 
			break;
		case 2:
			res = -1;
			dist = WIDTH_IN_COL_Z; 
			break;
		case 3:
			res = -1;
			dist = WIDTH_IN_COL_Z; 
			break;
		default:
			res = cell_number-4;
			dist = WIDTH_IN_COL_Z; 
			break;
	}
	if (cn>11) dist = WIDTH_BTW_LAYER;
	cell_neighbour_plus_z_dir = res;
	dist_plus_z = dist;
	return;
}


void Diffusion::find_cell_minus_z(){
	int res;
	double dist;
	int next_col_x = Diffusion::x; // To add: exception at boundary
	int next_col_y = Diffusion::y;
	int cell_number = Diffusion::cell_index;


	switch(cn){
		case 12:
			res = -1;
			dist = WIDTH_IN_COL_Z; 
			break;
		case 13:
			res = -1;
			dist = WIDTH_IN_COL_Z; 
			break;
		case 14:
			res = -1;
			dist = WIDTH_IN_COL_Z; 
			break;
		case 15:
			res = -1;
			dist = WIDTH_IN_COL_Z; 
			break;
		default:
			res = cell_number+4;
			dist = WIDTH_IN_COL_Z; 
		break;
	}
	if (cn>=8 && cn<12) dist = WIDTH_BTW_LAYER;

	cell_neighbour_minus_z_dir = res;
	dist_minus_z = dist;
	return;
}


