#ifndef DIFFUSION_H
#define DIFFUSION_H

#include <stdlib.h>
#include <stdio.h>

using namespace std;
class Diffusion
{
	public:
		Diffusion(); // default constructor
		Diffusion(int);
	 void	find_diffusion_neighbours();
		
	private:
		void find_cell_plus_x();
		void find_cell_minus_x();
		void find_cell_plus_y();
		void find_cell_minus_y();
		void find_cell_plus_z();
		void find_cell_minus_z();

	public:
		int x,y,cn,cell_index,cpc,ecpc;
		
		int cell_neighbour_plus_x_dir;
		int cell_neighbour_minus_x_dir;
		int cell_neighbour_plus_y_dir;
		int cell_neighbour_minus_y_dir;
		int cell_neighbour_plus_z_dir;
		int cell_neighbour_minus_z_dir;
		double dist_plus_x, dist_minus_x, dist_plus_y, dist_minus_y, dist_plus_z, dist_minus_z;

};
#endif /* DIFFUSION_H */

