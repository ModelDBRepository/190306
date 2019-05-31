#ifndef LFP_H
#define LFP_H

#include <stdlib.h>
#include <stdio.h>

using namespace std;
class LFP
{
	public:
		LFP(); // default constructor
		LFP(int x, int y, double radius);
		~LFP();
		LFP& operator=(const LFP& right_side); // assignment operator
		void count_and_fetch_neighbours();
		void pick_all_start_indicies();
	  void calculate_LFP(double* potentials, float& l23, float& l4, float& l5, float& l6, float&);

	private:
		int pick_start_index_from_one_xy(int x, int y);
		int count_neighbours();
		void fetch_neighbours();
	private:
		double radius;
		int x;
		int y;
		int n_neighbours;
		int* neighbour_x;
		int* neighbour_y;
		int* start_indicies;
		int neighbour_array_set;
		int start_indicies_array_set;
};
#endif /* LFP_H */

