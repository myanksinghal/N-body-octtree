#ifndef NBD_OBJECT_H
#define NBD_OBJECT_H

#include <vector>
#include "lin_alg.h"
#include <stdio.h>

using namespace std;
class nbd_object
{
	public:
		long id;
		float m;
		vector<double> r;
		vector<double> v;
		vector<double> F;

		nbd_object(long id_in,float m_in,vector<double> r_in, vector<double> v_in);
		void update_force(vector<double> const &F_new);
		void advance_timestep(double del_t);

		vector<double> calculate_force(nbd_object p2);
		void print_info();
};

#endif

