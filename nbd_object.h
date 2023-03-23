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
		double PE;
		double KE;
		double sugg_del_t;
		unsigned int t_block;
		vector<double> r;
		vector<double> v;
		vector<double> F_0;
		vector<double> F_1;
		vector<double> F_0_t0;
		vector<double> F_1_t0;


		nbd_object(long id_in,float m_in,vector<double> r_in, vector<double> v_in);
		void advance_timestep(double del_t);
		void primary_time_advance(double del_t,bool* start_flag, bool* in_current_block);
		void calculate_force(nbd_object *p2);
		void print_info();
	
	private:
		bool tblock_double_flip;
};

#endif

