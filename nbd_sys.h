#ifndef NBD_SYS_H
#define NBD_SYS_H

#include <vector>
#include "lin_alg.h"
#include "nbd_object.h"
#include <stdio.h>

using namespace std;
class nbd_sys
{
	public:
		long num_objects;
		vector<nbd_object> stars;
		double time;
		double max_size;

		nbd_sys(long num_objects_in,float mass_lower, float mass_upper,double max_size);
		void print_sys();
		void force_calculations();
		void apply_force_updates(double del_t);
		void apply_corrections(double del_t);
		void store_snapshot(FILE *outfile);
};

#endif

