#ifndef NBD_SYS_H
#define NBD_SYS_H

#include <vector>
#include "lin_alg.h"
#include "nbd_object.h"
#include <stdio.h>

using namespace std;
/**
 * @brief N body system class that contains information about all objects and the system
 * 
 */
class nbd_sys
{
	public:
		long num_objects;
		vector<nbd_object> stars;
		double time;
		double max_size;

		nbd_sys(long num_objects_in,float mass_lower, float mass_upper,double max_size);
		nbd_sys(FILE *infile);
		void print_sys();
		void force_calculations();
		void apply_force_updates(double del_t);
		void apply_corrections(double del_t);
		void store_snapshot(FILE *outfile);
		void external_potential();

};

#endif

