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
		
		nbd_sys(long num_objects_in,float mass_lower, float mass_upper);
		void print_sys();
		void force_calculations();
		void update(double del_t);
		void store_snapshot();

};

#endif

