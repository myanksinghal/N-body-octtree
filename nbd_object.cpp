#include "nbd_object.h"
#include "lin_alg.h"
#include <math.h>
#include <vector>
#include "consts.h"
nbd_object::nbd_object(long id_in,float m_in,vector<double> r_in, vector<double> v_in)
{
			id=id_in;
			m=m_in;
			r=r_in;
			v=v_in;
			F={0.0,0.0,0.0};

}

void nbd_object::update_force(vector<double> const &F_new)
{
			F=elementwise_sum(F,F_new);
}

vector<double> nbd_object::calculate_force(nbd_object p2)
{
	vector<double> R_vec=elementwise_sum(r, scaling_vector(p2.r, -1));
	double R_mag=norm(R_vec);
	vector<double> V_vec=elementwise_sum(v, scaling_vector(p2.v, -1));
	double V_mag=norm(V_vec);

	vector<double> F_0=scaling_vector(R_vec,-p2.m/pow(R_mag, 3));
	//double a=dot_product(R_vec, V_vec)/pow(R_mag,2);

	//vector<double> F_1=scaling_vector(V_vec,-p2.m/pow(R_mag, 3));
	//F_1=elementwise_sum(F_1, scaling_vector(F_0, -3*a));
	//
	return F_0;

}


void nbd_object::advance_timestep(double del_t)
{	

	F=scaling_vector(F,m*G);
	vector<double> v_new=scaling_vector(F,del_t);
	v_new=elementwise_sum(v_new, v);
	vector<double> r_new=scaling_vector(F,0.5*del_t*del_t);
	r_new=elementwise_sum(r_new, scaling_vector(v, del_t));
	r_new=elementwise_sum(r_new,r);

	r=r_new;
	v=v_new;
	F={0.0,0.0,0.0};
}


void nbd_object::print_info()
{
	printf("Particle ID %7d of mass %3.1f is at %3.3f,%3.3f,%3.3f\n",id,m,r[0],r[1],r[2]);
}
