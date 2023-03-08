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
			F_0={0.0,0.0,0.0};
			F_1={0.0,0.0,0.0};
			F_0_t0={0.0,0.0,0.0};
			F_1_t0={0.0,0.0,0.0};
			PE=0.0;
			KE=0.0;
			sugg_del_t=0.1;


}
void nbd_object::calculate_force(nbd_object *p2)
{
	vector<double> R_vec=elementwise_sum(r, scaling_vector(p2->r, -1));
	double R_mag=norm(R_vec);
	double softened_R_Mag=sqrt(pow(SOFTENING_PARAM,2)+pow(R_mag,2));
	vector<double> V_vec=elementwise_sum(v, scaling_vector(p2->v, -1));
	double V_mag=norm(V_vec);
	//printf("R_mag^3=%5.5f\n", pow(R_mag,3));

	//printf("R_vec is %3.3f,%3.3f,%3.3f\n",R_vec[0],R_vec[1],R_vec[2]);
	//printf("mass of second object is %2.2f",-p2->m);
	//printf("Scaling value is %5.5f",(-p2->m)/(pow(R_mag, 3)));
	
	vector<double> F_0=scaling_vector(R_vec,-p2->m/pow(softened_R_Mag, 3));
	double a=dot_product(R_vec, V_vec)/pow(softened_R_Mag,2);

	vector<double> F_1=scaling_vector(V_vec,-p2->m/pow(softened_R_Mag, 3));
	F_1=elementwise_sum(F_1, scaling_vector(F_0, -3*a));

	this->F_0=elementwise_sum(this->F_0, F_0);
	this->F_1=elementwise_sum(this->F_1, F_1);
	//
	//printf("Force in calc force %3.3f,%3.3f,%3.3f\n",F_0[0],F_0[1],F_0[2]);
	this->PE+=-G*this->m*p2->m/R_mag;

}


void nbd_object::primary_time_advance(double del_t,bool* start_flag)
{	
	vector<double> F_0=scaling_vector(this->F_0,G);
	vector<double> F_1=scaling_vector(this->F_1,G);
	//printf("Force in advance timestep is %3.10f,%3.10f,%3.10f\n",F[0],F[1],F[2]);
	vector<double> v_new=elementwise_sum(scaling_vector(F_0, del_t),scaling_vector(F_1,0.5*del_t*del_t));
	v_new=elementwise_sum(v_new, this->v);
	vector<double> r_new=scaling_vector(F_1,del_t*del_t*del_t/6);
	r_new=elementwise_sum(r_new,scaling_vector(F_0, 0.5*del_t*del_t));
	r_new=elementwise_sum(r_new, scaling_vector(this->v, del_t));
	r_new=elementwise_sum(r_new,this->r);

	//printf("new_r=%3.3f,%3.3f,%3.3f\n old_r=%3.3f,%3.3f,%3.3f\n",r_new[0],r_new[1],r_new[2],r[0],r[1],r[2]);
	this->r=r_new;
	this->v=v_new;
	if(!*start_flag)
	{
	vector<double> F_0_t0=scaling_vector(this->F_0_t0,G);
	vector<double> F_1_t0=scaling_vector(this->F_1_t0,G);

	vector<double> temp0=elementwise_sum(F_0_t0, scaling_vector(F_0, -1.0));
	temp0=scaling_vector(temp0, -6/(del_t*del_t));
	vector<double> temp1=elementwise_sum(F_1, scaling_vector(F_1_t0, 2));
	temp1=scaling_vector(temp1, -2/(del_t));

	vector<double> F_2_t0=elementwise_sum(temp0, temp1);

	temp0=elementwise_sum(F_0_t0, scaling_vector(F_0, -1.0));
	temp0=scaling_vector(temp0, 12/(del_t*del_t*del_t));
	temp1=elementwise_sum(F_1, F_1_t0);
	temp1=scaling_vector(temp1, 6/(del_t*del_t));

	vector<double> F_3_t0=elementwise_sum(temp0, temp1);


	//printf("Force in advance timestep is %3.10f,%3.10f,%3.10f\n",F[0],F[1],F[2]);
	vector<double> v_corr=elementwise_sum(scaling_vector(F_2_t0, del_t*del_t*del_t/6),scaling_vector(F_3_t0,del_t*del_t*del_t*del_t/24));
	this->v=elementwise_sum(v_corr, this->v);
	vector<double> r_corr=elementwise_sum(scaling_vector(F_2_t0,del_t*del_t*del_t*del_t/24),scaling_vector(F_3_t0,del_t*del_t*del_t*del_t*del_t/120));
	this->r=elementwise_sum(r_corr,this->r);
	//printf("new_r=%3.3f,%3.3f,%3.3f\n old_r=%3.3f,%3.3f,%3.3f\n",r_new[0],r_new[1],r_new[2],r[0],r[1],r[2]);
	//
	this->sugg_del_t=sqrt((neta*(norm(F_0)*norm(F_2_t0)+pow(norm(F_1),2)))/(norm(F_1)*norm(F_3_t0)+pow(norm(F_2_t0),2)));
	//printf("Suggested delta t is %4.8f\n",this->sugg_del_t);

	}
	this->F_0_t0=this->F_0;
	this->F_1_t0=this->F_1;
	this->F_0={0.0,0.0,0.0};
	this->F_1={0.0,0.0,0.0};
}

void nbd_object::corrections(double del_t)
{
	vector<double> F_0_t0=scaling_vector(this->F_0_t0,G);
	vector<double> F_1_t0=scaling_vector(this->F_1_t0,G);
	vector<double> F_0=scaling_vector(this->F_0, G);
	vector<double> F_1=scaling_vector(this->F_1, G);

	vector<double> temp0=elementwise_sum(F_0_t0, scaling_vector(F_0, -1.0));
	temp0=scaling_vector(temp0, -6/(del_t*del_t));
	vector<double> temp1=elementwise_sum(F_1, scaling_vector(F_1_t0, 2));
	temp1=scaling_vector(temp1, -2/(del_t));

	vector<double> F_2_t0=elementwise_sum(temp0, temp1);

	temp0=elementwise_sum(F_0_t0, scaling_vector(F_0, -1.0));
	temp0=scaling_vector(temp0, 12/(del_t*del_t*del_t));
	temp1=elementwise_sum(F_1, F_1_t0);
	temp1=scaling_vector(temp1, 6/(del_t*del_t));

	vector<double> F_3_t0=elementwise_sum(temp0, temp1);


	//printf("Force in advance timestep is %3.10f,%3.10f,%3.10f\n",F[0],F[1],F[2]);
	vector<double> v_corr=elementwise_sum(scaling_vector(F_2_t0, del_t*del_t*del_t/6),scaling_vector(F_3_t0,del_t*del_t*del_t*del_t/24));
	this->v=elementwise_sum(v_corr, this->v);
	vector<double> r_corr=elementwise_sum(scaling_vector(F_2_t0,del_t*del_t*del_t*del_t/24),scaling_vector(F_3_t0,del_t*del_t*del_t*del_t*del_t/120));
	this->r=elementwise_sum(r_corr,this->r);
	//printf("new_r=%3.3f,%3.3f,%3.3f\n old_r=%3.3f,%3.3f,%3.3f\n",r_new[0],r_new[1],r_new[2],r[0],r[1],r[2]);
	//
	this->sugg_del_t=sqrt((neta*(norm(F_0)*norm(F_2_t0)+pow(norm(F_1),2)))/(norm(F_1)*norm(F_3_t0)+pow(norm(F_2_t0),2)));
	this->F_0={0.0,0.0,0.0};
	this->F_1={0.0,0.0,0.0};
	this->F_0_t0={0.0,0.0,0.0};
	this->F_1_t0={0.0,0.0,0.0};


}

void nbd_object::print_info()
{
	printf("Particle ID %7d of mass %3.1f is at %3.3f,%3.3f,%3.3f with vel %3.3f,%3.3f,%3.3f\n",id,m,r[0],r[1],r[2],v[0],v[1],v[2]);
}
