#include "nbd_sys.h"
#include "nbd_object.h"
#include <random>
nbd_sys::nbd_sys(long num_objects_in,float mass_lower, float mass_upper)
{
	uniform_real_distribution<double> mass_dist(mass_lower,mass_upper); 

	uniform_real_distribution<double> r_dist(-1000.0,1000.0); 
	uniform_real_distribution<double> v_dist(-10.0,10.0); 

	default_random_engine re;
	num_objects=num_objects_in;

	for(int i=0;i<num_objects_in;i++)
	{	
			nbd_object temp_star(i,mass_dist(re),{r_dist(re),r_dist(re),r_dist(re)},{v_dist(re),v_dist(re),v_dist(re)});
			stars.push_back(temp_star);	
	}


}

void nbd_sys::print_sys()
{
	for(int i=0;i<num_objects;i++)
	{
		stars[i].print_info();	
	}
}

void nbd_sys::force_calculations()
{
	vector<double> temp_force;
	for(int particle=0; particle<num_objects; particle++)
	{
		for(int particle_2=0; particle_2<num_objects; particle_2++)
		{
			if (particle_2!=particle) 
			{
				temp_force=stars[particle].calculate_force(stars[particle_2]);
				stars[particle].update_force(temp_force);
			}

		}
	}

}

void nbd_sys::update(double del_t)
{

	for(int particle=0; particle<num_objects; particle++)
	{
		stars[particle].advance_timestep(del_t);
	}

}


void nbd_sys::store_snapshot()
{



}
