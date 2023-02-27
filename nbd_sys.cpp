#include "nbd_sys.h"
#include "nbd_object.h"
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>

using namespace std;

nbd_sys::nbd_sys(long num_objects_in,float mass_lower, float mass_upper,double max_size)
{
	uniform_real_distribution<double> mass_dist(mass_lower,mass_upper); 

	this->max_size=max_size;
	uniform_real_distribution<double> r_dist(-max_size,max_size); 
	uniform_real_distribution<double> v_dist(-1.5,1.5); 

	default_random_engine re;
	num_objects=num_objects_in;

	for(int i=0;i<num_objects_in;i++)
	{	
			nbd_object temp_star(i,mass_dist(re),{r_dist(re),r_dist(re),r_dist(re)},{v_dist(re),v_dist(re),v_dist(re)});
			stars.push_back(temp_star);	
	}

	this->time=0.0;


}

//TODO: Finish this to add input file for stars
nbd_sys::nbd_sys(FILE *infile)
{	
	if(infile==NULL)
		exit(EXIT_FAILURE);
	double mass;
	vector<double> r_in(3);
	vector<double> v_in(3);
	long id=0;
	double max_size=0.0;
	while (fscanf(infile, "%lf %lf %lf %lf %lf %lf %lf", &mass, &r_in[0],&r_in[1],&r_in[2],&v_in[0],&v_in[1],&v_in[2] )!=EOF) {
		nbd_object temp_star(id,mass,r_in,v_in);
		stars.push_back(temp_star);
		id++;
		if(max_size<norm(r_in))
		{
			max_size=norm(r_in);
		}
	}
	this->num_objects=id;
	this->time=0;
	this->max_size=max_size;
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
	for(int particle=0; particle<num_objects; particle++)
	{
		for(int particle_2=0; particle_2<num_objects; particle_2++)
		{
			if (particle_2!=particle) 
			{
				stars[particle].calculate_force(&stars[particle_2]);
			}

		}
	}

}


void nbd_sys::apply_force_updates(double del_t)
{
	double temp_max_r=0.0;
	for(int particle=0; particle<num_objects; particle++)
	{
		stars[particle].primary_time_advance(del_t);
		double temp_num=norm(stars[particle].r);
		if(temp_max_r<=temp_num)
		{temp_max_r=temp_num;}
	}
	this->max_size=temp_max_r;
}

void nbd_sys::apply_corrections(double del_t)
{
	double temp_max_r=0.0;
	for(int particle=0; particle<num_objects; particle++)
	{
		stars[particle].corrections(del_t);
		double temp_num=norm(stars[particle].r);
		if(temp_max_r<=temp_num)
		{temp_max_r=temp_num;}
	}
	this->max_size=temp_max_r;
	this->time+=del_t;

}


void nbd_sys::store_snapshot(FILE * outfile)
{
	for (auto it = this->stars.begin(); it != this->stars.end(); it++) 
		{
          nbd_object f=*it;
		  fprintf(outfile,"%8.3f,%7d,%7.8f,%7.12f,%7.12f,%8.10f\n",this->time, f.id, f.m, f.r[0], f.r[1],f.r[2]);
		}
}

void nbd_sys::external_potential()
{
}
