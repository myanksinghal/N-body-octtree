#include "nbd_sys.h"
#include "nbd_object.h"
#include "consts.h"
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <omp.h>

using namespace std;

nbd_sys::nbd_sys(long num_objects_in,float mass_lower, float mass_upper,double max_size)
{
	uniform_real_distribution<double> mass_dist(mass_lower,mass_upper); 

	this->max_size=max_size;
	uniform_real_distribution<double> r_dist(-max_size,max_size); 
	uniform_real_distribution<double> v_dist(-0.0,0.0); 

	default_random_engine re;
	default_random_engine ve;
	num_objects=num_objects_in;

	for(int i=0;i<num_objects_in;i++)
	{	
			nbd_object temp_star(i,mass_dist(re),{r_dist(re),r_dist(re),r_dist(re)},{v_dist(ve),v_dist(ve),v_dist(ve)});
			stars.push_back(temp_star);	
	}
	scale_standard_units();
	this->time=0.0;
	this->total_KE=0.0;
	this->total_PE=0.0;
}

void nbd_sys::scale_standard_units()
{
	vector<double> COM_r={0.0,0.0,0.0};
	vector<double> COM_v={0.0,0.0,0.0}; 
	double total_mass=0.0;
	for (auto it = this->stars.begin(); it != this->stars.end(); it++) 
	{
		total_mass+=it->m;
		COM_r=elementwise_sum(COM_r,scaling_vector(it->r, it->m));
		COM_v=elementwise_sum(COM_v,scaling_vector(it->v, it->m));
	}
	COM_r=scaling_vector(COM_r,-1/total_mass);
	COM_v=scaling_vector(COM_v,-1/total_mass);

	for (auto it = this->stars.begin(); it != this->stars.end(); it++) 
	{
		it->r=elementwise_sum(it->r,COM_r);
		it->v=elementwise_sum(it->v,COM_v);
		it->m=it->m/total_mass;
	}

	



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
	//scale_standard_units();
	this->num_objects=id;
	this->time=0.0;
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


void nbd_sys::apply_force_updates(bool* start_flag,unsigned int current_block)
{
	double temp_max_r=0.0;
	int count=0;
	#pragma omp parallel for
	for(int particle=0; particle<num_objects; particle++)
	{
		if(current_block%stars[particle].t_block==0)
		{count++;
		stars[particle].primary_time_advance(t1*pow(2,stars[particle].t_block-1),start_flag);
		}
	}
	for(int particle=0; particle<num_objects; particle++)
	{
		//printf("Particle num %d in block %d\n", particle, stars[particle].t_block);	
		double temp_num=norm(stars[particle].r);
		if(temp_max_r<=temp_num)
		{temp_max_r=temp_num;}
	}	
	printf("Number of stars in this block: %d \n",count);
	this->max_size=temp_max_r+5.0;
	*start_flag=false;
	this->time+=t1;
}



void nbd_sys::store_snapshot(FILE * outfile)
{
	for (auto it = this->stars.begin(); it != this->stars.end(); it++) 
	{
		nbd_object f=*it;
		fprintf(outfile,"%8.8f,%7d,%7.8f,%7.12f,%7.12f,%7.12f,%7.12f,%7.12f,%7.12f,%8.10f,%8.10f\n",this->time, f.id, f.m, f.r[0], f.r[1],f.r[2],f.v[0], f.v[1],f.v[2],f.KE,f.PE);
	}
}

void nbd_sys::external_potential()
{
}


void nbd_sys::system_energy()
{
	this->total_KE=0.0;
	this->total_PE=0.0;
	for (auto it = this->stars.begin(); it != this->stars.end(); it++) 	
	{
		this->total_KE+=it->KE;
		it->KE=0.0;
		this->total_PE+=it->PE;
		it->PE=0.0;
	}
	this->total_PE=0.5*this->total_PE;
}