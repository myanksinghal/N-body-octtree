#include "nbd_sys.h"
#include "nbd_object.h"
#include "consts.h"
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <omp.h>

using namespace std;
extern bool decrease_nstep;

/**
 * @brief Construct a new nbd sys::nbd sys object if input file is not provided
 *
 * @param num_objects_in Number of objects in the simulation
 * @param mass_lower Lower limit of the mass distribution
 * @param mass_upper Upper limit of the mass distribution
 * @param max_size Maximum size of the simulation cube
 */
nbd_sys::nbd_sys(long num_objects_in, float mass_lower, float mass_upper, double max_size, double max_vel)
{
	


	this->max_size = max_size;
	normal_distribution<double> r_dist(0.0, max_size);
	normal_distribution<double> v_dist(0.0, max_vel);

	default_random_engine re;
	default_random_engine ve;
	num_objects = num_objects_in;

	double mass_array[num_objects_in];
	mass_dist_power_law(mass_lower,mass_upper,num_objects_in,mass_array);


	for (int i = 0; i < num_objects_in; i++)
	{
		nbd_object temp_star(i,mass_array[i], {r_dist(re), r_dist(re), r_dist(re)}, {v_dist(ve), v_dist(ve), v_dist(ve)});
		stars.push_back(temp_star);
	}
	this->time = 0.0;
	this->total_KE = 0.0;
	this->total_PE = 0.0;

	if(external_body)
		{
			this->external_star=new nbd_object(-1,ext_mass,ext_r,ext_v);
		};

	scale_standard_units();

}

void nbd_sys::mass_dist_power_law(float mass_lower, float mass_upper, long num_objects_in, double* mass_array)
{
	uniform_real_distribution<double> unif_dist(0, 1);
	default_random_engine uni;
	//double mass_array[num_objects_in];
	for (int i=0; i<num_objects_in; i++)
	{	//printf("%d\n",i);
		*(mass_array+i)= pow(unif_dist(uni)*(pow(mass_upper,alpha+1)-pow(mass_lower,alpha+1)) + pow(mass_lower,alpha+1) , 1/(alpha+1)) ;
	}

	//x = [(x1^(n+1) - x0^(n+1))*y + x0^(n+1)]^(1/(n+1))
	

}



/**
 * @brief Scale the system to standard N body units
 *
 */
void nbd_sys::scale_standard_units()
{
	vector<double> COM_r = {0.0, 0.0, 0.0};
	vector<double> COM_v = {0.0, 0.0, 0.0};
	double total_mass = 0.0;
	for (auto it = this->stars.begin(); it != this->stars.end(); it++)
	{
		total_mass += it->m;
		COM_r = elementwise_sum(COM_r, it->r, it->m);
		COM_v = elementwise_sum(COM_v, it->v, it->m);
	}
	COM_r = scaling_vector(COM_r, -1 / total_mass);
	COM_v = scaling_vector(COM_v, -1 / total_mass);

	for (auto it = this->stars.begin(); it != this->stars.end(); it++)
	{
		it->r = elementwise_sum(it->r, COM_r, 1);
		it->v = elementwise_sum(it->v, COM_v, 1);
		it->m = it->m / total_mass;
	}

	//Should not be done since this would result in the separation of the external potential from the center of the system 
	//this->external_star->r=elementwise_sum(this->external_star->r,COM_r,1);
	//this->external_star->v=elementwise_sum(this->external_star->v,COM_v,1);
}

/**
 * @brief Construct a new nbd sys::nbd sys object from input file
 *
 * @param infile Pointer to the input file
 */
nbd_sys::nbd_sys(FILE *infile)
{
	if (infile == NULL)
		exit(EXIT_FAILURE);
	double mass;
	vector<double> r_in(3);
	vector<double> v_in(3);
	long id = 0;
	double max_size = 0.0;
	while (fscanf(infile, "%lf %lf %lf %lf %lf %lf %lf", &mass, &r_in[0], &r_in[1], &r_in[2], &v_in[0], &v_in[1], &v_in[2]) != EOF)
	{
		nbd_object temp_star(id, mass, r_in, v_in);
		stars.push_back(temp_star);
		id++;
		if (max_size < norm(r_in))
		{
			max_size = norm(r_in);
		}
	}
	// scale_standard_units();
	this->num_objects = id;
	this->time = 0.0;
	this->max_size = max_size;

	if(external_body)
		{
			this->external_star=new nbd_object(-1,ext_mass,ext_r,ext_v);
		};
}

/**
 * @brief Print the information of objects inside the system
 *
 */
void nbd_sys::print_sys()
{
	for (int i = 0; i < num_objects; i++)
	{
		stars[i].print_info();
	}
}

/**
 * @brief Calculate the force between all the objects in the system with N^2 operations. Used only for testing and small systems.
 *
 */
void nbd_sys::force_calculations()
{
	for (int particle = 0; particle < num_objects; particle++)
	{
		for (int particle_2 = 0; particle_2 < num_objects; particle_2++)
		{
			if (particle_2 != particle)
			{
				stars[particle].calculate_force(&stars[particle_2]);
			}
		}
	}
}

/**
 * @brief Apply the force updates to all the objects in the system, based on the time block the flag for corrections is also passed.
 *
 * @param start_flag Flag to indicate if this is the first time step
 * @param current_block Current time block
 */
void nbd_sys::apply_force_updates(bool *start_flag, unsigned int current_block)
{
	double temp_max_r = 0.0;
	unsigned int count = 0;
	bool is_current_block = false;
#pragma omp parallel for
	for (int particle = 0; particle < num_objects; particle++)
	{
		is_current_block = current_block % stars[particle].t_block == 0;
		stars[particle].primary_time_advance(t1 * pow(2, stars[particle].t_block - 1), start_flag, &is_current_block);
		// printf("Particle num %d in block %d\n", particle, stars[particle].t_block);
	}
	for (int particle = 0; particle < num_objects; particle++)
	{
		if (current_block % stars[particle].t_block == 0)
			count++;
		double temp_num = norm(stars[particle].r);
		if (temp_max_r <= temp_num)
		{
			temp_max_r = temp_num;
		}
	}

	if (current_block == 1)
	{
		if (count == 0)
			decrease_nstep = true;
	}
	// printf("Corrector calc for %d particles\n",count);
	this->max_size = temp_max_r + 5.0;
	*start_flag = false;
	this->time += t1;
}

/**
 * @brief Store the snapshot of the system in a file
 *
 * @param outfile File pointer to the output file
 */
void nbd_sys::store_snapshot(FILE *outfile)
{
	for (auto it = this->stars.begin(); it != this->stars.end(); it++)
	{
		nbd_object f = *it;
		fprintf(outfile, "%8.8f,%7d,%7.8f,%7.12f,%7.12f,%7.12f,%7.12f,%7.12f,%7.12f,%8.10f,%8.10f\n", this->time, f.id, f.m, f.r[0], f.r[1], f.r[2], f.v[0], f.v[1], f.v[2], f.KE, f.PE);
	}
}

void nbd_sys::external_potential(nbd_object *in_star)
{
	in_star->calculate_force(this->external_star,true);
}

/**
 * @brief Calculate the total energy of the system using the stored values of KE and PE of each object
 *
 */
void nbd_sys::system_energy()
{
	this->total_KE = 0.0;
	this->total_PE = 0.0;
	for (auto it = this->stars.begin(); it != this->stars.end(); it++)
	{
		this->total_KE += it->KE;
		it->KE = 0.0;
		this->total_PE += it->PE;
		it->PE = 0.0;
	}
	this->total_PE = 0.5 * this->total_PE;
}