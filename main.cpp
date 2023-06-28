#include <iterator>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include "consts.h"
#include "nbd_object.h"
#include "lin_alg.h"
#include "nbd_sys.h"
#include "octtree.h"
#include <omp.h>

bool increase_nstep = false;
bool decrease_nstep = false;
unsigned int n_time_blocks = starting_n_time_blocks;
double t1 = max_timestep / (pow(2, n_time_blocks));

using namespace std;

int main()
{

	//FILE *infile;
	//infile = fopen("n6.pos", "r");
	//nbd_sys temp_sys(infile);
	//fclose(infile);
	nbd_sys temp_sys(1000, 0.005, 10.0, 3.0,0.25);

	vector<double> cent = {0.0, 0.0, 0.0};
	FILE *outfile;

	outfile = fopen("test_data_file.csv", "w");
	fprintf(outfile, "time,id,mass,x,y,z,vx,vy,vz,KE,PE\n");
	double output_timestep=OUT_TIMESTEP;
	double integration_time = 100000;
	bool start_flag = true;
	double system_energy_start = 0.0;
	double previous_system_energy = 0.0;
	double current_system_energy = 0.0;
	double energy_error = 0.0;
	double del_t_min = t1;
	unsigned int current_block = 0;
	unsigned int energy_err_counter = 4;
	// for(int j=0;j<12;j++)
	while (integration_time > 0)
	{

		integration_time -= del_t_min;
		output_timestep -= del_t_min;
		OctTree *tree = new OctTree(cent, temp_sys.max_size);
		for (auto star_iterator = begin(temp_sys.stars); star_iterator != end(temp_sys.stars); ++star_iterator)
		{
			tree->insert(&*star_iterator);
		}

#pragma omp parallel for
		for (auto it = temp_sys.stars.begin(); it < temp_sys.stars.end(); ++it)
		{
			it->F_0_t0 = it->F_0;
			it->F_1_t0 = it->F_1;
			it->F_0 = {0.0, 0.0, 0.0};
			it->F_1 = {0.0, 0.0, 0.0};
			if(external_body)
				{temp_sys.external_potential(&*it);}
			tree->traverse(&*it);
		}
		delete tree;

		if (output_timestep < 0 && current_block == 0)
		{
			for (auto it = temp_sys.stars.begin(); it < temp_sys.stars.end(); ++it)
			{
				it->KE = (0.5) * (it->m) * pow(norm(it->v), 2);
			}
			temp_sys.store_snapshot(outfile);
			output_timestep = OUT_TIMESTEP;
		}

		if (true)
		{
			for (auto it = temp_sys.stars.begin(); it < temp_sys.stars.end(); ++it)
			{
				it->KE = (0.5) * (it->m) * pow(norm(it->v), 2);
			}
			temp_sys.system_energy();
			if (start_flag)
			{
				system_energy_start = temp_sys.total_KE + temp_sys.total_PE;
				previous_system_energy = system_energy_start;
			}
			current_system_energy = temp_sys.total_KE + temp_sys.total_PE;

			energy_error = (current_system_energy - system_energy_start);

			if (current_block == 0)
			{

				printf("system time is %3.10f\n", temp_sys.time);
				printf("Total Energy of system is %1.5f\n", temp_sys.total_KE + temp_sys.total_PE);
				printf("Energy Error of system is %1.10f\n", energy_error);
				printf("DE/E= %1.10f\n", energy_error / system_energy_start);
				printf("Q=%1.3f\n",2*temp_sys.total_KE/(-1*temp_sys.total_PE));
				if (abs(current_system_energy - previous_system_energy) >= DE_flag)
				{
					printf("energy error increase nstep\n");
					energy_err_counter++;
					if (energy_err_counter >= 5)
					{
						energy_err_counter = 0;
						increase_nstep = true;
					}
				}
			}
			previous_system_energy = temp_sys.total_KE + temp_sys.total_PE;
		}

		temp_sys.apply_force_updates(&start_flag, current_block);

		// printf("system time is %3.10f\n",temp_sys.time);
		// printf("Current TIME block is %d \n\n",current_block);
		current_block++;
		if (current_block % (n_time_blocks + 1) == 0)
		{
			current_block = 0;
			if (increase_nstep && n_time_blocks<max_t_blocks)
			{
				printf("Increase block \n");
				n_time_blocks++;
				t1 = max_timestep / (pow(2, n_time_blocks));
				printf("min del t is %1.7f\n", t1);
				increase_nstep = false;
				del_t_min = t1;
				for (auto it = temp_sys.stars.begin(); it < temp_sys.stars.end(); ++it)
				{
					it->t_block++;
				}
			}
			else if (decrease_nstep)
			{
				printf("Decreasing block \n");
				n_time_blocks--;
				t1 = max_timestep / (pow(2, n_time_blocks));
				printf("min del t is %1.7f\n", t1);
				decrease_nstep = false;
				del_t_min = t1;
			}
		}
	}
	fclose(outfile);

	return 0;
}
