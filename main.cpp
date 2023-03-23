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

using namespace std;


int main()
{

//nbd_sys temp_sys(2,1.0,1.0,10.0);
FILE *infile;
infile=fopen("n6.pos", "r");
nbd_sys temp_sys(infile);
fclose(infile);
vector<double> cent={0.0,0.0,0.0};
FILE *outfile;

outfile=fopen("test_data_file.csv","w" );
fprintf(outfile,"time,id,mass,x,y,z,vx,vy,vz,KE,PE\n");
double output_timestep=0.01;
double integration_time=10000;
bool start_flag=true;
double system_energy_start;
double energy_error=0.0;
double del_t_min=t1;
unsigned int current_block=0;
//for(int j=0;j<12;j++)
while(integration_time>0)
{


	//printf("Start value is %d\n",start_flag);	
	
	integration_time-=del_t_min;
	output_timestep-=del_t_min;
//	printf("Root size is %5.2f\n",temp_sys.max_size);
	OctTree *tree= new OctTree(cent,temp_sys.max_size);
	for(auto star_iterator=begin(temp_sys.stars); star_iterator!=end(temp_sys.stars); ++star_iterator)
	{
		tree->insert(&*star_iterator);

	}
//	printf("Tree Creation success\n");
		
	
	//unsigned int count=0;
	#pragma omp parallel for
	for(auto it=temp_sys.stars.begin();it<temp_sys.stars.end();++it)
	{
		//if(current_block%it->t_block==0)
		//{
	//	count++;
		it->F_0_t0=it->F_0;
		it->F_1_t0=it->F_1;
		it->F_0={0.0,0.0,0.0};
		it->F_1={0.0,0.0,0.0};
		tree->traverse(&*it);
		//}
	}
	//printf("Number of objects with force updates= %d\n", count);
	delete tree;
	
	if(output_timestep<0 && current_block==0)
	{
		for(auto it=temp_sys.stars.begin();it<temp_sys.stars.end();++it)
			{
					it->KE=(0.5)*(it->m)*pow(norm(it->v),2);
			}
		temp_sys.store_snapshot(outfile);
		output_timestep=0.01;
	}

	if(true)
	{
		for(auto it=temp_sys.stars.begin();it<temp_sys.stars.end();++it)
		{
				it->KE=(0.5)*(it->m)*pow(norm(it->v),2);
		}
		temp_sys.system_energy();
		if(start_flag)
		{system_energy_start=temp_sys.total_KE+temp_sys.total_PE;}

		energy_error=(system_energy_start-temp_sys.total_KE-temp_sys.total_PE);

		
		if(current_block==0)
		{
		printf("Total Energy of system is %1.5f\n",temp_sys.total_KE+temp_sys.total_PE);
		printf("Energy Error of system is %1.10f\n",energy_error);
		}
	}

	temp_sys.apply_force_updates(&start_flag,current_block);

	printf("system time is %3.10f\n",temp_sys.time);
	printf("Current TIME block is %d \n\n",current_block);
	current_block++;
	if (current_block%(n_time_blocks+1)==0)
	{
		current_block=0;
	}

}
fclose(outfile);


return 0;
}
