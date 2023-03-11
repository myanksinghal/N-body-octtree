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
double del_t=0.001;
double output_timestep=1;
double integration_time=10000;
bool start_flag=true;
double system_energy_start;
double energy_error=0.0;
while(integration_time>0)
{
	printf("Start value is %d\n",start_flag);	

	integration_time-=del_t;
	output_timestep-=del_t;
	printf("Root size is %5.2f\n",temp_sys.max_size);
	OctTree *tree= new OctTree(cent,temp_sys.max_size);
	for(auto star_iterator=begin(temp_sys.stars); star_iterator!=end(temp_sys.stars); ++star_iterator)
	{
		tree->insert(&*star_iterator);

	}
	printf("Tree Creation success\n");	
	
	#pragma omp parallel for
	for(auto it=temp_sys.stars.begin();it<temp_sys.stars.end();++it)
	{
		tree->traverse(&*it);
		it->KE=(0.5)*(it->m)*pow(norm(it->v),2);
	}
	delete tree;
	
	if(output_timestep<0)
	{
		temp_sys.store_snapshot(outfile);
		output_timestep=1;
	}

	temp_sys.system_energy();
	if(start_flag)
	{system_energy_start=temp_sys.total_KE+temp_sys.total_PE;}

	energy_error=(system_energy_start-temp_sys.total_KE-temp_sys.total_PE)/system_energy_start;

	

	printf("Total Energy of system is %1.5f\n",temp_sys.total_KE+temp_sys.total_PE);
	printf("Energy Error of system is %1.10f\n",energy_error);

	del_t=temp_sys.apply_force_updates(del_t,&start_flag);

	printf("Del t is %3.7f and system time is %3.4f\n \n",del_t,temp_sys.time);

}
fclose(outfile);


return 0;
}
