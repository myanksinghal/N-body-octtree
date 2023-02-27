#include <iterator>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include "consts.h"
#include "nbd_object.h"
#include "lin_alg.h"
#include "nbd_sys.h"
#include "octtree.h"

using namespace std;


int main()
{
//printf("Test print G=%1.15f",G);
//vector<double> r_int={0.0,10.0,-5.6};
//vector<double> v_int={1.0,13,-1.6};

//nbd_object test(123,12.4,r_int,v_int);
/*test.print_info();
vector<double> f_in={1,0,0};
test.update_force(f_in);

printf("F is %8.3f,%8.4f,%8.4f\n",test.F[0],test.F[1],test.F[2]);
test.advance_timestep(10);
test.print_info();

printf("Start of Test sys with delT=1\nInitial:\n");
nbd_sys temp_sys(2,1.0,1.0);
temp_sys.print_sys();
for(int i=0;i<5;i++){
temp_sys.force_calculations();
printf("Test sys with delT=1\nStep:%d\n",i);
temp_sys.update(1.0);
temp_sys.print_sys();
}*/
//vector<double> r_int2={12.0435,11.0033,-11.6};
//vector<double> v_int2={1.0,13,-1.6};

//nbd_object test2(13,12.4,r_int2,v_int);

//nbd_sys temp_sys(100,1.0,10.0,1000.0);
FILE *infile;
infile=fopen("n6.pos", "r");
nbd_sys temp_sys(infile);
fclose(infile);
vector<double> cent={0.0,0.0,0.0};
FILE *outfile;

outfile=fopen("test_data_file.csv","w" );
fprintf(outfile,"time,id,mass,x,y,z\n");
double del_t=1;
for(int num_updates=0; num_updates<100;num_updates++)
{	
	temp_sys.store_snapshot(outfile);
	printf("Root size is %5.2f\n",temp_sys.max_size);
	OctTree *tree= new OctTree(cent,temp_sys.max_size);
	printf("Tree Creation success\n");	
	for(auto star_iterator=begin(temp_sys.stars); star_iterator!=end(temp_sys.stars); ++star_iterator)
	{
		tree->insert(&*star_iterator);

	}
	//tree->print_tree();
	for(auto it=temp_sys.stars.begin();it<temp_sys.stars.end();++it)
	{
		tree->traverse(&*it);
	}
	delete tree;
//	auto t=norm(temp_sys.stars[0].F);
//	printf("force of particle 0 is %8.8f", t);
//	temp_sys.stars[0].print_info();

	temp_sys.apply_force_updates(del_t);

	//temp_sys.stars[0].print_info();
	OctTree *tree_corrections= new OctTree(cent,temp_sys.max_size);
	printf("Tree Creation success\n");	
	for(auto star_iterator=begin(temp_sys.stars); star_iterator!=end(temp_sys.stars); ++star_iterator)
	{
		tree_corrections->insert(&*star_iterator);

	}
	//temp_sys.store_snapshot(outfile);
	//tree->print_tree();
	for(auto it=temp_sys.stars.begin();it<temp_sys.stars.end();++it)
	{
		tree_corrections->traverse(&*it);
	}
	delete tree_corrections;
	temp_sys.apply_corrections(del_t);





}
fclose(outfile);


return 0;
}
