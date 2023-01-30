#include <stdio.h>
#include <vector>
#include "consts.h"
#include "nbd_object.h"
#include "lin_alg.h"
#include "nbd_sys.h"
#include "octtree.h"
using namespace std;


int main()
{
//printf("Test print G=%1.15f",G);
vector<double> r_int={1.0435,11.0033,-5.6};
vector<double> v_int={1.0,13,-1.6};

nbd_object test(123,12.4,r_int,v_int);
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
vector<double> r_int2={12.0435,11.0033,-11.6};
vector<double> v_int2={1.0,13,-1.6};

nbd_object test2(13,12.4,r_int2,v_int);

vector<double> cent={0.0,0.0,0.0};

bounding_box RootBox(cent,100.0);
printf("Contains point%d\n",RootBox.containsPoint(test));
OctTree *tree= new OctTree(cent,100.0);
printf("Tree Creation success\n");
printf("BB box in tree has HW = %5.3f\n",tree->boundary.halfDimension);
if(tree->oct_1==nullptr)
{
printf("Contains no child node");
}
printf("Adding one point\n");
tree->insert(&test);
printf("adding finished\n");

//tree->subdivide();
if(tree->oct_1==nullptr)
{
printf("Contains no child node");
}
printf("No. of points is now %d\n",tree->node_stars.size());
printf("Adding second point");
tree->insert(&test2);
tree->print_tree();

/*
printf("R is %8.3f,%8.4f,%8.4f\n",r_int[0],r_int[1],r_int[2]);
printf("V is %8.3f,%8.4f,%8.4f\n",v_int[0],v_int[1],v_int[2]);
vector<double> test_vec=elementwise_sum(r_int, v_int);
printf("Lin elementwise sum is %8.3f,%8.4f,%8.4f\n",test_vec[0],test_vec[1],test_vec[2]);
test_vec=scaling_vector(r_int, 10.0);
printf("Lin elementwise sscale is %8.3f,%8.4f,%8.4f\n",test_vec[0],test_vec[1],test_vec[2]);
*/

return 0;
}
