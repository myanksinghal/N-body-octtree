#include "octtree.h"
#include "lin_alg.h"
#include "nbd_object.h"
#include <cstdio>
#include <memory>
#include <vector>

using namespace std;


bounding_box::bounding_box(vector<double> center_in,float halfDimension_in)
{
center=center_in;
halfDimension=halfDimension_in;
}

bounding_box::bounding_box()
{
vector<double> tem={-1,-1,-1};
center=tem;
halfDimension=-1;
}

bool bounding_box::containsPoint(nbd_object P)
{
	P.print_info();
	printf("Cent %2.2f\n",center[0]);
	printf("Cent %2.2f\n",center[1]);
	printf("Cent %2.2f\n",center[2]);
	printf("r0 %2.2f\n",P.r[0]);
	printf("r1 %2.2f\n",P.r[1]);
	printf("r2 %2.2f\n",P.r[2]);
	printf("HW %2.2f\n",halfDimension);
	
	if(center[0]-halfDimension<=P.r[0])printf("tes tru");

	if(center[0]-halfDimension<=P.r[0] && P.r[0]<=center[0]+halfDimension && 
		center[1]-halfDimension<=P.r[1] && P.r[1]<=center[1]+halfDimension &&
		center[2]-halfDimension<=P.r[2] && P.r[2]<=center[2]+halfDimension)
	{
		printf("Truuu");
		return true;

	}
	
	printf("Falsee");
	return false;

}

OctTree::OctTree(vector<double> center_in,float halfDimension_in)
{	COM=nullptr;
	QT_NODE_CAPACITY=1;
	try{
		boundary=bounding_box(center_in,halfDimension_in);
		if(boundary.halfDimension<0.0)
		{
			throw(-1);
		}
	}
	catch(int num)
	{
		printf("Bounding box for octtree failed");
	}
}


bool OctTree::insert(nbd_object *P_in)
{
	printf("Insertion started\n");
	//obj doesnt belong to this tree
	printf("bb box size is %5.3f\n",boundary.halfDimension);
	printf("Cent %2.2f\n",boundary.center[0]);

	P_in->print_info();

	bool cond=boundary.center[0]-boundary.halfDimension<=P_in->r[0] && P_in->r[0]<=boundary.center[0]+boundary.halfDimension && 
		boundary.center[1]-boundary.halfDimension<=P_in->r[1] && P_in->r[1]<=boundary.center[1]+boundary.halfDimension &&
		boundary.center[2]-boundary.halfDimension<=P_in->r[2] && P_in->r[2]<=boundary.center[2]+boundary.halfDimension;
	if(!cond)
	{
	printf("Inside copy works");
	return false;

	}



	//printf("contain point inside insert %d",boundary.containsPoint(P_in));
/*	if(!boundary.containsPoint(P_in))
	{	printf("Object doesnt belong in boundary");
		return false;
	}
*/
	//if there is space and no subdivision

	if(node_stars.size()<QT_NODE_CAPACITY && oct_1==nullptr)
	{	printf("No subdivision needed");
		node_stars.push_back(P_in);
		return true;
	}

	if(oct_1==nullptr)
	{	printf("Subdiv needed");
		this->subdivide();
		
		for(int i=0; i<QT_NODE_CAPACITY;i++)
		{	nbd_object *temp_star=node_stars.back();
			node_stars.pop_back();
			if(COM==nullptr)
			{
				COM=new nbd_object(-1,temp_star->m,temp_star->r,temp_star->v);
				COM->r=scaling_vector(COM->r, COM->m);
			}
			else
			{
				COM->m=COM->m+temp_star->m;
				COM->r=elementwise_sum(COM->r,scaling_vector(temp_star->r, temp_star->m));
			}
			if(this->oct_1->insert(temp_star)) continue;	
			if(this->oct_2->insert(temp_star)) continue;	
			if(this->oct_3->insert(temp_star)) continue;	
			if(this->oct_4->insert(temp_star)) continue;	
			if(this->oct_5->insert(temp_star)) continue;	
			if(this->oct_6->insert(temp_star)) continue;	
			if(this->oct_7->insert(temp_star)) continue;	
			if(this->oct_8->insert(temp_star)) continue;	
		}
	}

	COM->m=COM->m+P_in->m;
	COM->r=elementwise_sum(COM->r,scaling_vector(P_in->r, P_in->m));

	if(this->oct_1->insert(P_in)) return true;	
	if(this->oct_2->insert(P_in)) return true;	
	if(this->oct_3->insert(P_in)) return true;	
	if(this->oct_4->insert(P_in)) return true;	
	if(this->oct_5->insert(P_in)) return true;	
	if(this->oct_6->insert(P_in)) return true;	
	if(this->oct_7->insert(P_in)) return true;	
	if(this->oct_8->insert(P_in)) return true;	
	printf("Insert Failed");
	return false;

}

void OctTree::subdivide()
{	printf("Subdivide Occuring\n");
	vector<double> c=boundary.center;
	float hw=boundary.halfDimension/2;
	this->oct_1= new OctTree(vector<double> {c[0]+hw,c[1]+hw,c[2]+hw},hw);
printf("Subdivide 1 done\n");

	this->oct_2= new OctTree(vector<double> {c[0]+hw,c[1]-hw,c[2]+hw},hw);

printf("Subdivide 2 done\n");
	this->oct_3= new OctTree(vector<double> {c[0]-hw,c[1]+hw,c[2]+hw},hw);

printf("Subdivide 3 done\n");
	this->oct_4= new OctTree(vector<double> {c[0]-hw,c[1]-hw,c[2]+hw},hw);

printf("Subdivide 4 done\n");
	this->oct_5= new OctTree(vector<double> {c[0]+hw,c[1]+hw,c[2]-hw},hw);

printf("Subdivide 5 done\n");
	this->oct_6= new OctTree(vector<double> {c[0]+hw,c[1]-hw,c[2]-hw},hw);

printf("Subdivide 6 done\n");
	this->oct_7= new OctTree(vector<double> {c[0]-hw,c[1]+hw,c[2]-hw},hw);

printf("Subdivide 7 done\n");
	this->oct_8= new OctTree(vector<double> {c[0]-hw,c[1]-hw,c[2]-hw},hw);

	printf("Subdivide 8 done\n");
	printf("Subdiv Done");
}

void OctTree::print_tree()
{
	printf("Printing this node\n");
	if(oct_1==nullptr) {
	printf("This node has %d points and no children nodes\n",node_stars.size());
	for(int i=0;i<node_stars.size();i++)
		node_stars[i]->print_info();
	}
	else {
	printf("This node has children nodes and the COM is %2.4f\n",COM->m);
	printf("oct_1\n");
	oct_1->print_tree();
	printf("oct_2\n");
	oct_2->print_tree();
	printf("oct_3\n");
	oct_3->print_tree();
	printf("oct_4\n");
	oct_4->print_tree();
	printf("oct_5\n");
	oct_5->print_tree();
	printf("oct_6\n");
	oct_6->print_tree();
	printf("oct_7\n");
	oct_7->print_tree();
	printf("oct_8\n");
	oct_8->print_tree();
	}

	printf("This node is over");

}


