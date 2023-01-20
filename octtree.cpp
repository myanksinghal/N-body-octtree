#include "octtree.h"
#include "lin_alg.h"
#include "nbd_object.h"
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
	if(center[0]-halfDimension<=P.r[0] && P.r[0]<=center[0]+halfDimension && 
		center[1]-halfDimension<=P.r[1] && P.r[1]<=center[1]+halfDimension &&
		center[2]-halfDimension<=P.r[2] && P.r[2]<=center[2]+halfDimension)
	{
		return true;

	}

	return false;

}

OctTree::OctTree(bounding_box boundary_in)
{	try{
		bounding_box boundary=boundary_in;
		if(boundary.halfDimension==-1)
		{
			throw(-1);
		}
	}
	catch(int num)
	{
		printf("Bounding box for octtree failed");
	}
}

bool OctTree::insert(nbd_object P_in)
{
	//obj doesnt belong to this tree
	if(!boundary.containsPoint(P_in))
		return false;

	//if there is space and no subdivision
	if(node_stars.size()<QT_NODE_CAPACITY && oct_1==nullptr)
	{
		node_stars.push_back(P_in);
		return true;
	}

	if(oct_1==nullptr)
	{
		subdivide();
		
		for(int i=0; i<QT_NODE_CAPACITY;i++)
		{	nbd_object temp_star=node_stars.back();
			node_stars.pop_back();
			if(COM==nullptr)
			{
				*COM=temp_star;
				COM->r=scaling_vector(COM->r, COM->m);
			}
			else
			{
				COM->m=COM->m+temp_star.m;
				COM->r=elementwise_sum(COM->r,scaling_vector(temp_star.r, temp_star.m));
			}
			if(oct_1->insert(temp_star)) continue;	
			if(oct_2->insert(temp_star)) continue;	
			if(oct_3->insert(temp_star)) continue;	
			if(oct_4->insert(temp_star)) continue;	
			if(oct_5->insert(temp_star)) continue;	
			if(oct_6->insert(temp_star)) continue;	
			if(oct_7->insert(temp_star)) continue;	
			if(oct_8->insert(temp_star)) continue;	
		}
	}

	COM->m=COM->m+P_in.m;
	COM->r=elementwise_sum(COM->r,scaling_vector(P_in.r, P_in.m));

	if(oct_1->insert(P_in)) return true;	
	if(oct_2->insert(P_in)) return true;	
	if(oct_3->insert(P_in)) return true;	
	if(oct_4->insert(P_in)) return true;	
	if(oct_5->insert(P_in)) return true;	
	if(oct_6->insert(P_in)) return true;	
	if(oct_7->insert(P_in)) return true;	
	if(oct_8->insert(P_in)) return true;	

	return false;

}
bool OctTree::subdivide()
{

return false;
}
//vector<point> query();


