#include "octtree.h"
#include "lin_alg.h"
#include "nbd_object.h"
#include "consts.h"
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

bool bounding_box::containsPoint(nbd_object *P_in)
{


	bool cond=this->center[0]-this->halfDimension<=P_in->r[0] && P_in->r[0]<=this->center[0]+this->halfDimension && 
		this->center[1]-this->halfDimension<=P_in->r[1] && P_in->r[1]<=this->center[1]+this->halfDimension &&
		this->center[2]-this->halfDimension<=P_in->r[2] && P_in->r[2]<=this->center[2]+this->halfDimension;
	if(cond)
	{
	return true;
	}
	else 
	{
	return false;
	}


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
	//obj doesnt belong to this tree

//	P_in->print_info();




	//printf("contain point inside insert %d",boundary.containsPoint(P_in));
	if(!boundary.containsPoint(P_in))
	{//	printf("Object doesnt belong in boundary\n");
		return false;
	}

	//if there is space and no subdivision

	if(node_stars.size()<QT_NODE_CAPACITY && oct_1==nullptr)
	{	
		node_stars.push_back(P_in);
		this->COM=new nbd_object(-1,P_in->m,P_in->r,P_in->v);
		this->COM->r=scaling_vector(this->COM->r, this->COM->m);

		return true;
	}

	if(oct_1==nullptr)
	{	
		this->subdivide();
		delete this->COM;
		this->COM=nullptr;
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
{	
	vector<double> c=boundary.center;
	float hw=boundary.halfDimension/2;
	this->oct_1= new OctTree(vector<double> {c[0]+hw,c[1]+hw,c[2]+hw},hw);

	this->oct_2= new OctTree(vector<double> {c[0]+hw,c[1]-hw,c[2]+hw},hw);

	this->oct_3= new OctTree(vector<double> {c[0]-hw,c[1]+hw,c[2]+hw},hw);

	this->oct_4= new OctTree(vector<double> {c[0]-hw,c[1]-hw,c[2]+hw},hw);

	this->oct_5= new OctTree(vector<double> {c[0]+hw,c[1]+hw,c[2]-hw},hw);

	this->oct_6= new OctTree(vector<double> {c[0]+hw,c[1]-hw,c[2]-hw},hw);

	this->oct_7= new OctTree(vector<double> {c[0]-hw,c[1]+hw,c[2]-hw},hw);

	this->oct_8= new OctTree(vector<double> {c[0]-hw,c[1]-hw,c[2]-hw},hw);

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

OctTree::~OctTree()
{
	delete this->COM;
	if(this->oct_1!=nullptr)
	{
		delete this->oct_1;
		delete this->oct_2;
		delete this->oct_3;
		delete this->oct_4;
		delete this->oct_5;
		delete this->oct_6;
		delete this->oct_7;
		delete this->oct_8;
	}
	
}

bool OctTree::traverse(nbd_object *P_in)
{
	if(this->node_stars.size()==0 && this->oct_1==nullptr)
	{
		return true;
	}
	float s=this->boundary.halfDimension;
	vector<double> corrected_r_COM=scaling_vector(this->COM->r,1/this->COM->m);
	double d=norm(elementwise_sum(P_in->r,scaling_vector(corrected_r_COM,-1)));
	if(s/d<THETA_THRESHOLD)
	{
		//force evaluation
		nbd_object *temporary_obj= new nbd_object(-1,this->COM->m,corrected_r_COM,vector<double>{-1,-1,-1});
		P_in->calculate_force(temporary_obj);

	//	printf("Force calculated from COM is %3.3f,%3.3f,%3.3f\n",F_0[0],F_0[1],F_0[2]);
		delete temporary_obj;	
		return true;
	}
	else if (this-> oct_1==nullptr) 
	{
		//Going through each particle in the box

		for(int i=0; i<QT_NODE_CAPACITY;i++)
		{
			nbd_object *temp_star=node_stars.back();
			if(temp_star->id!=P_in->id)
			{
				//force eval
				P_in->calculate_force(temp_star);

			//	printf("Force calculated from ind is %3.3f,%3.3f,%3.3f\n",F_0[0],F_0[1],F_0[2]);
				
			//	printf("Force before update %3.3f,%3.3f,%3.3f\n",P_in->F[0],P_in->F[1],P_in->F[2]);

			//	printf("Force after update %3.3f,%3.3f,%3.3f\n",P_in->F[0],P_in->F[1],P_in->F[2]);
			}

		}

		return true;

	}
	else{
			this->oct_1->traverse(P_in);	
			this->oct_2->traverse(P_in);	
			this->oct_3->traverse(P_in);	
			this->oct_4->traverse(P_in);	
			this->oct_5->traverse(P_in);	
			this->oct_6->traverse(P_in);	
			this->oct_7->traverse(P_in);	
			this->oct_8->traverse(P_in);	
			return true;
	}


}

