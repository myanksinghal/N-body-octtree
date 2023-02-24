#ifndef OCTTREE_H
#define OCTTREE_H

#include "nbd_object.h"
#include <memory>
#include <vector>

using namespace std;

struct bounding_box
{

	vector<double> center;
	float halfDimension;
	bounding_box();
	bounding_box(vector<double> center_in,float halfDimension_in);
	bool containsPoint(nbd_object *P_in);
//	bool intersects_BB(bounding_box other_BB);

};


class OctTree
{
	public:	
		int QT_NODE_CAPACITY;
		bounding_box boundary;

		vector<nbd_object*> node_stars;

		nbd_object *COM=nullptr;

		
		//Children Nodes
		
		OctTree* oct_1=nullptr;		
		OctTree* oct_2=nullptr;
		OctTree* oct_3=nullptr;
		OctTree* oct_4=nullptr;

		OctTree* oct_5=nullptr;
		OctTree* oct_6=nullptr;
		OctTree* oct_7=nullptr;
		OctTree* oct_8=nullptr;

		OctTree(vector<double> center_in,float halfDimension_in);
		OctTree();
		~OctTree();
		bool insert(nbd_object *P_in);
		void subdivide();
		bool traverse(nbd_object *P_in);
		void print_tree();

};


#endif
