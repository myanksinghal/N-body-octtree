#ifndef OCTTREE_H
#define OCTTREE_H

#include "nbd_object.h"
#include <memory>
#include <vector>

using namespace std;


struct point
{
	vector<double> pos;
	long star_ind;

	point(vector<double> pos_in, long star_in);

};

struct bounding_box
{

	vector<double> center;
	float halfDimension;
	bounding_box();
	bounding_box(vector<double> center_in,float halfDimension_in);
	bool containsPoint(nbd_object P);
//	bool intersects_BB(bounding_box other_BB);

};


class OctTree
{
	public:	
		const int QT_NODE_CAPACITY=1;
		bounding_box boundary;

		vector<nbd_object> node_stars;

		nbd_object *COM=nullptr;

		//Children Nodes
		
		OctTree* oct_1;
		OctTree* oct_2;
		OctTree* oct_3;
		OctTree* oct_4;

		OctTree* oct_5;
		OctTree* oct_6;
		OctTree* oct_7;
		OctTree* oct_8;

		OctTree(bounding_box boundary_in);
		bool insert(nbd_object P_in);
		bool subdivide();
		vector<point> query();

};


#endif
