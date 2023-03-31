#include "octtree.h"
#include "lin_alg.h"
#include "nbd_object.h"
#include "consts.h"
#include <cstdio>
#include <memory>
#include <vector>

using namespace std;

/**
 * @brief Construct a new bounding box with known parameters
 *
 * @param center_in The center of the bounding box
 * @param halfDimension_in The half length of the bounding box
 */
bounding_box::bounding_box(vector<double> center_in, float halfDimension_in)
{
	center = center_in;
	halfDimension = halfDimension_in;
}
/**
 * @brief Construct a new bounding box with temporary parameters
 *
 */
bounding_box::bounding_box()
{
	vector<double> tem = {-1, -1, -1};
	center = tem;
	halfDimension = -1;
}

/**
 * @brief Checks whether a point is within the bounding box
 *
 * @param P_in Pointer to the n body object to check
 * @return true Point is in the bounding box
 * @return false Point is outside the bounding box
 */
bool bounding_box::containsPoint(nbd_object *P_in)
{

	bool cond = this->center[0] - this->halfDimension <= P_in->r[0] && P_in->r[0] <= this->center[0] + this->halfDimension &&
				this->center[1] - this->halfDimension <= P_in->r[1] && P_in->r[1] <= this->center[1] + this->halfDimension &&
				this->center[2] - this->halfDimension <= P_in->r[2] && P_in->r[2] <= this->center[2] + this->halfDimension;
	if (cond)
	{
		return true;
	}
	else
	{
		return false;
	}
}
/**
 * @brief Construct a new Oct Tree object
 *
 * @param center_in Position of the center of the OctTree root bounding box
 * @param halfDimension_in Half length of the root bounding box
 */
OctTree::OctTree(vector<double> center_in, float halfDimension_in)
{
	COM = nullptr;
	this->OT_NODE_CAPACITY = NODE_CAPACITY;
	// Checks weather the bounding box created is valid or not
	try
	{
		boundary = bounding_box(center_in, halfDimension_in);
		if (boundary.halfDimension < 0.0)
		{
			throw(-1);
		}
	}
	catch (int num)
	{
		printf("Bounding box for octtree failed");
	}
}

/**
 * @brief Destroy the Oct Tree object
 *
 */
OctTree::~OctTree()
{
	delete this->COM;
	if (this->oct_1 != nullptr)
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

/**
 * @brief Insert a point into the OctTree
 *
 * @details Inserts the point in the OctTree if it belongs in the bounding box, if the bounding box is full then it subdivides the tree, and repopulates the new branches.
 *
 * @param P_in The point to be inserted
 * @return true Insertion succeeded
 * @return false Insertion failed
 */
bool OctTree::insert(nbd_object *P_in)
{

	// If the point does not belong in this tree, exit
	if (!boundary.containsPoint(P_in))
	{ //	printf("Object doesnt belong in boundary\n");
		return false;
	}

	// If there is space and no subdivision of nodes
	if (node_stars.size() < this->OT_NODE_CAPACITY && oct_1 == nullptr)
	{
		node_stars.push_back(P_in);
		this->COM = new nbd_object(-1, P_in->m, P_in->r, P_in->v);
		this->COM->r = scaling_vector(this->COM->r, this->COM->m);
		this->COM->r = scaling_vector(this->COM->r, this->COM->m);
		return true;
	}
	// If there is no space and no subdivision of nodes, it subdivides
	if (oct_1 == nullptr)
	{
		// subdivides node
		this->subdivide();
		delete this->COM;
		this->COM = nullptr;
		// Go through each star existing in the node already and then add its mass and radius to COM vector and then insert it into the subdivided nodes
		for (int i = 0; i < this->OT_NODE_CAPACITY; i++)
		{
			nbd_object *temp_star = node_stars.back();
			node_stars.pop_back();
			if (COM == nullptr)
			{
				COM = new nbd_object(-1, temp_star->m, temp_star->r, temp_star->v);
				COM->r = scaling_vector(COM->r, COM->m);
				COM->v = scaling_vector(COM->v, COM->m);
			}
			else
			{
				COM->m = COM->m + temp_star->m;
				COM->r = elementwise_sum(COM->r, temp_star->r, temp_star->m);
				COM->r = elementwise_sum(COM->r, temp_star->r, temp_star->m);
			}
			if (this->oct_1->insert(temp_star))
				continue;
			if (this->oct_2->insert(temp_star))
				continue;
			if (this->oct_3->insert(temp_star))
				continue;
			if (this->oct_4->insert(temp_star))
				continue;
			if (this->oct_5->insert(temp_star))
				continue;
			if (this->oct_6->insert(temp_star))
				continue;
			if (this->oct_7->insert(temp_star))
				continue;
			if (this->oct_8->insert(temp_star))
				continue;
		}
	}
	// Add the new inserting point to the COM values for the node
	COM->m = COM->m + P_in->m;
	COM->r = elementwise_sum(COM->r, P_in->r, P_in->m);
	COM->v = elementwise_sum(COM->v, P_in->v, P_in->m);
	if (this->oct_1->insert(P_in))
		return true;
	if (this->oct_2->insert(P_in))
		return true;
	if (this->oct_3->insert(P_in))
		return true;
	if (this->oct_4->insert(P_in))
		return true;
	if (this->oct_5->insert(P_in))
		return true;
	if (this->oct_6->insert(P_in))
		return true;
	if (this->oct_7->insert(P_in))
		return true;
	if (this->oct_8->insert(P_in))
		return true;
	printf("Insert Failed");
	return false;
}
/**
 * @brief Subdivides the node into 8 equal sized cubes
 */
void OctTree::subdivide()
{
	vector<double> c = boundary.center;
	float hw = boundary.halfDimension / 2;
	this->oct_1 = new OctTree(vector<double>{c[0] + hw, c[1] + hw, c[2] + hw}, hw);

	this->oct_2 = new OctTree(vector<double>{c[0] + hw, c[1] - hw, c[2] + hw}, hw);

	this->oct_3 = new OctTree(vector<double>{c[0] - hw, c[1] + hw, c[2] + hw}, hw);

	this->oct_4 = new OctTree(vector<double>{c[0] - hw, c[1] - hw, c[2] + hw}, hw);

	this->oct_5 = new OctTree(vector<double>{c[0] + hw, c[1] + hw, c[2] - hw}, hw);

	this->oct_6 = new OctTree(vector<double>{c[0] + hw, c[1] - hw, c[2] - hw}, hw);

	this->oct_7 = new OctTree(vector<double>{c[0] - hw, c[1] + hw, c[2] - hw}, hw);

	this->oct_8 = new OctTree(vector<double>{c[0] - hw, c[1] - hw, c[2] - hw}, hw);
}

/**
 * @brief Prints the tree
 *
 */
void OctTree::print_tree()
{
	printf("Printing this node\n");
	if (oct_1 == nullptr)
	{
		printf("This node has %d points and no children nodes\n", node_stars.size());
		for (int i = 0; i < node_stars.size(); i++)
			node_stars[i]->print_info();
	}
	else
	{
		printf("This node has children nodes and the COM is %2.4f\n", COM->m);
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
/**
 * @brief Traverses through the tree through the perspective of input star and its force calculations
 *
 * @param P_in The input star for which we calculate the force
 * @return true If force calculation was successful
 * @return false If force calculation was not successful
 */
bool OctTree::traverse(nbd_object *P_in)
{
	// skip if node has no particles and no subnodes
	if (this->node_stars.size() == 0 && this->oct_1 == nullptr)
	{
		return true;
	}

	// Calculates if node calculation does not need to be precise
	float s = this->boundary.halfDimension;
	vector<double> corrected_r_COM = scaling_vector(this->COM->r, 1 / this->COM->m);
	vector<double> corrected_v_COM = scaling_vector(this->COM->v, 1 / this->COM->m);
	double d = norm(elementwise_sum(P_in->r, corrected_r_COM, -1));
	if (s / d < THETA_THRESHOLD)
	{
		// Creation of temporary object with the COM params and force calculation
		nbd_object *temporary_obj = new nbd_object(-1, this->COM->m, corrected_r_COM, corrected_v_COM);
		P_in->calculate_force(temporary_obj);

		delete temporary_obj;
		return true;
	}
	// If there is no subnodes then go through each particle
	else if (this->oct_1 == nullptr)
	{
		// Going through each particle in the box

		for (int i = 0; i < this->OT_NODE_CAPACITY; i++)
		{
			nbd_object *temp_star = node_stars.back();
			if (temp_star->id != P_in->id)
			{
				P_in->calculate_force(temp_star);
			}
		}

		return true;
	}
	// Else go through each subnode
	else
	{
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
	// Should never reach this false statement
	printf("\nERROR: ISSUE WITH TRAVERSAL OF OCT TREE\n");
	return false;
}
