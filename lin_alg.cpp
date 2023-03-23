#include "lin_alg.h"
#include <functional>
#include <algorithm>
#include <math.h>
#include <vector>
#include <cmath>
using namespace std;

vector<double> elementwise_sum(vector<double> const &V1, vector<double> const &V2,double scale_factor)
{
	vector<double> result;
	for(int i=0;i<V1.size();i++)
	{
		result.push_back(V1[i]+scale_factor*V2[i]);
	}
	return result;
}


vector<double> scaling_vector(vector<double> const &V1, const double S1)
{
	vector<double> result;
	for(int i=0;i<V1.size();i++)
	{
		result.push_back(V1[i]*S1);
	}

	return result;

}


double norm(vector<double> const &V1)
{
	double result=0.0;
	for(int i=0; i<V1.size();i++)
	{
		result+=pow(V1[i],2);
	}
	result=sqrt(result);

	return result;


}

double dot_product(vector<double> const &V1, vector<double> const &V2)
{
	double result=0.0;
	for(int i=0; i<V1.size();i++)
	{
		result+=V1[i]*V2[i];
	}
	return result;


}
