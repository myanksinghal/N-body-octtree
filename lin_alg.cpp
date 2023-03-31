#include "lin_alg.h"
#include <functional>
#include <algorithm>
#include <math.h>
#include <vector>
#include <cmath>
using namespace std;

/**
 * @brief Elementwise sum of two vectors with a scaling factor for the second vector
 *
 * @param V1 First vector
 * @param V2 Second vector
 * @param scale_factor Scaling factor for the second vector
 * @return vector<double> Sum of the two vectors with the second vector scaled by scale_factor
 */
vector<double> elementwise_sum(vector<double> const &V1, vector<double> const &V2, double scale_factor)
{
	vector<double> result;
	for (int i = 0; i < V1.size(); i++)
	{
		result.push_back(V1[i] + scale_factor * V2[i]);
	}
	return result;
}

/**
 * @brief Scaling of a vector by a scalar
 *
 * @param V1 Input vector
 * @param S1 Scaling factor
 * @return vector<double> Scaled vector
 */
vector<double> scaling_vector(vector<double> const &V1, const double S1)
{
	vector<double> result;
	for (int i = 0; i < V1.size(); i++)
	{
		result.push_back(V1[i] * S1);
	}

	return result;
}

/**
 * @brief Calculates the norm of a vector
 *
 * @param V1 Input vector
 * @return double Norm of the vector
 */
double norm(vector<double> const &V1)
{
	double result = 0.0;
	for (int i = 0; i < V1.size(); i++)
	{
		result += pow(V1[i], 2);
	}
	result = sqrt(result);

	return result;
}
/**
 * @brief Calculates the dot product of two vectors
 *
 * @param V1 First vector
 * @param V2 Second vector
 * @return double Dot product of the two vectors
 */
double dot_product(vector<double> const &V1, vector<double> const &V2)
{
	double result = 0.0;
	for (int i = 0; i < V1.size(); i++)
	{
		result += V1[i] * V2[i];
	}
	return result;
}
