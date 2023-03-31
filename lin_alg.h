#ifndef LIN_ALG_H
#define LIN_ALG_H

#include <vector>

using namespace std;
vector<double> elementwise_sum(vector<double> const &V1, vector<double> const &V2, double scale_factor);
vector<double> scaling_vector(vector<double> const &V1, double S1);
double norm(vector<double> const &V1);
double dot_product(vector<double> const &V1, vector<double> const &V2);
#endif
