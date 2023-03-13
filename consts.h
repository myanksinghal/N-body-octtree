#ifndef CONSTS_H
#define CONSTS_H
#include <math.h>
const float G=1;
const double THETA_THRESHOLD=0.4;
const double SOFTENING_PARAM=0.01;
const int NODE_CAPACITY=1;
const double neta=0.01;
const double max_timestep=0.1;
const unsigned int n_time_blocks=12; 
const double t1=max_timestep/(pow(2,n_time_blocks));
#endif
