#ifndef CONSTS_H
#define CONSTS_H
#include <math.h>
const float G=1;
const double THETA_THRESHOLD=0.4;
const double SOFTENING_PARAM=0.0001;
const int NODE_CAPACITY=1;
const double neta=1;
const double max_timestep=0.1;
extern bool increase_nstep; 
extern bool decrease_nstep;
const unsigned int starting_n_time_blocks=5;
const float DE_flag=0.001;
extern unsigned int n_time_blocks; 
extern double t1;
#endif
