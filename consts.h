#ifndef CONSTS_H
#define CONSTS_H
#include <math.h>
#include <vector>
const float G = 1.0;
const double THETA_THRESHOLD = 0.4;
const double SOFTENING_PARAM = 0.0001;
const int NODE_CAPACITY = 1;
const double neta = 1;
const double max_timestep = 0.1;
extern bool increase_nstep;
extern bool decrease_nstep;
const unsigned int starting_n_time_blocks = 5;
const float DE_flag = 0.01;
extern unsigned int n_time_blocks;
extern double t1;

const bool external_body=true;
const float ext_mass = 1;
const std::vector<double> ext_r={0,0,0};
const std::vector<double> ext_v={0,0,0};

#endif
