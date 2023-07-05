#ifndef CONSTS_H
#define CONSTS_H
#include <math.h>
#include <vector>
const float G = 1.0;

//G=1, using Dist=1pc, total mass M = 1 Msun then V=0.06558km/s

const double THETA_THRESHOLD = 0.4;
const double SOFTENING_PARAM = 0.1;
const float alpha=-2.35; //Power for powerlaw profile
const int NODE_CAPACITY = 1;
const double neta = 0.01;
const double max_timestep = 0.01;
const double OUT_TIMESTEP =0.001;
extern bool increase_nstep;
extern bool decrease_nstep;
const unsigned int starting_n_time_blocks = 7;
const unsigned int max_t_blocks=20;
const float DE_flag = 0.001;
extern unsigned int n_time_blocks;
extern double t1;

const bool external_body=true;
const float ext_mass = 1.0;
const float ext_halfmass_rad=5.0;
const std::vector<double> ext_r={0,0,0};
const std::vector<double> ext_v={0,0,0};

#endif
