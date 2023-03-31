#include "nbd_object.h"
#include "lin_alg.h"
#include <math.h>
#include <vector>
#include "consts.h"

/**
 * @brief Construct a new nbd object::nbd object object
 *
 * @param id_in id of the object
 * @param m_in mass of the object
 * @param r_in position of the object
 * @param v_in velocity of the object
 */
nbd_object::nbd_object(long id_in, float m_in, vector<double> r_in, vector<double> v_in)
{
	id = id_in;
	m = m_in;
	r = r_in;
	v = v_in;
	F_0 = {0.0, 0.0, 0.0};
	F_1 = {0.0, 0.0, 0.0};
	F_0_t0 = {0.0, 0.0, 0.0};
	F_1_t0 = {0.0, 0.0, 0.0};
	PE = 0.0;
	KE = 0.0;
	sugg_del_t = 0.1;
	t_block = 1;
	tblock_double_flip = true;
}

/**
 * @brief Calculates the force on the object due to another object
 *
 * @param p2 Second object
 */
void nbd_object::calculate_force(nbd_object *p2)
{
	// Softening the radius for close encounters
	vector<double> R_vec = elementwise_sum(r, p2->r, -1);
	double R_mag = norm(R_vec);
	double softened_R_Mag = sqrt(pow(SOFTENING_PARAM, 2) + pow(R_mag, 2));
	vector<double> V_vec = elementwise_sum(v, p2->v, -1);
	double V_mag = norm(V_vec);

	vector<double> F_0 = scaling_vector(R_vec, -p2->m / pow(softened_R_Mag, 3));
	double a = dot_product(R_vec, V_vec) / pow(softened_R_Mag, 2);

	vector<double> F_1 = scaling_vector(V_vec, -p2->m / pow(softened_R_Mag, 3));
	F_1 = elementwise_sum(F_1, F_0, -3 * a);

	this->F_0 = elementwise_sum(this->F_0, F_0, 1);
	this->F_1 = elementwise_sum(this->F_1, F_1, 1);

	this->PE += -G * this->m * p2->m / R_mag;
}

/**
 * @brief Advances the position and velocity of the object by one timestep
 *
 * @param del_t The value of the timestep
 * @param start_flag Flag to indicate if the timestep is the first timestep
 * @param in_current_block flag to indicate if the object is in the current block for 4th order correction
 */
void nbd_object::primary_time_advance(double del_t, bool *start_flag, bool *in_current_block)
{
	vector<double> F_0 = scaling_vector(this->F_0, G);
	vector<double> F_1 = scaling_vector(this->F_1, G);

	vector<double> v_new = elementwise_sum(scaling_vector(F_0, del_t), scaling_vector(F_1, 0.5 * del_t * del_t), 1);
	v_new = elementwise_sum(v_new, this->v, 1);
	vector<double> r_new = scaling_vector(F_1, del_t * del_t * del_t / 6);
	r_new = elementwise_sum(r_new, F_0, 0.5 * del_t * del_t);
	r_new = elementwise_sum(r_new, this->v, del_t);
	r_new = elementwise_sum(r_new, this->r, 1);

	this->r = r_new;
	this->v = v_new;
	if (!*start_flag && *in_current_block)
	{
		vector<double> F_0_t0 = scaling_vector(this->F_0_t0, G);
		vector<double> F_1_t0 = scaling_vector(this->F_1_t0, G);

		vector<double> temp0 = elementwise_sum(F_0_t0, F_0, -1.0);
		temp0 = scaling_vector(temp0, -6 / (del_t * del_t));
		vector<double> temp1 = elementwise_sum(F_1, F_1_t0, 2);
		temp1 = scaling_vector(temp1, -2 / (del_t));

		vector<double> F_2_t0 = elementwise_sum(temp0, temp1, 1);

		temp0 = elementwise_sum(F_0_t0, F_0, -1.0);
		temp0 = scaling_vector(temp0, 12 / (del_t * del_t * del_t));
		temp1 = elementwise_sum(F_1, F_1_t0, 1);
		temp1 = scaling_vector(temp1, 6 / (del_t * del_t));

		vector<double> F_3_t0 = elementwise_sum(temp0, temp1, 1);

		// printf("Force in advance timestep is %3.10f,%3.10f,%3.10f\n",F[0],F[1],F[2]);
		vector<double> v_corr = elementwise_sum(scaling_vector(F_2_t0, del_t * del_t * del_t / 6), scaling_vector(F_3_t0, del_t * del_t * del_t * del_t / 24), 1);
		this->v = elementwise_sum(v_corr, this->v, 1);
		vector<double> r_corr = elementwise_sum(scaling_vector(F_2_t0, del_t * del_t * del_t * del_t / 24), scaling_vector(F_3_t0, del_t * del_t * del_t * del_t * del_t / 120), 1);
		this->r = elementwise_sum(r_corr, this->r, 1);
		// printf("new_r=%3.3f,%3.3f,%3.3f\n old_r=%3.3f,%3.3f,%3.3f\n",r_new[0],r_new[1],r_new[2],r[0],r[1],r[2]);
		//
		this->sugg_del_t = sqrt((neta * (norm(F_0) * norm(F_2_t0) + pow(norm(F_1), 2))) / (norm(F_1) * norm(F_3_t0) + pow(norm(F_2_t0), 2)));

		// If the suggested timestep is smaller than the current and not already in the smallest block, decrease the block number
		if (this->sugg_del_t < del_t && this->t_block != 1)
		{
			this->t_block--;
		}
		// If the suggested timestep is larger than the next time block and not already in the largest block, increase the block number
		else if (this->tblock_double_flip && this->sugg_del_t > 2 * del_t && this->t_block < n_time_blocks - 1)
		{
			this->t_block++;
		}
		// Makes sure the block number is decreased if the number of blocks is decreased.
		if (this->t_block > n_time_blocks - 1)
		{
			this->t_block = n_time_blocks - 1;
		}

		// Makes sure the increase in block happens every other "current" timestep.
		this->tblock_double_flip = !this->tblock_double_flip;
	}
}

/**
 * @brief Prints the information of the object to the screen
 *
 */
void nbd_object::print_info()
{
	printf("Particle ID %7d of mass %3.1f is at %3.3f,%3.3f,%3.3f with vel %3.3f,%3.3f,%3.3f\n", id, m, r[0], r[1], r[2], v[0], v[1], v[2]);
}
