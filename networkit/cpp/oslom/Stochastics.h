#pragma once

#ifndef STOCHASTICS_HPP
#define STOCHASTICS_HPP

#include "utils/OslomRandom.h"
#include "LogTable.h"
#include "DataTypes.h"
#include "Parameters.h"

# define sqrt_two 1.41421356237
# define bisection_precision 1e-2

class Stochastics {

public:

	Stochastics() = delete;

	static void init();

	/**
	 * r1 \in (a,b)
	 * r2 \in (c,d)
	 * compute the probability p(r2<r1)
	 */
	static double compare_r_variables(double a, double b, double c, double d);

	static double right_error_function(double x);

	static double log_together(double minus_log_total, int number);

	static double fitted_exponent(int N);

	static inline double order_statistics_left_cumulative(int N, int pos, double x) {
		// the routine computes the probality c_pos = p(X_pos <= x)
		// N is the total number of variables, pos is from 1 to N. N is the smallest.
		return log_table->cum_binomial_right(N - pos + 1, N, x);
	}

	static double inverse_order_statistics(int sample_dim, int pos, const double &zerof, double lo,
	                                       double hi);

	static inline double pron_min_exp(int N, double xi) {
		// this should return the probability that the minimum of the quantiles p(c_min<=xi)
		//cout<<"-> "<<fitted_exponent(N)<<std::endl;
		return 1 - exp(-fitted_exponent(N) * xi);
	}

	static inline double
	compute_probability_to_stop(const double &a, const double &b, const double &critical_xi,
	                            int Nstar, int pos) {
		/*	this routine is to compute the bootstrap probaility that the node with extremes a and b will be below threshold
			when I already know that the average is below critical_xi */
		if (order_statistics_left_cumulative(Nstar, pos, b) <= critical_xi)
			return 1;

		return (inverse_order_statistics(Nstar, pos, critical_xi, (a + b) * 0.5, b) - a) /
		       (b - a);
	}

	static bool equivalent_check(int pos_first, int pos_last, double &A_average, double &B_average,
	                             int equivalents, int Nstar, const double &critical_xi);

	static bool
	equivalent_check_gather(CupDataStruct &a, int &add_nodes, const double &probability_a,
	                        const double &probability_b, int Nstar,
	                        const double &critical_xi);

	static inline double hyper_table(int kin_node, int kout_g, int tm, int degree_node) {
		return log_table->hyper(kin_node, kout_g, tm, degree_node);
	}

	/**
	 * Calculate the cumulative probability that the node has kin_node or more edges into the
	 * module in the null model (?) = r-Score ?
	 * @param kin_node Number of edges between the node and the module.
	 * @param kout_g Number of outgoing stubs of the module.
	 * @param tm Number of stubs outside of the group.
	 * @param degree_node Degree of the node.
	 * @return
	 */
	static inline double topological_05(int kin_node, int kout_g, int tm, int degree_node) {
		return log_table->right_cumulative_function(degree_node, kout_g, tm, kin_node + 1) +
		       ran4() * hyper_table(kin_node, kout_g, tm, degree_node);
	}

	// TODO: Does this function compute the r-Score?
	static double compute_global_fitness(int kin_node, int kout_g, int tm, int degree_node,
	                                     double minus_log_total, int number_of_neighs, int Nstar,
	                                     double &boot_interval);

	static double compute_global_fitness_step(int kin_node, int kout_g, int tm, int degree_node,
	                                          double minus_log_total, int number_of_neighs,
	                                          int Nstar, double _step_);

	static inline double
	compute_global_fitness_ofive(int kin_node, int kout_g, int tm, int degree_node,
	                             double minus_log_total, int number_of_neighs, int Nstar) {
		return compute_global_fitness_step(kin_node, kout_g, tm, degree_node, minus_log_total,
		                                   number_of_neighs, Nstar, 0.5);
	}


	static inline double
	compute_global_fitness_randomized(int kin_node, int kout_g, int tm, int degree_node,
	                                  double minus_log_total, int number_of_neighs,
	                                  int Nstar) {
		return compute_global_fitness_step(kin_node, kout_g, tm, degree_node, minus_log_total,
		                                   number_of_neighs, Nstar, ran4());
	}

	static double
	compute_global_fitness_randomized_short(int kin_node, int kout_g, int tm, int degree_node,
	                                        double minus_log_total);

private:
	static LogFactTable *log_table;
	static Parameters *paras;
};

#endif //STOCHASTICS_HPP
