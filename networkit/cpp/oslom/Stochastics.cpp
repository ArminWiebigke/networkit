#include <cassert>

#include "Stochastics.h"


void Stochastics::init(int num_edges) {
	log_table = LogFactTable::get_instance();
	log_table->set(num_edges * 2);
	paras = Parameters::get_instance();
}

double Stochastics::compare_r_variables(double a, double b, double c, double d) {
	//
	if (c < a)
		return (1 - compare_r_variables(c, d, a, b));
	if (c > b)
		return 0;
	if (d > b)
		return 0.5 * (b - c) * (b - c) / ((d - c) * (b - a));
	else
		return (b - 0.5 * (d + c)) / (b - a);
}

double Stochastics::right_error_function(double x) {
	return 0.5 * std::erfc(x / sqrt_two);
}

double Stochastics::log_together(double minus_log_total, int number) {
	if (number < 11) {
		double fa = 1;
		double zum = 1;
		for (int i = 1; i < number; i++) {
			fa *= minus_log_total / i;
			zum += fa;
		}
		return std::max(zum * exp(-minus_log_total), 1e-100);
	}
	double mu = number;
	return std::max(right_error_function((minus_log_total - mu) / sqrt(mu)), 1e-100);
}

double Stochastics::fitted_exponent(int N) {
	double l = log(double(N));

	if (N > 100)
		return 4.2 * l - 8.5;

	if (N > 30)
		return 3.5 * l - 5.5;

	if (N > 7)
		return 2.5 * l - 2;

	if (N > 1)
		return 1.3 * l + 0.1;

	return 1;
}

double
Stochastics::inverse_order_statistics(int sample_dim, int pos, const double &zerof, double lo,
                                      double hi) {
	// anyway it's possible that I will tabulated this stuff for N small like sample_dim < 10000 which can already help a lot.
	// moreover, i should not use bisection but something better
	// finally there is the gaussian approx which could be studied more...
	//file: /Users/admin/Desktop/ambiente/right_binomial_cum/approx.cpp

	//double zerof= - log(1-threshold)/fitted_exponent(sample_dim);
	double mid = (hi + lo) / 2;

	while ((mid != lo) && (mid != hi)) {
		double fmid = order_statistics_left_cumulative(sample_dim, pos, mid);

		if (fabs(fmid - zerof) < bisection_precision * zerof)
			break;

		if ((fmid - zerof) <= 0)
			lo = mid;
		else
			hi = mid;

		mid = (hi + lo) / 2;
	}
	return mid;
}

bool
Stochastics::equivalent_check(int pos_first, int pos_last, double &A_average, double &B_average,
                              int equivalents, int Nstar, const double &critical_xi) {
	// returns true is the test was passed
	/*small_simulation(pos_first, pos_last, A_average, B_average, equivalents, Nstar, critical_xi);
	cout<<pos_first<<" "<<pos_last<<" A, B "<<A_average<<" "<<B_average<<" "<<equivalents<<" "<<Nstar<<" "<<critical_xi<<std::endl;*/

	int pos = pos_first;
	double cr_previous = A_average;

	for (int i = equivalents; i >= 1; i--) {
		// loop which starts from the best node
		//cout<<"i... "<<i<<std::endl;
		if (order_statistics_left_cumulative(Nstar, pos, cr_previous) <= critical_xi) {
			if (order_statistics_left_cumulative(Nstar, pos, B_average) <= critical_xi)
				return true;

			double cr = inverse_order_statistics(Nstar, pos, critical_xi, cr_previous,
			                                     B_average);
			//cout<<i<<" cr: "<<cr<<" "<<order_statistics_left_cumulative(equivalents, i, (cr-A_average)/(B_average-A_average))<<std::endl;
			//cout<<" this,  "<<order_statistics_left_cumulative(Nstar, pos, cr)<<" "<<Nstar<<" "<<pos<<std::endl;
			if (order_statistics_left_cumulative(equivalents, i,
			                                     (cr - A_average) /
			                                     (B_average - A_average)) > 0.5)
				return true;
			cr_previous = cr;
		}

		--pos;
	}
	return false;
}

bool
Stochastics::equivalent_check_gather(CupDataStruct &a, int &add_nodes, const double &probability_a,
                                     const double &probability_b, int Nstar,
                                     const double &critical_xi) {
	int nodes_added = 0;
	double A_average = 0;
	double B_average = 0;
	int pos_first = -1;
	int pos_last = -1;

	auto itl = a.begin();
	for (; itl != a.end(); ++itl) {
		if (nodes_added == add_nodes)
			break;

		if (compare_r_variables(probability_a, probability_b,
		                        itl->first - itl->second.second,
		                        itl->first + itl->second.second) >
		    paras->equivalence_parameter) {

			if (pos_first == -1)
				pos_first = Nstar - nodes_added;

			A_average += itl->first - itl->second.second;
			B_average += itl->first + itl->second.second;
			pos_last = Nstar - nodes_added;
			//cout<<"pos_first: "<<pos_first<<" ... "<<pos_last<<std::endl;
		}
		++nodes_added;
	}

	int equivalents = pos_first - pos_last + 1;
	A_average /= equivalents;
	B_average /= equivalents;

	if (equivalents == 1)
		return true;

	if (equivalent_check(pos_first, pos_last, A_average, B_average, equivalents, Nstar,
	                     critical_xi)) {
		//cout<<"check passed"<<std::endl;
		return true;
	} else {
		add_nodes = 0;
		//cout<<"check not passed"<<std::endl;
		return false;
	}
}

/**
 * Compute the probability that a node from the null model has k_in or more edges into the group.
 * @param k_in Edges from node to group
 * @param gr_out Outgoing stubs of the group
 * @param open_stubs Stubs in the graph, excluding the node and the group
 * @param k_degree Degree of the node
 * @param boot_interval
 * @return
 */
double Stochastics::compute_simple_fitness(int k_in, int gr_out, int open_stubs, int k_degree) {
	double topologic = log_table->right_cumulative_function(k_degree, gr_out, open_stubs, k_in);

	if (paras->weighted)
		throw std::runtime_error("Don't use weighted graphs!");
	if (topologic > 1)
		topologic = 1;
	return std::max(topologic, 1e-100);
}

/**
 *
 * @param k_in Edges from node to group.
 * @param gr_out Outgoing stubs of the group.
 * @param tm External stubs. Stubs from (G - {group and node})
 * @param k_degree Degree of the node
 * @param minus_log_total ???
 * @param number_of_neighs ???
 * @param Nstar Number of external nodes (#G - #group - 1)
 * @param boot_interval (out): Half the domain of the fitness (?)
 * @return
 */
double Stochastics::compute_global_fitness(int k_in, int gr_out, int tm, int k_degree,
                                           double minus_log_total, int number_of_neighs, int Nstar,
                                           double &boot_interval) {
	/* k_in is referred to the node and not to the module */
	/* boot_interval is half the domain of the fitness. We assume that also in the weighted case the fitness can be linearized */
	/* which is true if the boot_interval is small enough	*/
	boot_interval = (0.5 + 1e-6 * (ran4() - 0.5)) *
			hypergeom_dist(k_in, gr_out, tm, k_degree);
	double topologic =
			log_table->right_cumulative_function(k_degree, gr_out, tm, k_in + 1) +
			boot_interval;

	//cout<<"----------> "<<log_table->right_cumulative_function(k_degree, gr_out, tm, k_in+1)<<"  boot_interval "<<boot_interval<<" k_in: "<<k_in<<" / "<<k_degree<<std::endl;
	//cout<<"<> "<<log_table->right_cumulative_function(k_degree, gr_out, tm, k_in)<<"  boot_interval "<<boot_interval<<" k_in: "<<k_in<<" / "<<k_degree<<std::endl;
	if (paras->weighted)
		throw std::runtime_error("Don't use weighted graphs!");
	if (topologic > 1)
		topologic = 1;
	boot_interval = std::min(boot_interval, 1. - topologic);
	boot_interval = std::min(boot_interval, topologic);
	return std::max(topologic, 1e-100);
}

double Stochastics::compute_global_fitness_step(int k_in, int gr_out, int tm, int k_degree,
                                                double minus_log_total, int number_of_neighs,
                                                int Nstar, double _step_) {
	double topologic =
			log_table->right_cumulative_function(k_degree, gr_out, tm, k_in + 1) +
			_step_ * (hypergeom_dist(k_in, gr_out, tm, k_degree));
	if (paras->weighted)
		throw std::runtime_error("Don't use weighted graphs!");
	return std::max(topologic, 1e-100);
}

double Stochastics::compute_global_fitness_randomized_short(int k_in, int gr_out, int tm,
                                                            int k_degree,
                                                            double minus_log_total) {
	// this function is used in try_to_assign_homeless.
	// the usual problem is that we don't know the number of neighbors of the module.
	// this could be taken in account with some more thinking...

	double b2 = log_table->right_cumulative_function(k_degree, gr_out, tm,
	                                                 k_in + 1);

	double topologic = b2 + 0.5 * (hypergeom_dist(k_in, gr_out, tm, k_degree));

	if (!paras->weighted)
		return std::max(topologic, 1e-100);

	double weight_part = log_together(minus_log_total, k_in);

	if (topologic <= 1e-100 || weight_part <= 1e-100)
		return 1e-100;

	return std::min(topologic, weight_part);
}

double Stochastics::calc_score(int node_degree, int k_in, int gr_out, int ext_stubs, int ext_nodes,
                               int position) {
	assert(gr_out <= ext_stubs);
	double boot_interval;
	double r_score = Stochastics::compute_simple_fitness(k_in, gr_out, ext_stubs,
	                                                     node_degree);
	double ordered_stat = Stochastics::order_statistics_left_cumulative(ext_nodes,
	                                                                    ext_nodes - position,
	                                                                    r_score);
	return ordered_stat;
}


LogFactTable *Stochastics::log_table;
Parameters *Stochastics::paras;