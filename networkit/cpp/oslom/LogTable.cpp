#include <iostream>
#include <cmath>

#include "LogTable.h"

LogFactTable *LogFactTable::instance = nullptr;

LogFactTable *LogFactTable::get_instance() {
	if (instance == nullptr)
		instance = new LogFactTable();
	return instance;
}

void LogFactTable::set(int size) {
	if (size < lnf.size())
		return;

	std::cout << "allocating " << size << " factorials..." << std::endl;
	lnf.reserve(size + 1);

	double f = 0.0;
	size_t old_size = lnf.size();

	if (old_size > 0) {
		f = lnf.back();
	} else {
		lnf.push_back(0.0);
		old_size = 1;
	}

	for (int i = old_size; i <= size; i++) {
		f += std::log(i);
		lnf.push_back(f);
	}
	std::cout << "done" << std::endl;
	//prints(lnf);
}

double
LogFactTable::cum_hyper_right(int k_in, int gr_out, int open_stubs, int k_degree) {
	// TODO
	//std::cout<<"k_in... "<<k_in<<" "<<gr_out<<" "<<open_stubs<<" "<<k_degree<<std::endl;
	// this is bigger or equal p(x >= k_in)   *** EQUAL ***
	if (k_in > std::min(k_degree, gr_out))
		throw std::runtime_error("k_in too large!");

#define W(x) #x << "=" << x << ", "
	if (open_stubs - gr_out - k_degree + k_in <= 0) {
		std::cout << W(open_stubs) << W(gr_out) << W(k_degree) << W(k_in) << std::endl;
		throw std::runtime_error("open_stubs too small!");
	}

	if (k_in <= 0)
		return 1;

	// k_in is smaller than the expected value for a random node
	if (k_in < double(gr_out + 1) / double(open_stubs + 2) * double(k_degree + 1)) {
		return (1. - cum_hyper_left(k_in, gr_out, open_stubs, k_degree));
	}

	int x = k_in;
	double pzero = hypergeom_dist(x, gr_out, open_stubs, k_degree);

	if (pzero <= 1e-40)
		return 0;

	double ga = open_stubs - gr_out - k_degree; // TODO: open_stubs with or without k_degree?
	int kout_g_p = gr_out + 1;
	double degree_node_p = k_degree + 1;
	double z_zero = 1.;
	double sum = z_zero;

	while (true) {
		++x;
		z_zero *= double(kout_g_p - x) / (x * (ga + x)) * (degree_node_p - x);
		if (z_zero < log_table_pr * sum)
			break;
		if (pzero * sum > 1)
			return pzero;
		sum += z_zero;
	}

	return pzero * sum;
}

//double
//LogFactTable::cum_hyper_right2(int k_in, int gr_out, int ext_stubs, int k_degree) {
//	// k_in is smaller than the expected value for a random node, so it is for sure (?) not
//	// significant
//	if (k_in < gr_out / static_cast<double>(ext_stubs) * (k_degree))
//		return 1.0;
//
//	int x = k_in;
//	double pzero = hypergeom_dist(x, gr_out, ext_stubs, k_degree);
//
//	if (pzero <= 1e-40)
//		return 0;
//
//	double z_zero = 1.;
//	double sum = z_zero;
//
//	while (true) {
//		++x;
//		int kout_g_p = gr_out + 1 - x;
//		double degree_node_p = k_degree + 1 - x;
//		double ga = ext_stubs - gr_out + x; // TODO: ext_stubs with or without k_degree?
//		z_zero *= kout_g_p / (x * ga) * degree_node_p;
//		if (z_zero < precision * sum)
//			break;
//		if (pzero * sum > 1)
//			return pzero;
//		sum += z_zero;
//	}
//
//	return pzero * sum;
//}

double LogFactTable::cum_hyper_left(int k_in, int gr_out, int open_stubs, int k_degree) {
	// this is strictly less  p(x < k_in)   *** NOT EQUAL ***
	//std::cout<<k_in<<" node: "<<degree_node<<" group: "<<open_stubs<<" "<<degree_node<<std::endl;
	if (k_in <= 0)
		return 0;

	if (open_stubs - gr_out - k_degree + k_in <= 0)
		return 0;

	if (k_in > std::min(k_degree, gr_out))
		return 1;

	if (k_in > double(gr_out + 1) / double(open_stubs + 2) * double(k_degree + 1))
		return (1. - cum_hyper_right(k_in, gr_out, open_stubs, k_degree));

	int x = k_in - 1;
	double pzero = hypergeom_dist(x, gr_out, open_stubs, k_degree);

	//std::cout<<"pzero: "<<pzero<<" "<<log_hyper(x, gr_out, open_stubs, degree_node)<<" gsl: "<<(gsl_ran_hypergeometric_pdf(x, gr_out, open_stubs - gr_out,  degree_node))<<std::endl;
	if (pzero <= 1e-40)
		return 0;

	double ga = open_stubs - gr_out - k_degree;
	int kout_g_p = gr_out + 1;
	double degree_node_p = k_degree + 1;
	double z_zero = 1.;
	double sum = z_zero;
	while (true) {
		z_zero *= (ga + x) / ((degree_node_p - x) * (kout_g_p - x)) * x;
		--x;
		if (z_zero < log_table_pr * sum)
			break;
		if (pzero * sum > 1)
			return pzero;
		sum += z_zero;
	}

	return pzero * sum;
}

double LogFactTable::cum_binomial_right(int x, int N, double prob) {
	// this is bigger  or equal p(x >= kin_node)   *** EQUAL ***
	//std::cout<<"x "<<x<<" N "<<N <<"  prob "<<prob<<std::endl;
	if (x <= 0)
		return 1;

	if (x > N)
		return 0;

	if (prob - 1 > -1e-11)
		return 1;

	if (x < N * prob)
		return 1 - cum_binomial_left(x, N, prob);

	double pzero = binom(x, N, prob);
	if (pzero <= 1e-40)
		return 0;

	double z_zero = 1.;
	double sum = z_zero;
	while (true) {
		z_zero *= prob * double(N - x) / ((x + 1) * (1 - prob));
		x++;
		//std::cout<<"zzero sum "<<z_zero<<" "<<sum<<" "<<std::endl;
		if (z_zero < log_table_pr * sum)
			break;
		sum += z_zero;
	}
	return pzero * sum;
}

double LogFactTable::cum_binomial_left(int x, int N, double prob) {
	// this is less strictly p(x < kin_node)   *** NOT EQUAL ***
	if (x <= 0)
		return 0;

	if (x > N)
		return 1;

	if (prob < 1e-11)
		return 1;

	if (x > N * prob)
		return 1 - cum_binomial_right(x, N, prob);

	--x;
	double pzero = binom(x, N, prob);
	if (pzero <= 1e-40)
		return 0;

	double z_zero = 1.;
	double sum = z_zero;
	while (true) {
		--x;
		z_zero *= (1 - prob) * double(x + 1) / ((N - x) * prob);
		//std::cout<<"zzero sum "<<z_zero<<" "<<sum<<" "<<(ga + x)<<std::endl;
		if (z_zero < log_table_pr * sum)
			break;
		sum += z_zero;
	}
	return pzero * sum;
}

double LogFactTable::slow_symmetric_eq(int k1, int k2, int H, int x) {
	// k1, k2 and k3 are the three colors
	//std::cout<<"k3: "<<k3<<std::endl;
	int l1 = std::max(0, -H);
	int l2 = std::min(k1, k2);
	//std::cout<<"l1: "<<l1<<" l2: "<<l2<<std::endl;
	if (x < l1)
		return 0;

	if (x > l2)
		return 0;

	double p = 0;
	for (int ix = l1; ix <= l2; ++ix) {
		//std::cout<<ix<<" "<<(log_symmetric_eq(k1, k2, H, ix))<<std::endl;
		p += std::exp(log_symmetric_eq(k1, k2, H, ix));
	}
	//std::cout<<"p: "<<p<<std::endl;
	return std::exp(log_symmetric_eq(k1, k2, H, x)) / p;
}

inline double
LogFactTable::fast_right_cum_symmetric_eq(int k1, int k2, int H, int x, int mode,
                                          int ext_stubs) {
	// I want k1 to be the smaller between the two
	//std::cout<<"k1 "<<k1<<" "<<k2<<" "<<H<<" "<<x<<" "<<mode<<" "<<2*H+k1+k2<<std::endl;
	if (k1 > k2)
		return fast_right_cum_symmetric_eq(k2, k1, H, x, mode, ext_stubs);

	double ri = 1;
	double q1 = 0;
	double q2 = 0;
	double ratio;

	if (x == mode)
		++q2;
	else
		++q1;

	int l1 = std::max(0, -H);

	double ii = mode - 1;

	while (ii >= l1) {
		ratio = sym_ratio(k1, k2, H, ii + 1);
		ri /= ratio;
		q1 += ri;
		if (q1 > 1e280)
			return cum_hyper_right(x, k2, ext_stubs, k1);
		if (ri < log_table_pr * q1)
			break;
		--ii;
	}

	/*double cum1=exp(log_symmetric_eq(k1, k2, H, l1));
	double lg0=log_symmetric_eq(k1, k2, H, l1);
	std::cout<<"x: "<<l1<<" "<<exp(log_symmetric_eq(k1, k2, H, l1))<<" "<<exp(lg0)<<std::endl;*/
	ri = 1;
	ii = mode + 1;
	//for(double i=mode+1; i<x; i++)
	while (ii < x) {
		ratio = sym_ratio(k1, k2, H, ii);
		ri *= ratio;
		q1 += ri;
		if (q1 > 1e280)
			return cum_hyper_right(x, k2, ext_stubs, k1);

		if (ri < log_table_pr * q1)
			break;
		//std::cout<<ii<<" "<<ratio<<" "<<ri/q1<<" b"<<std::endl;;
		++ii;
		//std::cout<<"dx-->: "<<ii<<" "<<exp(log_symmetric_eq(k1, k2, H, ii))<<" "<<exp(lg0) * ri<<" "<<sym_ratio(k1, k2, H, ii+1)<<std::endl;
	}
	ii = std::max(x, mode + 1);
	ri = exp(log_symmetric_eq(k1, k2, H, cast_int(ii - 1)) -
	         log_symmetric_eq(k1, k2, H, mode));

	//for(double i=max(x, mode+1); i<=k1; i++)
	while (ii <= k1) {
		ratio = sym_ratio(k1, k2, H, ii);
		ri *= ratio;
		q2 += ri;
		if (q2 > 1e280)
			return cum_hyper_right(x, k2, ext_stubs, k1);
		//std::cout<<ii<<" "<<ratio<<" "<<ri/q2<<" c "<<x<<" "<<q2<<std::endl;;
		++ii;
		if (ri < log_table_pr * q2)
			break;
		//std::cout<<"ddx-->: "<<i<<" "<<exp(log_symmetric_eq(k1, k2, H, i))<<" "<<exp(lg0) * ri<<" "<<sym_ratio(k1, k2, H, i+1)<<std::endl;
	}
	/* std::cout<<"fast q12: "<<q1<<" "<<q2<<std::endl;*/

	return std::max(q2 / (q1 + q2), 1e-100);
}

// k1 = node_degree
// k2 = g_out
// x = k_in
double LogFactTable::right_cumulative_function(int k1, int k2, int open_stubs, int x) {
	if (x > k1 || x > k2)
		return 0;

	if (k1 * k1 < open_stubs)
		return cum_hyper_right(x, k2, open_stubs, k1);

	//	std::cout << "Calculate propability directly" << std::endl;
	// k1 is the degree of the node
	// k2 is the degree of the other node (the bigger)
	// k3 is 2m - k1 - k2
	int k3 = open_stubs - k1;
	int H = (k3 - k1 - k2) / 2;
	int l1 = std::max(0, -H);

	if (x == l1)
		return 1;

	int mode = std::max(cast_int(k2 / double(k1 + k3) * k1),
	                    l1); // this mode is underestimated anyway
	if (mode > k2)
		mode = k2;

//	std::cout << "mode: " << mode << std::endl;
	if (x < mode)
		return cum_hyper_right(x, k2, open_stubs, k1);
	return fast_right_cum_symmetric_eq(k1, k2, H, x, mode, open_stubs);
}

double LogFactTable::log_factorial(int a) { return lnf[a]; }

void LogFactTable::delete_instance() {
	delete instance;
	instance = nullptr;
}
