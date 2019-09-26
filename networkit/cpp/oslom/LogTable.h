#pragma once

#ifndef LOG_TABLE_HPP
#define LOG_TABLE_HPP

#include <vector>
#include <cmath>
#include <memory>

#include "utils/Cast.h"

#define log_table_pr 1e-5

class LogFactTable {

private:

    LogFactTable() = default;

    static LogFactTable *instance;

public:

    static LogFactTable *get_instance();

	static void delete_instance();

    ~LogFactTable() = default;

    double log_factorial(int a);

    /**
     * Memoize the log sums from 0 to size.
     * log sum (x) = sum(log(i)) for i = 0 to x
     */
    void set(int size);

    inline double log_hyper(int k_in, int g_out, int open_stubs, int k_degree) {
        return log_choose(g_out, k_in)
               + log_choose(open_stubs - g_out, k_degree - k_in)
               - log_choose(open_stubs, k_degree);
    };

	/**
	 * Calculate the probability that the node has exactly k_in edges into the group, according to
	 * a hypergeometric distribution.
	 * P(X = kin_node)
	 * @param k_in Edges from node to group
	 * @param g_out Outgoing stubs of group
	 * @param open_stubs Open stubs, g_out + M (sum of degree of all external nodes)
	 *                           or: 2 * m - g_in - k
	 * @param k_degree Degree of the node
	 * @return the probability
	 */
    inline double hypergeom_dist(int k_in, int g_out, int open_stubs, int k_degree) {
        return std::max(0., std::exp(log_hyper(k_in, g_out, open_stubs, k_degree)));
    };

    inline double binom(int x, int N, double p) {
        return std::exp(log_choose(N, x) + x * std::log(p) + (N - x) * std::log(1 - p));
    };

    /**
     * Calculate the number of possible configurations to choose k out of n
     * == binomial coefficient n choose k
     * @param n
     * @param k
     * @return natural logarithm of the binomiala coefficent
     */
    inline double log_choose(int n, int k) {
        return lnf[n] - lnf[n - k] - lnf[k];
    };

    double cum_binomial_right(int x, int N, double prob);

    double cum_binomial_left(int x, int N, double prob);

    inline double log_symmetric_eq(int k1, int k2, int H, int x) {
        return -x * lnf[2] - lnf[k1 - x] - lnf[k2 - x] - lnf[x + H] - lnf[x];
    };

    double slow_symmetric_eq(int k1, int k2, int H, int x);

    double fast_right_cum_symmetric_eq(int kTotal, int cOut, int MInEdges, int kIn, int mode, int extStubs);

    double right_cumulative_function(int kTotal, int cOut, int extStubs, int kIn);

private:

    std::vector<double> lnf;

    inline double sym_ratio(int &kTotal, int &cOut, int &MInEdges, double kIn) {
        return 0.5 * (kTotal - kIn + 1) / ((kIn + MInEdges) * kIn) * (cOut - kIn + 1);
    };

    double cum_hyper_right(int k_in, int gr_out, int open_stubs, int k_degree);

    double cum_hyper_left(int k_in, int gr_out, int open_stubs, int k_degree);
};

#endif // LOG_TABLE_HPP