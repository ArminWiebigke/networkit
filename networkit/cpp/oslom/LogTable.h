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

    void set(int);
    //void set_small_tab_right_hyper(int, int);

    inline double log_hyper(int kin_node, int kout_g, int tm, int degree_node) {
        return log_choose(kout_g, kin_node) +
               log_choose(tm - kout_g, degree_node - kin_node) -
               log_choose(tm, degree_node);
    };

    inline double hyper(int kin_node, int kout_g, int tm, int degree_node) {
        return std::max(0., std::exp(log_hyper(kin_node, kout_g, tm, degree_node)));
    };

    inline double binom(int x, int N, double p) {
        return std::exp(log_choose(N, x) + x * std::log(p) + (N - x) * std::log(1 - p));
    };

    inline double log_choose(int tm, int degree_node) {
        return lnf[tm] - lnf[tm - degree_node] - lnf[degree_node];
    };

    double cum_binomial_right(int x, int N, double prob);

    double cum_binomial_left(int x, int N, double prob);

    inline double log_symmetric_eq(int k1, int k2, int H, int x) {
        return -x * lnf[2] - lnf[k1 - x] - lnf[k2 - x] - lnf[x + H] - lnf[x];
    };

    double slow_symmetric_eq(int k1, int k2, int H, int x);

    //vector<vector<vector<vector<double> > > > small_rh;
    // /*********/// small_rh[tm][kout][k][kin]		tm>=kout>=k>=kin
    //double right_cum_symmetric_eq(int k1, int k2, int H, int x);
    double fast_right_cum_symmetric_eq(int k1, int k2, int H, int x, int mode, int tm);

    double right_cumulative_function(int k1, int k2, int k3, int x);

private:

    std::vector<double> lnf;

    inline double sym_ratio(int &k1, int &k2, int &H, double i) {
        return 0.5 * (k1 - i + 1) / ((i + H) * i) * (k2 - i + 1);
    };

    double cum_hyper_right(int kin_node, int kout_g, int tm, int degree_node);

    double cum_hyper_left(int kin_node, int kout_g, int tm, int degree_node);
};

#endif // LOG_TABLE_HPP