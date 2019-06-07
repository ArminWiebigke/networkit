#pragma once

#ifndef WEIGHTED_TABDEG_HPP
#define WEIGHTED_TABDEG_HPP

#include <map>
#include <deque>
#include <iostream>
#include <cmath>

# define num_up_to 5


class facts {

public:
    facts(int a, double b, std::multimap<double, int>::iterator c, int d) {
        internal_degree = a;
        minus_log_total_wr = b;
        fitness_iterator = c;
        degree = d;
    };

    ~facts() = default;;

    int degree;
    int internal_degree;
    double minus_log_total_wr; // wr is the right part of the exponential for the weights, this is the sum over the internal stubs of that
    std::multimap<double, int>::iterator fitness_iterator;
};

class WeightedTabdeg {

public:

    WeightedTabdeg() = default;;

    ~WeightedTabdeg() = default;;

    void _set_(WeightedTabdeg &);

    void clear();

    void edinsert(int a, int kp, int kt, double mtlw, double fit);

    bool erase(int a);

    void set_deque(std::deque<int> &);


    int size() { return lab_facts.size(); };

    void print_nodes(std::ostream &, std::deque<int> &);

    void
    set_and_update_group(int nstar, int nn, int kout_g, int tm, WeightedTabdeg &one);

    void
    set_and_update_neighs(int nstar, int nn, int kout_g, int tm, WeightedTabdeg &one);

    bool update_group(int a, int delta_degree, double delta_mtlw, int nstar, int nn,
                      int kout_g, int tm, int kt, std::deque<int> &to_be_erased);

    bool update_neighs(int a, int delta_degree, double delta_mtlw, int nstar, int kout_g,
                       int tm, int kt);

    int
    best_node(int &lab, double &best_fitness, int kout_g, int Nstar, int nneighs, int tm);

    int worst_node(int &lab, double &worst_fitness, int kout_g, int Nstar, int nneighs,
                   int tm);

    bool is_internal(int a);

    std::map<int, facts> lab_facts;         // maps the label into the facts
    std::multimap<double, int> fitness_lab; // maps the fitness into the label  (this can be optimized)
};

#endif //WEIGHTED_TABDEG_HPP