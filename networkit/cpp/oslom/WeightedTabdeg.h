#pragma once

#ifndef WEIGHTED_TABDEG_HPP
#define WEIGHTED_TABDEG_HPP

#include <map>
#include <deque>
#include <iostream>
#include <cmath>

# define num_up_to 5


/**
 * This class stores the properties (degree etc.) of one node.
 */
class facts {

public:
    facts(int a, double b, std::multimap<double, int>::iterator c, int d) {
        internal_degree = a;
        internal_edgeweight = b;
        fitness_iterator = c;
        degree = d;
    };

    ~facts() = default;;

    int degree; // Total degree of the node
    int internal_degree;  // Number of edges between the node and the group.

    // The sum of the edgeweights into the group. Identical to the internal degree if the graph
    // is unweighted.
    double internal_edgeweight; // wr is the right part of the exponential for the weights, this is the sum over the internal stubs of that

    std::multimap<double, int>::iterator fitness_iterator;
};

/**
 * This class stores information about a set of nodes, e.g. (1) the current group or (2) the
 * neighbors of the current group.
 */
class WeightedTabdeg {

public:

    WeightedTabdeg() = default;;

    ~WeightedTabdeg() = default;;

    void _set_(WeightedTabdeg &);

    void clear();

    void insert_node(int a, int kp, int kt, double mtlw, double fit);

    bool erase_node(int a);

    void set_deque(std::deque<int> &);

    int size();

    void print_nodes(std::ostream &, std::deque<int> &);

    void
    set_and_update_group(int nstar, int nn, int kout_g, int tm, WeightedTabdeg &one);

    void
    set_and_update_neighs(int nstar, int nn, int kout_g, int tm, WeightedTabdeg &one);

    bool update_group(int node, int delta_degree, double delta_mtlw, int nstar, int nn,
                      int kout_g, int tm, int kt, std::deque<int> &to_be_erased);

    bool update_neighs(int a, int delta_degree, double delta_mtlw, int nstar, int kout_g,
                       int tm, int kt);

    int
    best_node(int &lab, double &best_fitness, int kout_g, int Nstar, int nneighs, int tm);

    int worst_node(int &lab, double &worst_fitness, int kout_g, int Nstar, int nneighs,
                   int tm);

    bool is_internal(int a);

    std::map<int, facts> labels_to_facts;        // maps the label into the facts
    std::multimap<double, int> fitness_to_label; // maps the fitness into the label (this can be optimized)
};

#endif //WEIGHTED_TABDEG_HPP