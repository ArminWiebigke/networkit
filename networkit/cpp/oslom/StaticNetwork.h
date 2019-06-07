#pragma once

#ifndef STATIC_NETWORK_HPP
#define STATIC_NETWORK_HPP

#include <string>
#include <deque>
#include <set>
#include <map>
#include <algorithm>
#include <iostream>

#include "utils/Wsarray.h"
#include "ModuleCollection.h"
#include "Parameters.h"

class StaticNetwork {

public:

    StaticNetwork();

    ~StaticNetwork();

    int draw(std::string);

    void print_id(const std::deque<int> &a, std::ostream &);

    void print_id(const std::deque<std::deque<int>> &, std::ostream &);

    void print_id(const std::deque<std::set<int>> &, std::ostream &);

    void print_id(const std::set<int> &, std::ostream &);

    void deque_id(std::deque<int> &);

    void set_subgraph(std::deque<int> &group, std::deque<std::deque<int>> &link_per_node,
                      std::deque<std::deque<std::pair<int, double>>> &weights_per_node);

    int translate(std::deque<std::deque<int>> &);

    int translate_anyway(std::deque<std::deque<int>> &ten);

    void get_id_label(std::map<int, int> &);

    int id_of(int a) { return vertices[a]->id_num; };

    int size() { return dim; };

    double stubs() { return oneM; };

    int kin_m(const std::deque<int> &);

    int kin_m(const std::set<int> &);

    int ktot_m(const std::deque<int> &);

    int ktot_m(const std::set<int> &);

    void set_graph(std::map<int, std::map<int, std::pair<int, double> >> &A);

    bool set_graph(const std::string &file_name);

    void set_graph(std::deque<std::deque<int>> &link_per_node,
                   std::deque<std::deque<std::pair<int, double>>> &weights_per_node,
                   std::deque<int> &label_rows);

    void clear();

    void set_proper_weights();

    void set_connected_components(std::deque<std::deque<int>> &);

    int propagate_distances(std::deque<int> &new_shell, std::set<int> &already_gone,
                            std::deque<std::pair<int, int>> &distances_node, int shell,
                            std::deque<double> &ML, int &, int);

    void same_component(int, std::set<int> &);

    int set_upper_network(std::map<int, std::map<int, std::pair<int, double> >> &neigh_weight_f,
                          ModuleCollection &Mcoll);

    void print_degree_of_homeless(std::deque<int> &homel, std::ostream &outt);

    int draw_with_weight_probability(std::string file_name);

protected:

    class vertex {

    public:

        vertex(int, int, int);

        ~vertex();

        void kplus_global_and_quick(std::deque<int> &a, int &stubs_in, double &strength_in);

        int kplus_m(const std::deque<int> &);

        double kplus_w(const std::deque<int> &);

        int kplus_m(const std::set<int> &);

        int id_num;       // id
        double strength;  // sum of the weights
        int stub_number;  // number of stubs
        Wsarray *links;   // array with label of neighbor, multiple links, sm of the weights towards it
        // links->l[i] = label of neighbor i, links->w[i] = number of edges, sum of weights to the neighbor i
        std::deque<double> original_weights;
    };

    int dim;   // number of nodes
    int oneM;  // number of stubs

    std::deque<vertex *> vertices;

    Parameters *paras;
};

#endif




