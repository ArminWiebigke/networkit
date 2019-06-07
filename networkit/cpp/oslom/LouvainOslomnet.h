#pragma once

#ifndef LOUVAIN_OSLOMNET_HPP
#define LOUVAIN_OSLOMNET_HPP

#include <map>
#include <iostream>

#include "StaticNetwork.h"
#include "Stochastics.h"

struct oslom_module {

    explicit oslom_module(int a) {
        nc = 1;
        kout = a;
        ktot = a;
    };

    ~oslom_module() = default;;

    int nc;
    int kout;
    int ktot;
};

typedef std::map<int, std::pair<int, double> > mapip;
typedef std::map<int, oslom_module> map_int_om;

void prints(map_int_om &M);

class oslomnet_louvain : public StaticNetwork {

public:

    oslomnet_louvain() : StaticNetwork() {};

    ~oslomnet_louvain() = default;;

    int collect_raw_groups_once(std::deque <std::deque<int>> &);

private:

    void weighted_favorite_of(const int &node, int &fi, int &kp, int &kop);

    void unweighted_favorite_of(const int &node, int &new_label, int &stubs_into_new, int &stubs_into_old);

    void single_pass_unweighted();

    void single_pass_weighted();

    inline void update_modules(const int &node, const int &new_label, const int &stubs_into_new, const int &stubs_into_old);

    void module_initializing();

    void set_partition_collected(std::deque <std::deque<int>> &M);

    //int check_all();
    std::map<int, oslom_module> label_module;
    std::deque<int> vertex_label;
    std::deque<int> vertex_order;
    std::deque<bool> vertex_to_check;
    std::deque<bool> vertex_to_check_next;
    int nodes_changed;
};

#endif // LOUVAIN_OSLOMNET_HPP