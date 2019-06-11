#pragma once

#ifndef OSLOMNET_EVALUATE_HPP
#define OSLOMNET_EVALUATE_HPP

#include <deque>
#include <set>
#include <map>
#include <string>

#include "LouvainOslomnet.h"
#include "WeightedTabdeg.h"
#include "DataTypes.h"

// NOTE: CUP = Clean Up Procedure, cgroup = input cluster/community

inline double log_zero(double a) {

    if (a <= 0)
        return -1e20;
    else
        return std::log(a);
}

class oslom_net_evaluate : public oslomnet_louvain {

public:

    oslom_net_evaluate(std::deque<std::deque<int> > &b,
                       std::deque<std::deque<std::pair<int, double> > > &c,
                       std::deque<int> &d);;

    explicit oslom_net_evaluate(std::string a);;

    explicit oslom_net_evaluate(std::map<int, std::map<int, std::pair<int, double> > > &A);;

    ~oslom_net_evaluate() = default;

    double CUP_both(const std::deque<int> &group_to_clean, std::deque<int> &cleaned_group,
                    int number_of_runs = -1);

    double CUP_check(const std::deque<int> &_c_, std::deque<int> &gr_cleaned,
                     int number_of_runs = -1);

    double CUP_search(const std::deque<int> &_c_, std::deque<int> &gr_cleaned,
                      int number_of_runs = -1);

    double group_inflation(const std::deque<int> &_c_, std::deque<int> &gr_cleaned,
                           int number_of_runs = -1);

private:

    void erase_from_cgroup(int wnode);

    void insert_into_cgroup(int node);

    bool erase_the_worst(int &wnode);

    bool erase_the_worst();

    int set_maxbord();

    void set_cgroup_and_neighs(const std::deque<int> &G);

    double CUP_from_list(CupDataStruct &a, std::deque<int> &gr_cleaned);

    double CUP_from_list_simple(CupDataStruct &a, std::deque<int> &gr_cleaned);

    void
    get_external_scores(WeightedTabdeg &neighbors, CupDataStruct &fitness_label_to_sort,
                        int kout_g, int tm, int Nstar, int nneighs,
                        const double &max_r, bool only_c,
                        WeightedTabdeg &previous_tab_c);

    double CUP_runs(WeightedTabdeg &previous_tab_c, WeightedTabdeg &previous_tab_n,
                    int kin_cgroup_prev, int ktot_cgroup_prev, std::deque<int> &gr_cleaned,
                    bool only_c, int runs);

    void initialize_for_evaluation(const std::deque<int> &eval_group, WeightedTabdeg &previous_tab_c,
                                   WeightedTabdeg &previous_tab_n, int &kin_cgroup_prev,
                                   int &ktot_cgroup_prev);

    void initialize_for_evaluation(WeightedTabdeg &previous_tab_c,
                                   WeightedTabdeg &previous_tab_n, int &kin_cgroup_prev,
                                   int &ktot_cgroup_prev);

    double
    clean_up_procedure(WeightedTabdeg &previous_tab_c, WeightedTabdeg &previous_tab_n,
                       int kin_cgroup_prev, int ktot_cgroup_prev,
                       std::deque<int> &border_group,
                       bool only_c);

    void set_changendi_cum();

    void insertion(int changendi);

    bool insert_the_best();

    double CUP_iterative(const std::deque<int> &_c_, std::deque<int> &gr_cleaned,
                         int number_of_runs = -1);

    /* DATA ***************************************************/

    double max_r_bord;            // this is the maximum r allowed for the external nodes (we don't want to look at all the graph, it would take too long)
    int maxb_nodes;               // this is the maximum number of nodes allowed in the border (similar as above)
    std::deque<double> changendi_cum;  // this is the cumulative distribution of the number of nodes to add to the cluster in the group_inflation function

    // ************* things to update *************************
    WeightedTabdeg cgroup;  // Cluster
    WeightedTabdeg neighs;  // Neighbors
    int kin_cgroup;
    int ktot_cgroup;
    /*********************************************************/
};



#endif // OSLOMNET_EVALUATE_HPP