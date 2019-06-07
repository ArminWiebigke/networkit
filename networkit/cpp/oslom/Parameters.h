#pragma once

#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <string>
#include <map>
#include <deque>
#include <iostream>
#include <algorithm>
#include <vector>

#include "utils/Types.h"

struct Parameters {

private:

    Parameters();

    static Parameters* instance;

public:

    ~Parameters() = default;

    static Parameters* get_instance();

	static void delete_instance();

    void set(const std::vector<std::string> &args);

    void print();

    double threshold;

    int clean_up_runs;
    int inflate_runs;
    int inflate_stopper;
    double equivalence_parameter;
    int cut_off;

    int maxborder_nodes;
    double maxbg_ordinary;
    int iterative_stopper;
    int minimality_stopper;

    double coverage_inclusion_module_collection;
    double coverage_percentage_fusion_left;
    double check_inter_p;
    double coverage_percentage_fusion_or_submodules;

    bool print_flag_subgraph;

    bool fast;
    bool weighted;
    bool assign_homeless;

    int max_iteration_convergence;

    bool keep_bad_groups;
    bool check_unions;
    bool check_minimality;
    bool simple_cleanup;
    std::string cleanup_strategy;
    double max_group_extend;
    std::string bad_groups_filename;
    bool merge_discarded;
    bool discard_max_extend_groups;
};

#endif // PARAMETERS_HPP