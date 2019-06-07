#pragma once

#ifndef OSLOM_NET_GLOBAL_HPP
#define OSLOM_NET_GLOBAL_HPP

#include "OslomnetEvaluate.h"
#include "utils/Types.h"

class OslomNetGlobal : public oslom_net_evaluate {

public:

    OslomNetGlobal(IntMatrix &b,
                     std::deque<std::deque<std::pair<int, double> > > &c, std::deque<int> &d);

    explicit OslomNetGlobal(std::string a);

    explicit OslomNetGlobal(std::map<int, std::map<int, std::pair<int, double> > > &A);

    ~OslomNetGlobal() = default;

    ModuleCollection clean_up(const std::vector<std::deque<int>> &modules, int upper_node_id);
    void print_modules(bool not_homeless, const std::string &tp, ModuleCollection &Mcoll);

    void print_modules(bool not_homeless, std::ostream &out1, ModuleCollection &Mcoll);

    int try_to_assign_homeless(ModuleCollection &Mcoll, bool anyway);

    void print_statistics(std::ostream &outt, ModuleCollection &Mcoll);

private:

	void try_add_good_group(std::deque<int> &group, double &b_score,
	                        const std::deque<int> &original_group,
	                        IntMatrix &good_modules, std::deque<double> &bscores_good,
	                        IntMatrix &bad_groups);

    int try_to_merge_discarded(IntMatrix &discarded,
                               IntMatrix &good_modules,
                               std::deque<double> &bscores_good,
                               IntMatrix &new_discarded);

    void get_single_trial_partition(IntMatrix &good_modules_to_prune,
                                    std::deque<double> &bscores_good);

    void
    single_gather(IntMatrix &good_modules_to_prune, std::deque<double> &bscores_good, int);

    void check_minimality_all(IntMatrix &A, std::deque<double> &bss, ModuleCollection &minimal_modules);

    void
    check_minimality_matrix(IntMatrix &A, std::deque<double> &bss, ModuleCollection &minimal_modules,
                            IntMatrix &suggestion_matrix, std::deque<double> &suggestion_bs,
                            int counter);

    bool check_minimality(std::deque<int> &group, double &bs_group,
                          ModuleCollection &minimal_modules,
                          IntMatrix &suggestion_matrix, std::deque<double> &suggestion_bs);

    int check_unions_and_overlap(ModuleCollection &mall, bool only_similar = false);

    void check_existing_unions(ModuleCollection &mall);

    bool fusion_module_its_subs(const std::deque<int> &A, IntMatrix &its_submodules);

    bool fusion_with_empty_A(IntMatrix &its_submodules, std::deque<int> &grc1, double &bs);

    bool check_fusion_with_gather(ModuleCollection &mall);

    int check_intersection(ModuleCollection &Mcoll);

    int check_intersection(std::deque<int> &to_check, ModuleCollection &Moll);

    int
    fusion_intersection(std::set<std::pair<int, int> > &pairs_to_check, ModuleCollection &Mcoll);

    bool decision_fusion_intersection(int ai1, int ai2, std::deque<int> &new_insertions,
                                      ModuleCollection &Mcoll,
                                      double prev_over_percentage);
};

#endif // OSLOM_NET_GLOBAL_HPP