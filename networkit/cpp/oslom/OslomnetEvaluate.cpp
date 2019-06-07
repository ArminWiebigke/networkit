#include "OslomnetEvaluate.h"
#include "Stochastics.h"

void oslom_net_evaluate::set_changendi_cum() {
    if (dim != 0 && oneM != 0) {
        int flat_until = cast_int(oneM / dim * 3);
        flat_until = std::min(dim / 2, flat_until);

        int max_p = std::max(paras->cut_off,
                             flat_until); // this is something which might be optimized
        max_p = std::min(dim / 2, max_p);

        powerlaw(max_p, flat_until + 1, 3, changendi_cum);
        std::deque<double> distr;
        distribution_from_cumulative(changendi_cum, distr);
        double ac = 1;

        if (!distr.empty())
            ac = distr[0];

        for (int i = 0; i < flat_until; i++)
            distr.push_front(ac);

        deque_numeric::normalize_one(distr);
        cumulative_from_distribution(changendi_cum, distr);
    }
}

int oslom_net_evaluate::set_maxbord() {
    max_r_bord = paras->maxbg_ordinary;
    maxb_nodes = paras->maxborder_nodes;
    return 0;
}

/**
 * Move a node out of a module and add it to the neighbors. Update all relevant things.
 * @param wnode The node to be moved
 */
void oslom_net_evaluate::erase_cgroup(int wnode) {
    auto itm = cgroup.lab_facts.find(wnode);
    if (itm != cgroup.lab_facts.end()) {
        int kp = itm->second.internal_degree;
        int kt = itm->second.degree;
        double mtlw = itm->second.minus_log_total_wr;

        kin_cgroup -= 2 * kp;
        ktot_cgroup -= kt;
        int kout_g = ktot_cgroup - kin_cgroup;
        int tm = oneM - ktot_cgroup;

        double fi = Stochastics::compute_global_fitness_ofive(kp, kout_g, tm, kt, mtlw,
                                                              neighs.size() + 1,
                                                              dim - cgroup.size() + 1);
        neighs.edinsert(wnode, kp, kt, mtlw, fi);
        cgroup.erase(wnode);

        std::deque<int> tobe;
        for (int i = 0; i < vertices[wnode]->links->size(); i++) {
            if (!cgroup.update_group(vertices[wnode]->links->l[i],
                                     -vertices[wnode]->links->w[i].first,
                                     -vertices[wnode]->links->w[i].second,
                                     dim - cgroup.size(), neighs.size(), kout_g, tm,
                                     vertices[vertices[wnode]->links->l[i]]->stub_number,
                                     tobe))
                neighs.update_neighs(vertices[wnode]->links->l[i],
                                     -vertices[wnode]->links->w[i].first,
                                     -vertices[wnode]->links->w[i].second,
                                     dim - cgroup.size(), kout_g, tm,
                                     vertices[vertices[wnode]->links->l[i]]->stub_number);
        }
        for (int i : tobe)
            erase_cgroup(i);
    }
}

/**
 * Find the worst node in the module and remove it from the module.
 * @param wnode (out): The node that was removed
 * @return True iff a node was removed
 */
bool oslom_net_evaluate::erase_the_worst(int &wnode) {
    // this function is to look for the worst node in cgroup and to erase it
    int Nstar = dim - cgroup.size();
    int nn = neighs.size();
    int kout_g = ktot_cgroup - kin_cgroup;
    int tm = oneM - ktot_cgroup;

    double wf;
    cgroup.worst_node(wnode, wf, kout_g, Nstar, nn, tm);

    if (cgroup.size() == 0) {
        return false;
    }
    erase_cgroup(wnode);
    return true;
}

/**
 * Insert a node into the module. Update all relevant things.
 * @param wnode The node that should be inserted
 */
void oslom_net_evaluate::insert_cgroup(int wnode) {
    // it needs to be differenciated between weighted and unweighted
    int kp, kt;
    double mtlw;
    {
        auto itm = neighs.lab_facts.find(wnode);
        if (itm != neighs.lab_facts.end()) {
            kp = itm->second.internal_degree;
            kt = itm->second.degree;
            mtlw = itm->second.minus_log_total_wr;
        } else {

            kp = 0;
            kt = vertices[wnode]->stub_number;
            mtlw = 0;
        }
    }

    int kout_g = ktot_cgroup - kin_cgroup;
    int tm = oneM - ktot_cgroup;
    double fi = Stochastics::compute_global_fitness_ofive(kp, kout_g, tm, kt, mtlw, neighs.size(),
                                                          dim - cgroup.size());
    kin_cgroup += 2 * kp;
    ktot_cgroup += kt;
    kout_g = ktot_cgroup - kin_cgroup;
    tm = oneM - ktot_cgroup;

    cgroup.edinsert(wnode, kp, kt, mtlw, fi);
    neighs.erase(wnode);

    std::deque<int> tobe;
    for (int i = 0; i < vertices[wnode]->links->size(); i++) {
        if (!cgroup.update_group(vertices[wnode]->links->l[i],
                                 vertices[wnode]->links->w[i].first,
                                 vertices[wnode]->links->w[i].second,
                                 dim - cgroup.size(), neighs.size(), kout_g, tm,
                                 vertices[vertices[wnode]->links->l[i]]->stub_number,
                                 tobe))
            neighs.update_neighs(vertices[wnode]->links->l[i],
                                 vertices[wnode]->links->w[i].first,
                                 vertices[wnode]->links->w[i].second,
                                 dim - cgroup.size(), kout_g, tm,
                                 vertices[vertices[wnode]->links->l[i]]->stub_number);
    }
}

/**
 * Try to throw nodes out of the module.
 * @param _c_ The module to check
 * @param gr_cleaned The cleaned module
 * @param number_of_runs Number of iterations of the cleanup
 * @return
 */
double oslom_net_evaluate::CUP_check(const std::deque<int> &_c_, std::deque<int> &gr_cleaned,
                                     int number_of_runs) {
    if (number_of_runs == -1)
        number_of_runs = paras->clean_up_runs;
//    std::cout << "CUP_check" << std::endl;
    /*_c_ is the module to clean up and gr_cleaned is the result */
    WeightedTabdeg previous_tab_c;
    WeightedTabdeg previous_tab_n;
    int kin_cgroup_prev;
    int ktot_cgroup_prev;
    double bscore;

    initialize_for_evaluation(_c_, previous_tab_c, previous_tab_n, kin_cgroup_prev,
                              ktot_cgroup_prev);
    bscore = CUP_runs(previous_tab_c, previous_tab_n, kin_cgroup_prev, ktot_cgroup_prev,
                      gr_cleaned, true, number_of_runs);

    return bscore;
}

/**
 * Try to add neighbors to the module.
 * @param _c_
 * @param gr_cleaned
 * @param number_of_runs
 * @return
 */
double oslom_net_evaluate::CUP_search(const std::deque<int> &_c_, std::deque<int> &gr_cleaned,
                                      int number_of_runs) {
    if (number_of_runs == -1)
        number_of_runs = paras->clean_up_runs;
    /*_c_ is the module to clean up and gr_cleaned is the result */
//    std::cout << "CUP_search" << std::endl;
    WeightedTabdeg previous_tab_c;
    WeightedTabdeg previous_tab_n;
    int kin_cgroup_prev;
    int ktot_cgroup_prev;
    double bscore;

    initialize_for_evaluation(_c_, previous_tab_c, previous_tab_n, kin_cgroup_prev,
                              ktot_cgroup_prev);
    bscore = CUP_runs(previous_tab_c, previous_tab_n, kin_cgroup_prev, ktot_cgroup_prev,
                      gr_cleaned, false, number_of_runs);

    return bscore;
}

double oslom_net_evaluate::CUP_both(const std::deque<int> &_c_, std::deque<int> &gr_cleaned,
                                    int number_of_runs) {
    if (number_of_runs == -1)
        number_of_runs = paras->clean_up_runs;
//    std::cout << "CUP_both" << std::endl;
    /*_c_ is the module to clean up and gr_cleaned is the result */
    WeightedTabdeg previous_tab_c;
    WeightedTabdeg previous_tab_n;
    int kin_cgroup_prev;
    int ktot_cgroup_prev;
    double bscore;

    initialize_for_evaluation(_c_, previous_tab_c, previous_tab_n, kin_cgroup_prev,
                              ktot_cgroup_prev);
    bscore = CUP_runs(previous_tab_c, previous_tab_n, kin_cgroup_prev, ktot_cgroup_prev,
                      gr_cleaned, false, number_of_runs);

    initialize_for_evaluation(gr_cleaned, previous_tab_c, previous_tab_n, kin_cgroup_prev,
                              ktot_cgroup_prev);
    bscore = CUP_runs(previous_tab_c, previous_tab_n, kin_cgroup_prev, ktot_cgroup_prev, gr_cleaned,
                      true, /*this "true" means I can only look at nodes in previous_tab_c */
                      number_of_runs);
    return bscore;
}

/**
 * Run multiple cleanups on a single module. If enough runs return that the module is significant,
 * the module is considered significant.
 * @param previous_tab_c
 * @param previous_tab_n
 * @param kin_cgroup_prev
 * @param ktot_cgroup_prev
 * @param gr_cleaned
 * @param only_c
 * @param number_of_runs
 * @return B-Score of the module
 */
double oslom_net_evaluate::CUP_runs(WeightedTabdeg &previous_tab_c,
                                    WeightedTabdeg &previous_tab_n,
                                    int kin_cgroup_prev,
                                    int ktot_cgroup_prev,
                                    std::deque<int> &gr_cleaned,
                                    bool only_c,
                                    int number_of_runs) {
    // These if statements are here to speed up the program if there are big clusters
    if (previous_tab_c.size() > 100000)
        number_of_runs = std::min(3, number_of_runs);
    else if (previous_tab_c.size() > 10000)
        number_of_runs = std::min(5, number_of_runs);
    else if (previous_tab_c.size() > 1000)
        number_of_runs = std::min(10, number_of_runs);

    gr_cleaned.clear();

    if (previous_tab_c.size() == 0)
        return 1;

    int max_gr_size = 0;
    double bscore = 1;
    int good_runs = 0;

    for (int i = 0; i < number_of_runs; i++) {
        std::deque<int> gr_run_i;
        double score_i = clean_up_procedure(previous_tab_c, previous_tab_n, kin_cgroup_prev,
                                            ktot_cgroup_prev, gr_run_i, only_c);

        if (cgroup.size() + int(gr_run_i.size()) > max_gr_size) {
            bscore = score_i;
            cgroup.set_deque(gr_cleaned);
            for (int j : gr_run_i)
                gr_cleaned.push_back(j);

            max_gr_size = gr_cleaned.size();
            sort(gr_cleaned.begin(), gr_cleaned.end());
        }

        if (!gr_run_i.empty()) {
            ++good_runs;
            if (good_runs >= 0.55 * number_of_runs)
                return bscore;
        }
    }

    if (good_runs < 0.55 * number_of_runs) {
        gr_cleaned.clear();
        bscore += paras->threshold;
        bscore = std::min(1., bscore);
    }
    return bscore;
}

/**
 * Cleanup a single module (one iteration).
 * @param previous_tab_c
 * @param previous_tab_n
 * @param kin_cgroup_prev
 * @param ktot_cgroup_prev
 * @param border_group (out): The neighbors that were added to the module.
 * @param only_c (in): If true, only consider nodes inside the module, so no neighbors are added
 * @return The B-Score of the module
 */
double oslom_net_evaluate::clean_up_procedure(
        WeightedTabdeg &previous_tab_c,  // Cluster
        WeightedTabdeg &previous_tab_n,  // Neighbors
        int kin_cgroup_prev,
        int ktot_cgroup_prev,
        std::deque<int> &border_group,
        bool only_c) {
    // still there is some stochasticity due to possible ties
    /*	previous_stuff is the module-stuff before the CUP (Clean Up Procedure)
        cgroup + border_group is the module cleaned	*/
    border_group.clear();
    cgroup._set_(previous_tab_c);
    neighs._set_(previous_tab_n);
    kin_cgroup = kin_cgroup_prev;
    ktot_cgroup = ktot_cgroup_prev;

    if (cgroup.size() == dim) {
        return 1;
    }

    double bscore = 1;

    while (true) {
        bscore = all_external_test(ktot_cgroup - kin_cgroup,
                                   oneM - ktot_cgroup,
                                   dim - cgroup.size(),
                                   neighs.size(),
                                   maxb_nodes / double(dim - cgroup.size()),
                                   max_r_bord,
                                   border_group, only_c, previous_tab_c);
        // If the group is not empty, the module is considered significant
        if (!border_group.empty())
            break;
        // If all nodes have been removed, the module is insignificant
        if (cgroup.size() == 0)
            break;

        int wnode;
        erase_the_worst(wnode);
    }
//    if (only_c && !border_group.empty())
//        throw std::runtime_error("Added neighbors!");
    return bscore;
}

double oslom_net_evaluate::all_external_test(int kout_g,
                                             int tm,
                                             int Nstar,
                                             int nneighs,
                                             const double &max_r_one,
                                             const double &maxr_two,
                                             std::deque<int> &gr_cleaned,
                                             bool only_c,
                                             WeightedTabdeg &previous_tab_c) {
    double max_r = std::min(max_r_one, maxr_two);
    CupDataStruct fitness_label_to_sort;
    get_external_scores(neighs, fitness_label_to_sort, kout_g, tm, Nstar, nneighs, max_r,
                        only_c, previous_tab_c);
    double bscore;
    if (paras->simple_cleanup)
        bscore = CUP_from_list_simple(fitness_label_to_sort, gr_cleaned);
    else
        bscore = CUP_from_list(fitness_label_to_sort, gr_cleaned);
    return bscore;
}

/**
 * Get the scores of all nodes of the module and all neighbors.
 * @param neighbors
 * @param fitness_label_to_sort
 * @param kout_g
 * @param tm
 * @param Nstar
 * @param nneighs
 * @param max_r
 * @param only_c If true, only get the scores of the original group, ignoring neighbors.
 * @param previous_tab_c
 */
void oslom_net_evaluate::get_external_scores(WeightedTabdeg &neighbors,
                                             CupDataStruct &fitness_label_to_sort,
                                             int kout_g, int tm, int Nstar, int nneighs,
                                             const double &max_r, bool only_c,
                                             WeightedTabdeg &previous_tab_c) {
    auto bit = neighbors.fitness_lab.begin();
    int counter = 0;

    while (bit != neighbors.fitness_lab.end()) {
        auto itm = neighbors.lab_facts.find(bit->second);
        double interval;
        double F = Stochastics::compute_global_fitness(itm->second.internal_degree, kout_g, tm,
                                                       itm->second.degree,
                                                       itm->second.minus_log_total_wr, nneighs,
                                                       Nstar,
                                                       interval);
        if (F > max_r) {
            // r-Score is too bad, so we ignore this node.
            /*if(only_c == false || previous_tab_c.is_internal(itm->first))
                cout<<"no node: "<<vertices[itm->first]->id_num<<"  "<<itm->second.internal_degree<<" / "<< itm->second.degree<<" fitness: "<<F<<endl;*/
            counter++;
            if (counter > num_up_to)
                break;
        } else {
            /*if(only_c == false || previous_tab_c.is_internal(itm->first))
                cout<<"node: "<<"  "<<itm->second.internal_degree<<" / "<< itm->second.degree<<" fitness: "<<F<<" "<<interval<<endl;*/
            if (previous_tab_c.is_internal(itm->first) || !only_c)
                fitness_label_to_sort.insert(make_pair(F, std::make_pair(itm->first, interval)));
        }
        bit++;
    }
}

/**
 * Run the cleanup procedure from a list of nodes/scores.
 * @param a (in): Contains the scores of the nodes
 * @param gr_cleaned  (out): The nodes in the cleaned module
 * @return
 */
double oslom_net_evaluate::CUP_from_list(CupDataStruct &a, std::deque<int> &gr_cleaned) {
    int Nstar; // Number of neighbors or number of external nodes
    if (paras->weighted)
        Nstar = neighs.size();
    else
        Nstar = dim - cgroup.size();

    double critical_xi = -log(1 - paras->threshold) / Stochastics::fitted_exponent(Nstar);
    int pos = Nstar;
    int add_nodes = 0;                      // add_nodes tells how many nodes should be included into the cluster
    double probability_a = 0, probability_b = 0; // these are the two extremes of a possible good node I could have found
    double c_min = 1;                    // this is the score we give to the border we are evaluating here
    //cout<<"critical_xi: "<<critical_xi<<" --------------------------------------- "<<neighs.size()<<" cgroup "<<cgroup.size()<<endl<<endl<<endl;

    auto itl = a.begin();

    while (itl != a.end()) {
        double score = itl->first;
        double score_interval = itl->second.second;
        /*
         c_pos is the probability that Omega_{pos} <= score == Phi(c_pos)
         If this is low, the score of the node is better than expected in the null model,
         so we can assume that the node is part of the community
        */
        double c_pos = Stochastics::order_statistics_left_cumulative(Nstar, pos, score);
        //cout<<"position .... "<<pos<<" "<<order_statistics_left_cumulative(Nstar, pos, itl->first + itl->second.second)<<" >>AAA<< "<<Nstar<<" +++ "<<itl->first + itl->second.second<<" "<<itl->first - itl->second.second<<endl;
        c_min = std::min(c_pos, c_min);

        if (c_pos < critical_xi) {
            /*
            this is the basic condition of the order statistics test
            it's saying: look, this guy (itl->second.first) has an average fitness (itl->first)
            whose order_statistics_left_cumulative is below the threshold
            */
            if (add_nodes == 0) {
                // this node is the first node to be below the threshold
                add_nodes = Nstar - pos + 1;
                c_min = c_pos;
                probability_a = score - score_interval;
                probability_b = score + score_interval;
            } else {
                /*
                the previous node was already below the threshold.
                In this case I need to know if I should stop now or go on.
                The condition is related to the probability_to_overtake the previous guy
                */
                double probability_to_overtake = Stochastics::compare_r_variables(
                        probability_a, probability_b, score - score_interval,
                        score + score_interval);
                if (probability_to_overtake > 0.4999) {
                    /*preliminary check: this node is basically equivalent to the previous guy, I consider it good*/
                    add_nodes = Nstar - pos + 1;
                    c_min = c_pos;
                    probability_a = score - score_interval;
                    probability_b = score + score_interval;
                } else {
                    /*now I need to compute the bootstrap probability that the previous guy would have stopped the process*/
                    auto calc_probability = [&]() {
                        auto probability_to_stop = Stochastics::compute_probability_to_stop(
                                probability_a, probability_b, critical_xi, Nstar, pos + 1);
                        return (1. - probability_to_overtake) * probability_to_stop;
                    };
                    if (probability_to_overtake == 0 || calc_probability() > 0.5001) {
                        if (Stochastics::equivalent_check_gather(a, add_nodes, probability_a,
                                                                 probability_b, Nstar, critical_xi))
                            break;
                    }
                    add_nodes = Nstar - pos + 1;
                    c_min = c_pos;
                    probability_a = score - score_interval;
                    probability_b = score + score_interval;
                }
            }
        } else {  /* this node is not below the threshold */
            if (add_nodes != 0) {
                /* this means that this node is not good and the previous one was good. So, I stop here */
                if (Stochastics::equivalent_check_gather(a, add_nodes, probability_a, probability_b,
                                                         Nstar, critical_xi))
                    break;
            }
        }

        --pos;
        ++itl;
    }
    // equalizer check
    // this check is important to see if the procedure stopped just because there were a lot of equivalents nodes
    if (add_nodes != 0 && itl == a.end())
        Stochastics::equivalent_check_gather(a, add_nodes, probability_a, probability_b, Nstar,
                                             critical_xi);

    // Insert nodes in the cleaned group
    itl = a.begin();
    for (int i = 0; i < add_nodes; ++i) {
        gr_cleaned.push_back(itl->second.first);
        ++itl;
    }
    return Stochastics::pron_min_exp(Nstar, c_min);
}

/**
 * Run the cleanup procedure from a list of nodes/scores.
 * @param a (in): Contains the scores of the nodes
 * @param gr_cleaned  (out): The nodes in the cleaned module
 * @return
 */
double oslom_net_evaluate::CUP_from_list_simple(CupDataStruct &a, std::deque<int> &gr_cleaned) {
    int Nstar; // Number of neighbors or number of external nodes
    if (paras->weighted)
        Nstar = neighs.size();
    else
        Nstar = dim - cgroup.size();

    double critical_xi = -log(1 - paras->threshold) / Stochastics::fitted_exponent(Nstar);
    int pos = Nstar;
    int add_nodes = 0; // add_nodes tells how many nodes should be included into the cluster
    double c_min = 1;  // this is the score we give to the border we are evaluating here

    auto itl = a.begin();
//    std::cout << Nstar << "(" << critical_xi << ") ----------------" << std::endl;
    while (itl != a.end()) {
        double score = itl->first;
        /*
         c_pos is the probability that Omega_{pos} <= score == Phi(c_pos)
         If this is low, the score of the node is better than expected in the null model,
         so we can assume that the node is part of the community
        */
        double c_pos = Stochastics::order_statistics_left_cumulative(Nstar, pos, score);
        c_min = std::min(c_pos, c_min);
//        if (c_pos < 0.1)
//            std::cout << score << " - " << c_pos << std::endl;

        if (c_pos < critical_xi) {
            /*
            this is the basic condition of the order statistics test
            it's saying: look, this guy (itl->second.first) has an average fitness (itl->first)
            whose order_statistics_left_cumulative is below the threshold
            */
            add_nodes = Nstar - pos + 1;
            c_min = c_pos;
        } else {  /* this node is not below the threshold */
            if (add_nodes != 0) {
                /* this means that this node is not good and the previous one was good. So, I stop here */
                break;
            }
        }

        --pos;
        ++itl;
    }

    // Insert nodes in the cleaned group
    itl = a.begin();
    for (int i = 0; i < add_nodes; ++i) {
        gr_cleaned.push_back(itl->second.first);
        ++itl;
    }
    return Stochastics::pron_min_exp(Nstar, c_min);
}

double oslom_net_evaluate::CUP_iterative(const std::deque<int> &_c_, std::deque<int> &gr_cleaned,
                                         int number_of_runs) {
    if (number_of_runs == -1)
        number_of_runs = paras->clean_up_runs;
    double bs = CUP_both(_c_, gr_cleaned, number_of_runs);
    int stopp = 0;
    do {
        std::deque<int> _c_temp = gr_cleaned;
        bs = CUP_search(_c_temp, gr_cleaned, number_of_runs);
        ++stopp;
        if (stopp == paras->iterative_stopper)
            break;
    } while (gr_cleaned.size() > _c_.size());
    return bs;
}

bool oslom_net_evaluate::insert_the_best() {
    int Nstar = dim - cgroup.size();
    int nn = neighs.size();
    int kout_g = ktot_cgroup - kin_cgroup;
    int tm = oneM - ktot_cgroup;

    double lowest_r;
    int benode;
    neighs.best_node(benode, lowest_r, kout_g, Nstar, nn, tm);

    if (benode == -1)
        return false;

    insert_cgroup(benode);

    return true;
}

void oslom_net_evaluate::insertion(int changendi) {
    for (int i = 0; i < changendi; i++)
        insert_the_best();
}

double oslom_net_evaluate::group_inflation(const std::deque<int> &_c_, std::deque<int> &gr_cleaned,
                                           int number_of_runs) {
    if (number_of_runs == -1)
        number_of_runs = paras->inflate_runs;
//    std::cout << "group_inflation" << std::endl;
    /* preliminary check */
    double bscore = CUP_iterative(_c_, gr_cleaned);

    if (!gr_cleaned.empty()) {
        return bscore;
    }
    /* preliminary check */

    //cout<<"group inflating... "<<endl;
    WeightedTabdeg _c_tab_c;
    WeightedTabdeg _c_tab_n;
    int kin_cgroup_c;
    int ktot_cgroup_c;

    initialize_for_evaluation(_c_, _c_tab_c, _c_tab_n, kin_cgroup_c, ktot_cgroup_c);

    WeightedTabdeg previous_tab_c;
    WeightedTabdeg previous_tab_n;
    int kin_cgroup_prev;
    int ktot_cgroup_prev;

    int stopper = 0;
    while (true) {
        cgroup._set_(_c_tab_c);
        neighs._set_(_c_tab_n);
        kin_cgroup = kin_cgroup_c;
        ktot_cgroup = ktot_cgroup_c;

        int changendi = lower_bound(changendi_cum.begin(), changendi_cum.end(), ran4()) -
                        changendi_cum.begin() + 1;
        changendi = std::min(changendi, neighs.size());
        insertion(changendi);

        if (cgroup.size() == dim)
            return 1;

        /*here it make a CUP_search using c_group with the nodes added*/
        initialize_for_evaluation(previous_tab_c, previous_tab_n, kin_cgroup_prev,
                                  ktot_cgroup_prev);
        CUP_runs(previous_tab_c, previous_tab_n, kin_cgroup_prev, ktot_cgroup_prev,
                 gr_cleaned, false, number_of_runs);

        if (!gr_cleaned.empty()) {
            /*the first clean up passed. now it makes the CUP_check*/
            initialize_for_evaluation(gr_cleaned, previous_tab_c, previous_tab_n,
                                      kin_cgroup_prev,
                                      ktot_cgroup_prev);
            bscore = CUP_runs(previous_tab_c, previous_tab_n, kin_cgroup_prev,
                              ktot_cgroup_prev, gr_cleaned,
                              true, paras->clean_up_runs);

            //cout<<"exiting... "<<gr_cleaned.size()<<endl;
            if (!gr_cleaned.empty()) {
                return bscore;
            }
        }
        ++stopper;
        if (stopper == paras->inflate_stopper)
            break;
    }
    //cout<<"bad group"<<endl;
    return 1;
}

void oslom_net_evaluate::initialize_for_evaluation(WeightedTabdeg &previous_tab_c,
                                                   WeightedTabdeg &previous_tab_n,
                                                   int &kin_cgroup_prev,
                                                   int &ktot_cgroup_prev) {
    int Nstar = dim - cgroup.size();
    int nn = neighs.size();
    int kout_g = ktot_cgroup - kin_cgroup;
    int tm = oneM - ktot_cgroup;

    previous_tab_c.set_and_update_group(Nstar, nn, kout_g, tm, cgroup);
    previous_tab_n.set_and_update_neighs(Nstar, nn, kout_g, tm, neighs);
    kin_cgroup_prev = kin_cgroup;
    ktot_cgroup_prev = ktot_cgroup;
}

void oslom_net_evaluate::initialize_for_evaluation(const std::deque<int> &_c_,
                                                   WeightedTabdeg &previous_tab_c,
                                                   WeightedTabdeg &previous_tab_n,
                                                   int &kin_cgroup_prev,
                                                   int &ktot_cgroup_prev) {
    set_cgroup_and_neighs(_c_);

    int Nstar = dim - cgroup.size();
    int nn = neighs.size();
    int kout_g = ktot_cgroup - kin_cgroup;
    int tm = oneM - ktot_cgroup;

    previous_tab_c.set_and_update_group(Nstar, nn, kout_g, tm, cgroup);
    previous_tab_n.set_and_update_neighs(Nstar, nn, kout_g, tm, neighs);
    kin_cgroup_prev = kin_cgroup;
    ktot_cgroup_prev = ktot_cgroup;
}

/**
 * Initialize the data structures for the group and its neighbors
 * @param G
 */
void oslom_net_evaluate::set_cgroup_and_neighs(const std::deque<int> &G) {
    kin_cgroup = 0;
    ktot_cgroup = 0;
    cgroup.clear();
    neighs.clear();

    for (int i : G)
        insert_cgroup(i);
}

oslom_net_evaluate::oslom_net_evaluate(std::string a) : oslomnet_louvain() {
    set_graph(a);
    set_maxbord();
    set_changendi_cum();
}

oslom_net_evaluate::oslom_net_evaluate(std::map<int, std::map<int, std::pair<int, double> > > &A)
: oslomnet_louvain() {
    set_graph(A);
    set_maxbord();
    set_changendi_cum();
}

oslom_net_evaluate::oslom_net_evaluate(std::deque<std::deque<int> > &b,
                                       std::deque<std::deque<std::pair<int, double> > > &c,
std::deque<int> &d) : oslomnet_louvain() {
    set_graph(b, c, d);
    set_maxbord();
    set_changendi_cum();
}