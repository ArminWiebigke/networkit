#include "WeightedTabdeg.h"
#include "Stochastics.h"

void WeightedTabdeg::clear() {
    lab_facts.clear();
    fitness_lab.clear();
}

/**
 * Insert an element. If it was already inserted, edit the element instead.
 * @param a
 * @param kp
 * @param kt
 * @param mtlw
 * @param fit
 */
void WeightedTabdeg::edinsert(int a, int kp, int kt, double mtlw, double fit) {
    // this function inserts element a (or edit it if it was already inserted)
    erase(a);
    auto fiit = fitness_lab.insert(std::make_pair(fit, a));
    facts F(kp, mtlw, fiit, kt);
    lab_facts.insert(std::make_pair(a, F));
}

/**
 * Erase an element if it exists.
 * @param a The element to earase
 * @return True iff the element exists
 */
bool WeightedTabdeg::erase(int a) {
    auto itm = lab_facts.find(a);
    if (itm != lab_facts.end()) {
        fitness_lab.erase(itm->second.fitness_iterator);
        lab_facts.erase(itm);
        return true;
    }
    return false;
}

bool WeightedTabdeg::is_internal(int a) {
    auto itm = lab_facts.find(a);
    if (itm == lab_facts.end())
        return false;
    return true;
}

void WeightedTabdeg::set_deque(std::deque<int> &vv) {
    vv.clear();
    for (auto &lab_fact : lab_facts)
        vv.push_back(lab_fact.first);
}

void WeightedTabdeg::print_nodes(std::ostream &outb, std::deque<int> &lab_id) {
    std::cout << "printing nodes:.. (lab intk mtlw fitness degree) " << size() << std::endl;
    for (auto &lab_fact : lab_facts)
        std::cout << lab_id[lab_fact.first] << " " << lab_fact.second.internal_degree << " "
                  << lab_fact.second.minus_log_total_wr << " "
                  << (lab_fact.second.fitness_iterator)->first << " " << lab_fact.second.degree
                  << std::endl;
}

int WeightedTabdeg::worst_node(int &lab, double &worst_fitness, int kout_g, int Nstar,
                                int nneighs, int tm) {
    //std::cout<<"worst_node fitness - lab - (cgroup)"<<endl;
    //prints(fitness_lab);
    lab = -1;
    worst_fitness = -1;
    auto bit = fitness_lab.end();
    if (bit == fitness_lab.begin())
        return -1;

    int stopper = 0;
    while (bit != fitness_lab.begin()) {
        bit--;
        auto itm = lab_facts.find(bit->second);
        double F = Stochastics::compute_global_fitness_randomized(
                itm->second.internal_degree,
                kout_g + 2 * itm->second.internal_degree - itm->second.degree,
                tm + itm->second.degree,
                itm->second.degree,
                itm->second.minus_log_total_wr,
                nneighs + 1, Nstar + 1);

        if (F > worst_fitness) {
            worst_fitness = F;
            lab = itm->first;
        }

        stopper++;
        if (stopper == num_up_to)
            break;
    }
    return 0;
}

int WeightedTabdeg::best_node(int &lab, double &best_fitness, int kout_g, int Nstar,
                               int nneighs, int tm) {
    // I can try to compute the fitness here
    /*std::cout<<"NE BEST NODE "<<endl;
    std::cout<<"fitness_lab  "<<endl;
    prints(fitness_lab);*/
    lab = -1;
    best_fitness = 1;

    auto bit = fitness_lab.begin();
    if (bit == fitness_lab.end()) {
        return -1;
    }

    int stopper = 0;
    while (bit != fitness_lab.end()) {
        auto itm = lab_facts.find(bit->second);
        double F = Stochastics::compute_global_fitness_randomized(
                itm->second.internal_degree, kout_g,
                tm, itm->second.degree,
                itm->second.minus_log_total_wr,
                nneighs, Nstar);
        //std::cout<<itm->first<<" "<<F<<" ... node-fit"<<endl;
        if (F < best_fitness) {
            best_fitness = F;
            lab = itm->first;
        }

        stopper++;
        if (stopper == num_up_to)
            break;

        bit++;
    }

    return 0;
}

void WeightedTabdeg::_set_(WeightedTabdeg &one) {
    clear();
    for (auto &lab_fact : one.lab_facts)
        edinsert(lab_fact.first, lab_fact.second.internal_degree, lab_fact.second.degree,
                 lab_fact.second.minus_log_total_wr, (lab_fact.second.fitness_iterator)->first);
}

bool WeightedTabdeg::update_group(int a, int delta_degree, double delta_mtlw, int nstar,
                                   int nn, int kout_g, int tm, int kt,
                                   std::deque<int> &to_be_erased) {
    // this function is to change the internal degree and mtlw of a certain node (to insert it or erase if necessary)
    auto itm = lab_facts.find(a);
    if (itm == lab_facts.end())
        return false;

    itm->second.minus_log_total_wr += delta_mtlw;
    itm->second.internal_degree += delta_degree;

    if (itm->second.internal_degree == 0 && size() > 1) {
        to_be_erased.push_back(a);
        return true;
    }

    //std::cout<<"UPdating... group "<<a<<endl;
    double fit = Stochastics::compute_global_fitness_ofive(
            itm->second.internal_degree,
            kout_g + 2 * itm->second.internal_degree -
            itm->second.degree,
            tm + itm->second.degree, itm->second.degree,
            itm->second.minus_log_total_wr, nn + 1,
            nstar + 1);

    fitness_lab.erase(itm->second.fitness_iterator);
    auto fiit = fitness_lab.insert(std::make_pair(fit, a));
    itm->second.fitness_iterator = fiit;

    return true;
}

bool WeightedTabdeg::update_neighs(int a, int delta_degree, double delta_mtlw, int nstar,
                                    int kout_g, int tm, int kt) {
    // this function is to change the internal degree and mtlw of a certain node (to insert it or erase if necessary)
    //std::cout<<"UPdating... neighs "<<a<<" "<<kt<<endl;
    auto itm = lab_facts.find(a);
    if (itm == lab_facts.end()) {
        edinsert(a, 0, kt, 0, 1);
        itm = lab_facts.find(a);
    }

    itm->second.internal_degree += delta_degree;
    if (itm->second.internal_degree == 0) {
        //std::cout<<"erased from neigh update "<<a<<std::endl;
        erase(a);
        return true;
    }

    itm->second.minus_log_total_wr += delta_mtlw;
    double fit = Stochastics::compute_global_fitness_ofive(
            itm->second.internal_degree, kout_g, tm,
            itm->second.degree,
            itm->second.minus_log_total_wr, size(),
            nstar);

    fitness_lab.erase(itm->second.fitness_iterator);
    auto fiit = fitness_lab.insert(std::make_pair(fit, a));
    itm->second.fitness_iterator = fiit;
    return true;
}

void WeightedTabdeg::set_and_update_group(int nstar, int nn, int kout_g, int tm,
                                           WeightedTabdeg &one) {
    /*this function is to set and update the fitnesses of all the nodes in cgroup*/
    clear();
    for (auto &lab_fact : one.lab_facts) {
        double fit = Stochastics::compute_global_fitness_ofive(
                lab_fact.second.internal_degree,
                kout_g + 2 * lab_fact.second.internal_degree - lab_fact.second.degree,
                tm + lab_fact.second.degree,
                lab_fact.second.degree,
                lab_fact.second.minus_log_total_wr, nn + 1,
                nstar + 1);
        edinsert(lab_fact.first, lab_fact.second.internal_degree, lab_fact.second.degree,
                 lab_fact.second.minus_log_total_wr, fit);
    }
}

void WeightedTabdeg::set_and_update_neighs(int nstar, int nn, int kout_g, int tm,
                                            WeightedTabdeg &one) {
    /*this function is to set and update the fitnesses of all the nodes in neighs*/
    clear();
    for (auto &lab_fact : one.lab_facts) {
        double fit = Stochastics::compute_global_fitness_ofive(
                lab_fact.second.internal_degree, kout_g, tm,
                lab_fact.second.degree,
                lab_fact.second.minus_log_total_wr, nn,
                nstar);
        edinsert(lab_fact.first, lab_fact.second.internal_degree, lab_fact.second.degree,
                 lab_fact.second.minus_log_total_wr, fit);
    }
}

