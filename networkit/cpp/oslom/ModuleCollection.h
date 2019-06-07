#pragma once

#include <deque>
#include <set>
#include <map>

#include "utils/OslomRandom.h"
#include "utils/Histograms.h"
#include "utils/DequeNumeric.h"
#include "utils/Types.h"
#include "Parameters.h"

class ModuleCollection {

    /* all the labels refers to the index in int_matrix modules */

public:

    explicit ModuleCollection(int d);

    ~ModuleCollection() = default;

    int size() const;

    bool insert(std::deque<int> &c, double bs, int &new_name);

    bool insert(std::deque<int> &c, double bs);

    bool erase(int);

    void print(std::ostream &outt, std::deque<int> &netlabels, bool);

    void fill_gaps();

    void put_gaps();

    void homeless(std::deque<int> &h);

    int coverage();

    int effective_groups();

    void set_partition(std::deque<std::deque<int> > &A);

    void set_partition(std::deque<std::deque<int> > &A, std::deque<double> &b);

    void compute_inclusions();

    void erase_included();

    bool almost_equal(int module_id, std::deque<int> &smaller);

    void compact();

    void sort_modules(std::deque<int> &);

    void merge(std::deque<int> &c);

    /*************************** DATA ***************************/

    std::deque<std::set<int> > memberships;
    IntMatrix modules;
    std::map<int, double> module_bs;   /* it maps the module id into the b-score */

    /***********************************************************/

private:
    Parameters *paras;

    void _set_(int dim);

    bool check_already(const std::deque<int> &c);

    bool erase_first_shell(std::map<int, std::deque<int> > &erase_net);

    bool egomodules_to_merge(std::deque<int> &egom, std::deque<int> &smaller);
};




