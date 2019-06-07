
#include <algorithm>

#include "ModuleCollection.h"

ModuleCollection::ModuleCollection(int dim) {
    _set_(dim);
    paras = Parameters::get_instance();
}

void ModuleCollection::_set_(int dim) {
    std::set<int> first;
    for (int i = 0; i < dim; i++)
        memberships.push_back(first);
}

bool ModuleCollection::insert(std::deque<int> &c, double bs) {
    int new_name;
    return insert(c, bs, new_name);
}

bool ModuleCollection::insert(std::deque<int> &c, double bs, int &new_name) {
    if (bs == 0)
        bs = ran4() * 1e-100;

    std::sort(c.begin(), c.end());
    new_name = -1;

    if (check_already(c)) {
        new_name = modules.size();
        for (int i : c)
            memberships[i].insert(new_name);
        modules.push_back(c);
        module_bs[new_name] = bs;
        return true;
    }
    return false;
}

bool ModuleCollection::erase(int a) {
    // it erases module a
    if (module_bs.find(a) == module_bs.end()) // it only erases not empty modules
        return false;

    std::deque<int> &nodes_a = modules[a];

    for (int i : nodes_a)
        memberships[i].erase(a);
    modules[a].clear();
    module_bs.erase(a);
    return true;
}

void ModuleCollection::print(std::ostream &outt, std::deque<int> &netlabels, bool not_homeless) {
    int nmod = 0;
    for (auto &module_b : module_bs)
        if (!not_homeless || modules[module_b.first].size() > 1) {
            nmod++;
            std::deque<int> &module_nodes = modules[module_b.first];
            outt << "#module " << module_b.first << " size: " << modules[module_b.first].size()
                 << " bs: " << module_bs[module_b.first] << std::endl;

            std::deque<int> labseq;
            for (int module_node : module_nodes) {
                labseq.push_back(netlabels[module_node]);
            }
            std::sort(labseq.begin(), labseq.end());
            for (int i : labseq) {
                outt << i << " ";
            }
            outt << std::endl;
        }
}

void ModuleCollection::fill_gaps() {
    for (int i = 0; i < int(memberships.size()); i++)
        if (memberships[i].empty()) {
            std::deque<int> new_d;
            new_d.push_back(i);
            insert(new_d, 1.);
        }
}

void ModuleCollection::put_gaps() {
    std::deque<int> to_erase;
    for (int i = 0; i < int(modules.size()); i++) {
        if (modules[i].size() == 1)
            to_erase.push_back(i);
    }
    for (int i : to_erase)
        erase(i);
}

void ModuleCollection::homeless(std::deque<int> &h) {
    h.clear();
    for (int i = 0; i < int(memberships.size()); i++)
        if (memberships[i].empty())
            h.push_back(i);
    for (auto &module : modules) {
        if (module.size() == 1)
            h.push_back(module[0]);
    }
    std::sort(h.begin(), h.end());
}

int ModuleCollection::coverage() {
    // this function returns the number of nodes which are covered by at least one module
    int cov = 0;
    for (auto &membership : memberships)
        if (!membership.empty())
            cov++;
    return cov;
}

int ModuleCollection::effective_groups() {
    int nmod = 0;
    for (auto &module_b : module_bs)
        if (modules[module_b.first].size() > 1)
            nmod++;
    return nmod;
}

void ModuleCollection::set_partition(std::deque<std::deque<int> > &A) {
    A.clear();
    for (auto &module_b : module_bs)
        if (modules[module_b.first].size() > 1)
            A.push_back(modules[module_b.first]);
}

void ModuleCollection::set_partition(std::deque<std::deque<int> > &A, std::deque<double> &b) {
    A.clear();
    b.clear();
    for (auto &module_b : module_bs)
        if (modules[module_b.first].size() > 1) {
            A.push_back(modules[module_b.first]);
            b.push_back(module_bs[module_b.first]);
        }
}

bool ModuleCollection::check_already(const std::deque<int> &c) {
    // returns false if the module is already present

    // it maps the index of the modules into the overlap (overlap=number of overlapping nodes)
    std::map<int, int> com_ol;
    for (int i : c) {
        for (int itj : memberships[i])
            int_histogram(itj, com_ol);
    }

    for (auto &itm : com_ol) {
        if (itm.second == int(c.size()) &&
            itm.second == int(modules[itm.first].size()))
            return false;
    }
    return true;
}

void ModuleCollection::compute_inclusions() {
    put_gaps();
    erase_included();
    compact();
}

void ModuleCollection::erase_included() {
    std::map<int, std::deque<int> > erase_net;
    for (auto &module_b : module_bs) {
        std::deque<int> smaller;
        almost_equal(module_b.first, smaller);
        erase_net[module_b.first] = smaller;
    }
    while (true) {
        if (!erase_first_shell(erase_net))
            break;
    }
}

bool ModuleCollection::erase_first_shell(std::map<int, std::deque<int> > &erase_net) {
    bool again = false;
    std::set<int> roots;

    for (auto &module_b : module_bs)
        roots.insert(module_b.first);

    for (auto &itm : erase_net) {
        std::deque<int> &smaller = itm.second;
        for (int i : smaller)
            roots.erase(i);
    }
    //cout<<"roots:"<<endl;
    //prints(roots);

    for (int root : roots) {
        std::deque<int> &smaller = erase_net[root];
        for (int i : smaller) {
            if (module_bs.find(i) != module_bs.end()) {
                erase(i);
                erase_net.erase(i);
                again = true;
            }
        }
    }
    return again;
}

bool ModuleCollection::almost_equal(int module_id, std::deque<int> &smaller) {
    // c is the module you want to know about
    // smaller is set to contain the module ids contained by module_id
    smaller.clear();
    std::deque<int> &c = modules[module_id];
    std::map<int, int> com_ol; // maps the index of the modules into the overlap (overlap=number of overlapping nodes)
    for (int i : c) {
        for (int itj : memberships[i])
            int_histogram(itj, com_ol);
    }

    for (auto &itm : com_ol)
        if (itm.first != module_id && modules[itm.first].size() <= c.size()) {
            const count &other_size = modules[itm.first].size();
            if (double(itm.second) / other_size >=
                paras->coverage_inclusion_module_collection) {
                if (c.size() > other_size)
                    smaller.push_back(itm.first);
                else if (c.size() == other_size &&
                         module_bs[module_id] < module_bs[itm.first])
                    smaller.push_back(itm.first);
            }
        }
    return true;
}

void ModuleCollection::compact() {
    /* this function is used to have continuos ids */
    put_gaps();
    std::map<int, int> from_old_index_to_new;
    {
        std::deque<std::deque<int> > modules2;
        std::map<int, double> module_bs2;
        for (auto &module_b : module_bs) {
            from_old_index_to_new.insert(
                    std::make_pair(module_b.first, from_old_index_to_new.size()));
            modules2.push_back(modules[module_b.first]);
            module_bs2[from_old_index_to_new.size() - 1] = module_b.second;
        }

        modules = modules2;
        module_bs = module_bs2;
    }
    for (auto &membership : memberships) {
        std::set<int> first;
        for (int its : membership)
            first.insert(from_old_index_to_new[its]);
        membership = first;
    }
}

void ModuleCollection::sort_modules(std::deque<int> &module_order) {
    module_order.clear();
    std::multimap<double, int> rank_id;
    /* modules are sorted from the biggest to the smallest. if they have equal size, we look at the score */
    for (auto &module_b : module_bs) {
        //cout<<modules[itm->first].size()<<" ... "<<endl;
        rank_id.insert(
                std::make_pair(-double(modules[module_b.first].size()) + 1e-2 * module_b.second,
                               module_b.first));
    }
    //cout<<"rank_id"<<endl;
    //prints(rank_id);
    for (auto &itm : rank_id)
        module_order.push_back(itm.second);
}

bool ModuleCollection::egomodules_to_merge(std::deque<int> &egom, std::deque<int> &smaller) {
    // egom is the module you want to know about
    // smaller is set to contain the module ids to merge with egom
    smaller.clear();
    // it maps the index of the modules into the overlap (overlap=numeber of overlapping nodes)
    std::map<int, int> com_ol;
    for (int i : egom) {
        for (int itj : memberships[i])
            int_histogram(itj, com_ol);
    }
    for (auto &itm : com_ol) {
        const count &other_size = std::min(modules[itm.first].size(), egom.size());
        if (double(itm.second) / other_size >=
            paras->coverage_inclusion_module_collection)
            smaller.push_back(itm.first);
    }
    return true;
}

void ModuleCollection::merge(std::deque<int> &c) {
    std::deque<int> to_merge;
    egomodules_to_merge(c, to_merge);
    if (to_merge.empty())
        insert(c, 1.);
    else {
        for (int i : to_merge) {
            std::set<int> si;
            deque_numeric::deque_to_set_app(modules[i], si);
            deque_numeric::deque_to_set_app(c, si);
            erase(i);
            std::deque<int> to_insert;
            deque_numeric::set_to_deque(si, to_insert);
            insert(to_insert, 1.);
        }
    }
    erase_included();
}

int ModuleCollection::size() const {
    return module_bs.size();
}
