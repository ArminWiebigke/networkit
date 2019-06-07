#include <fstream>

#include "StaticNetwork.h"
#include "utils/DequeNumeric.h"
#include "utils/Cast.h"
#include "utils/Histograms.h"

StaticNetwork::vertex::vertex(int b, int c, int preall) {
    id_num = b;
    links = new Wsarray(preall);
}

StaticNetwork::vertex::~vertex() {
    delete links;
    links = nullptr;
}

void StaticNetwork::vertex::kplus_global_and_quick(std::deque<int> &a, int &stubs_in,
                                                   double &strength_in) {
    // a must be sorted, otherwise it does not work
    // this function has never been used (so, it should be checked)
    stubs_in = 0;
    strength_in = 0;

    if (int(a.size()) > stub_number) {
        // in this case there must be a loop on links because a is longer
        for (int i = 0; i < links->size(); i++) {
            if (std::binary_search(a.begin(), a.end(), links->l[i])) {
                stubs_in += links->w[i].first;
                strength_in += links->w[i].second;
            }
        }
    } else {
        for (int i : a) {
            std::pair<int, double> A = links->posweightof(i).second;
            stubs_in += A.first;
            strength_in += A.second;
        }
    }
}

int StaticNetwork::vertex::kplus_m(const std::deque<int> &a) {
    // computes the internal degree of the vertex respect with a
    int s = 0;
    //double f=0;
    for (int i : a) {
        std::pair<int, double> A = links->posweightof(i).second;
        s += A.first;
        //f+=A.second;
    }
    return s;
}

double StaticNetwork::vertex::kplus_w(const std::deque<int> &a) {
    // computes the internal degree of the vertex respect with a
    double s = 0;
    //double f=0;
    for (int i : a) {
        std::pair<int, double> A = links->posweightof(i).second;
        s += A.second;
        //f+=A.second;
    }
    return s;
}

int StaticNetwork::vertex::kplus_m(const std::set<int> &a) {
    // computes the internal degree of the vertex respect with a (a is supposed to be sorted)
    int s = 0;
    //double f=0;
    for (int i = 0; i < links->size(); i++) {
        if (a.find(links->l[i]) != a.end()) {
            s += links->w[i].first;
            //f+=links->w[i].second;
        }
    }
    return s;
}

StaticNetwork::~StaticNetwork() {
    clear();
}

void StaticNetwork::clear() {
    for (auto &vertice : vertices) {
        delete vertice;
        vertice = nullptr;
    }
    vertices.clear();
    dim = 0;
    oneM = 0;
}

// node_id -> (neighbor_id -> (edge_id, edge_weight))
void StaticNetwork::set_graph(std::map<int, std::map<int, std::pair<int, double> > > &A) {
    // this maps the id into the usual stuff neighbors-weights
    std::deque<std::deque<int>> link_per_node;
    std::deque<std::deque<std::pair<int, double> > > weights_per_node;
    std::deque<int> label_rows;

    for (auto &itm : A) {
        label_rows.push_back(itm.first);
        std::deque<int> n;
        std::deque<std::pair<int, double>> w;
        for (auto &itm2 : itm.second)
            if (itm2.second.first > 0) {
                n.push_back(itm2.first);
                w.push_back(itm2.second);
            }
        link_per_node.push_back(n);
        weights_per_node.push_back(w);
    }
    /*
    prints(label_rows);
    printm(link_per_node);
    printm(weights_per_node);*/
    set_graph(link_per_node, weights_per_node, label_rows);
}

//   all this stuff here should be improved
bool StaticNetwork::set_graph(const std::string &file_name) {
    clear();

    char b[file_name.size() + 1];
    cast_string_to_char(file_name, b);

    std::map<int, int> newlabels;
    std::deque<int> link_i;

    bool good_file = true;
    {
        int label = 0;
        std::ifstream inb(b);
        std::string ins;
        while (getline(inb, ins)) {
            if (!ins.empty() && ins[0] != '#') {
                std::deque<double> ds;
                cast_string_to_doubles(ins, ds);

                if (ds.size() < 2) {
                    std::cerr << "From file " << file_name << ": string not readable " << ins
                              << " " << std::endl;
                    good_file = false;
                    break;
                } else {
                    int innum1 = cast_int(ds[0]);
                    int innum2 = cast_int(ds[1]);

                    auto itf = newlabels.find(innum1);
                    if (itf == newlabels.end()) {
                        newlabels.insert(std::make_pair(innum1, label++));
                        link_i.push_back(1);
                    } else
                        link_i[itf->second]++;

                    itf = newlabels.find(innum2);
                    if (itf == newlabels.end()) {
                        newlabels.insert(std::make_pair(innum2, label++));
                        link_i.push_back(1);
                    } else
                        link_i[itf->second]++;
                }
            }
        }
    }

    dim = newlabels.size();

    for (int i = 0; i < dim; i++)
        vertices.push_back(new vertex(0, 0, link_i[i]));

    for (auto &newlabel : newlabels)
        vertices[newlabel.second]->id_num = newlabel.first;

    if (good_file) {
        std::ifstream inb(b);
        std::string ins;
        while (getline(inb, ins))
            if (!ins.empty() && ins[0] != '#') {
                std::deque<double> ds;
                cast_string_to_doubles(ins, ds);

                int innum1 = cast_int(ds[0]);
                int innum2 = cast_int(ds[1]);

                double w = 1;
                int multiple_l = 1;

                if (ds.size() >= 4) {
                    if (paras->weighted == false) {
                        if (ds[2] > 0.99) {
                            multiple_l = cast_int(ds[2]);
                        }
                    }

                    if (paras->weighted == true) {
                        if (ds[2] > 0) {
                            w = ds[2];
                        } else {
                            std::cerr << "error: not positive weights" << std::endl;
                            return false;
                        }

                        if (ds[3] > 0.99) {
                            multiple_l = cast_int(ds[3]);
                        }
                    }
                }

                if (ds.size() == 3) {

                    if (paras->weighted) {
                        if (ds[2] > 0)
                            w = ds[2];
                        else {
                            std::cerr << "error: not positive weights" << std::endl;
                            return false;
                        }
                    } else {
                        if (ds[2] > 0.99)
                            multiple_l = cast_int(ds[2]);
                    }
                }

                int new1 = newlabels[innum1];
                int new2 = newlabels[innum2];
                if (new1 != new2 && w > 0) {        // no self loops!
                    vertices[new1]->links->push_back(new2, multiple_l, w);
                    vertices[new2]->links->push_back(new1, multiple_l, w);
                }
            }

        oneM = 0;

        for (int i = 0; i < dim; i++) {
            vertices[i]->links->freeze();
            int stub_number_i = 0;
            double strength_i = 0;
            for (int j = 0; j < vertices[i]->links->size(); j++) {
                stub_number_i += vertices[i]->links->w[j].first;
                strength_i += vertices[i]->links->w[j].second;
            }

            vertices[i]->stub_number = stub_number_i;
            vertices[i]->strength = strength_i;
            oneM += stub_number_i;
        }
    } else
        std::cerr << "File corrupted" << std::endl;

    //draw("before_set");
    if (paras->weighted)
        set_proper_weights();

    return good_file;
}

void StaticNetwork::set_graph(std::deque<std::deque<int>> &link_per_node,
                              std::deque<std::deque<std::pair<int, double>>> &weights_per_node,
                              std::deque<int> &label_rows) {
    clear();
    // there is no check between label_rows and link per node but they need to have the same labels
    // link_per_node and weights_per_node are the list of links and weights. label_rows[i] is the label corresponding to row i
    std::map<int, int> newlabels;        // this maps the old labels with the new one
    for (int &label_row : label_rows)
        newlabels.insert(std::make_pair(label_row, newlabels.size()));

    dim = newlabels.size();

    for (int i = 0; i < dim; i++)
        vertices.push_back(new vertex(0, 0, link_per_node[i].size()));

    for (auto &newlabel : newlabels)
        vertices[newlabel.second]->id_num = newlabel.first;

    for (int i = 0; i < int(link_per_node.size()); i++) {
        for (int j = 0; j < int(link_per_node[i].size()); j++) {
            int new2 = newlabels[link_per_node[i][j]];
            vertices[i]->links->push_back(new2, weights_per_node[i][j].first,
                                          weights_per_node[i][j].second);
        }
    }

    oneM = 0;

    for (int i = 0; i < dim; i++) {
        vertices[i]->links->freeze();
        int stub_number_i = 0;
        double strength_i = 0;
        for (int j = 0; j < vertices[i]->links->size(); j++) {
            stub_number_i += vertices[i]->links->w[j].first;
            strength_i += vertices[i]->links->w[j].second;
        }

        vertices[i]->stub_number = stub_number_i;
        vertices[i]->strength = strength_i;
        oneM += stub_number_i;
    }

    if (paras->weighted)
        set_proper_weights();
}

int StaticNetwork::kin_m(const std::deque<int> &seq) {
    if (seq.size() > double(oneM) / dim) {
        std::set<int> H;
        deque_numeric::deque_to_set(seq, H);
        return kin_m(H);
    }

    int k = 0;
    for (int i = 0; i < int(seq.size()); i++)
        k += vertices[seq[i]]->kplus_m(seq);

    return k;
}

int StaticNetwork::ktot_m(const std::deque<int> &seq) {
    int k = 0;
    for (int i : seq)
        k += vertices[i]->stub_number;
    return k;
}

int StaticNetwork::ktot_m(const std::set<int> &s) {
    int k = 0;
    for (int it : s)
        k += vertices[it]->stub_number;
    return k;
}

int StaticNetwork::kin_m(const std::set<int> &s) {
    int k = 0;
    for (int it : s)
        k += vertices[it]->kplus_m(s);

    return k;
}

int StaticNetwork::draw(std::string file_name) {
    int h = file_name.size();

    char b[h + 1];
    for (int i = 0; i < h; i++)
        b[i] = file_name[i];
    b[h] = '\0';

    std::ofstream graph_out(b);

    if (paras->weighted) {
        for (auto &vertice : vertices)
            for (int j = 0; j < vertice->links->size(); j++)
                if (vertice->id_num <= vertices[vertice->links->l[j]]->id_num)
                    graph_out << vertice->id_num << "\t"
                              << vertices[vertice->links->l[j]]->id_num << "\t"
                              << vertice->original_weights[j] << "\t"
                              << vertice->links->w[j].first << "\t" << std::endl;
    } else {
        for (auto &vertice : vertices)
            for (int j = 0; j < vertice->links->size(); j++)
                if (vertice->id_num <= vertices[vertice->links->l[j]]->id_num)
                    graph_out << vertice->id_num << "\t"
                              << vertices[vertice->links->l[j]]->id_num << "\t"
                              << vertice->links->w[j].first << std::endl;
    }
    return 0;
}

void StaticNetwork::get_id_label(std::map<int, int> &a) {

    for (int i = 0; i < dim; i++)
        a.insert(std::make_pair(vertices[i]->id_num, i));
}

void StaticNetwork::deque_id(std::deque<int> &a) {
    for (int &i : a)
        i = vertices[i]->id_num;
}

void StaticNetwork::print_id(const std::deque<int> &a, std::ostream &pout) {
    for (int i : a)
        pout << vertices[i]->id_num << "\t";
    pout << std::endl;
}

void StaticNetwork::print_id(const std::set<int> &a, std::ostream &pout) {

    for (int its : a)
        pout << vertices[its]->id_num << "\t";
    pout << std::endl;
}

void StaticNetwork::print_id(const std::deque<std::deque<int>> &a, std::ostream &pout) {
    for (const auto &i : a)
        print_id(i, pout);
}

void StaticNetwork::print_id(const std::deque<std::set<int>> &a,
                             std::ostream &pout) {
    for (const auto &i : a)
        print_id(i, pout);
}

int StaticNetwork::translate_anyway(std::deque<std::deque<int>> &ten) {
    std::map<int, int> A;
    get_id_label(A);
    std::deque<std::deque<int>> ten2;

    for (auto &i : ten) {
        std::deque<int> ff;
        for (int j : i) {
            auto itf = A.find(j);
            if (itf != A.end())
                ff.push_back(itf->second);
        }

        if (!ff.empty())
            ten2.push_back(ff);
    }
    ten = ten2;
    return 0;
}

int StaticNetwork::translate(std::deque<std::deque<int>> &ten) {
    std::map<int, int> A;
    get_id_label(A);
    std::deque<std::deque<int>> ten2;

    for (auto &i : ten) {
        std::deque<int> ff;
        for (int j : i) {
            auto itf = A.find(j);
            if (itf == A.end()) {
                std::cerr << "warning: the nodes in the communities are different from those ones"
                             " in the network!"
                          << std::endl;
            } else
                ff.push_back(itf->second);
        }

        if (!ff.empty())
            ten2.push_back(ff);
    }
    ten = ten2;
    return 0;
}

void StaticNetwork::set_subgraph(
        std::deque<int> &group,
        std::deque<std::deque<int>> &link_per_node,
        std::deque<std::deque<std::pair<int, double>>> &weights_per_node) {
    // in this function I'm not using id's... because I want to work with the same labels (don't want to translate)
    sort(group.begin(), group.end());

    weights_per_node.clear();

    link_per_node.clear();

    for (int i = 0; i < int(group.size()); i++) {
        int nodei = group[i];
        std::deque<int> link_i;
        std::deque<std::pair<int, double>> weight_i;

        for (int j = 0; j < vertices[nodei]->links->size(); j++)
            if (binary_search(group.begin(), group.end(), vertices[nodei]->links->l[j])) {
                link_i.push_back(vertices[nodei]->links->l[j]);
                //weight_i.push_back(make_pair(vertices[nodei]->links->w[j].first,   vertices[nodei]->original_weights[j]));
                if (paras->weighted)
                    weight_i.emplace_back(vertices[nodei]->links->w[j].first,
                                          vertices[nodei]->original_weights[j]);
                else
                    weight_i.emplace_back(vertices[nodei]->links->w[j].first, 1);
                //weight_i.push_back(make_pair(vertices[nodei]->links->w[j].first, vertices[nodei]->links->w[j].second));
            }

        link_per_node.push_back(link_i);
        weights_per_node.push_back(weight_i);
    }
}

void StaticNetwork::set_proper_weights() {
    // this function id to normalize the weights in order to have the -log(prob(Weiight>w)) which is simply w[i].second / <w_ij>
    if (dim == 0) {
        //cout<<"network empty"<<endl;
        //cherr();
    } else {
        for (int i = 0; i < dim; i++) {
            vertices[i]->original_weights.clear();
            for (int j = 0; j < vertices[i]->links->size(); j++) {
                vertices[i]->original_weights.push_back(vertices[i]->links->w[j].second);
            }
        }

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < vertices[i]->links->size(); j++) {
                double w1 = (vertices[i]->strength / vertices[i]->stub_number) *
                            vertices[i]->links->w[j].first;
                double w2 = (vertices[vertices[i]->links->l[j]]->strength /
                             vertices[vertices[i]->links->l[j]]->stub_number) *
                            vertices[i]->links->w[j].first;
                vertices[i]->links->w[j].second /= 2. / (1. / w1 + 1. / w2);
            }
        }
    }
}

void StaticNetwork::set_connected_components(std::deque<std::deque<int>> &comps) {
    comps.clear();
    std::set<int> not_assigned;
    for (int i = 0; i < dim; i++)
        not_assigned.insert(i);

    while (!not_assigned.empty()) {
        int source = *not_assigned.begin();
        std::set<int> mates;
        same_component(source, mates);

        std::deque<int> ccc;
        for (int mate : mates) {
            ccc.push_back(mate);
            not_assigned.erase(mate);
        }
        comps.push_back(ccc);
    }
}

void StaticNetwork::same_component(int source, std::set<int> &already_gone) {
    already_gone.clear();
    already_gone.insert(source);

    std::deque<int> this_shell;
    this_shell.push_back(source);

    while (!this_shell.empty()) {
        std::deque<int> next_shell;
        for (int i : this_shell) {
            for (int j = 0; j < vertices[i]->links->size(); j++) {
                if (already_gone.insert(vertices[i]->links->l[j]).second)
                    next_shell.push_back(vertices[i]->links->l[j]);
            }
        }
        this_shell = next_shell;
    }
}

int StaticNetwork::propagate_distances(std::deque<int> &new_shell, std::set<int> &already_gone,
                                       std::deque<std::pair<int, int>> &distances_node,
                                       int shell, std::deque<double> &ML, int &reached, int step) {
    shell++;
    std::deque<int> next_shell;
    for (int i : new_shell)
        for (int j = 0; j < vertices[i]->links->size(); j++) {
            if (already_gone.insert(vertices[i]->links->l[j]).second) {
                distances_node.emplace_back(shell, vertices[i]->links->l[j]);
                next_shell.push_back(vertices[i]->links->l[j]);
            }
        }
    /*
    cout<<"new shell "<<shell<<endl;
    print_id(next_shell, cout);
    prints(ML);
    */
    if (!next_shell.empty()) {
        if (shell >= int(ML.size()))
            ML.push_back(dim * step);
        reached += next_shell.size();
        ML[shell] += reached;
        return propagate_distances(next_shell, already_gone, distances_node, shell, ML, reached,
                                   step);
    }
    return 0;
}

int StaticNetwork::set_upper_network(
        std::map<int, std::map<int, std::pair<int, double> > > &neigh_weight_f,
        ModuleCollection &Mcoll) {
    // loop on all the edges of the network...
    neigh_weight_f.clear();

    if (Mcoll.size() == 0)
        return 0;

    std::map<int, std::map<int, std::pair<double, double> > > neigh_weight_s;

    for (auto & module_b : Mcoll.module_bs) {
        std::map<int, std::pair<double, double> > neigh_weight;
        std::map<int, std::pair<int, double> > ooo;
        neigh_weight_s.insert(make_pair(module_b.first, neigh_weight));
        neigh_weight_f.insert(make_pair(module_b.first, ooo));
    }

    for (int i = 0; i < dim; i++) {
        std::set<int> &mem1 = Mcoll.memberships[i];
        for (int j = 0; j < vertices[i]->links->size(); j++) {
            int &neigh = vertices[i]->links->l[j];
            std::set<int> &mem2 = Mcoll.memberships[neigh];
            double denominator = mem1.size() * mem2.size();
            // I add a link between all different modules

            if (paras->weighted) {
                for (int itk : mem1)
                    for (int itkk : mem2)
                        if (itk != itkk)
                            int_histogram(itkk, neigh_weight_s[itk],
                                          double(vertices[i]->links->w[j].first) / denominator,
                                          vertices[i]->original_weights[j] / denominator);
            } else {
                for (int itk : mem1)
                    for (int itkk : mem2)
                        if (itk != itkk)
                            int_histogram(itkk, neigh_weight_s[itk],
                                          double(vertices[i]->links->w[j].first) / denominator,
                                          double(vertices[i]->links->w[j].first) / denominator);
            }
        }
    }

    for (auto &neigh_weight_ : neigh_weight_s) {
        for (auto itm_ = neigh_weight_.second.begin(); itm_ != neigh_weight_.second.end(); itm_++) {
            int module_a = itm_->first;
            int module_b = neigh_weight_.first;
            if (module_a < module_b) { // So we don't count twice?
                int intml = cast_int(itm_->second.first);
                if (intml > 0) {
                    neigh_weight_f[module_b].insert(
                            std::make_pair(module_a, std::make_pair(intml, itm_->second.second)));
                    neigh_weight_f[module_a].insert(
                            std::make_pair(module_b, std::make_pair(intml, itm_->second.second)));
                }
            }
        }
    }
    return 0;
}

int StaticNetwork::draw_with_weight_probability(std::string file_name) {
    int h = file_name.size();

    char b[h + 1];
    for (int i = 0; i < h; i++)
        b[i] = file_name[i];
    b[h] = '\0';

    std::ofstream graph_out(b);
    if (paras->weighted) {
        for (auto & vertice : vertices)
            for (int j = 0; j < vertice->links->size(); j++)
                graph_out << vertice->id_num << "\t"
                          << vertices[vertice->links->l[j]]->id_num << "\t"
                          << vertice->original_weights[j] << "\t"
                          << vertice->links->w[j].first << std::endl;
    }
    return 0;
}


void StaticNetwork::print_degree_of_homeless(std::deque<int> &homel, std::ostream &outt) {
    std::deque<int> degree_of_homeless;
    for (int i : homel)
        degree_of_homeless.push_back(vertices[i]->stub_number);

    outt << "average degree of homeless nodes: " << average_func(degree_of_homeless)
         << " dev: " << sqrt(variance_func(degree_of_homeless)) << std::endl;


}

StaticNetwork::StaticNetwork() : dim(0), oneM(0) {
    paras = Parameters::get_instance();
}
