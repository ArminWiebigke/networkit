#include "LouvainOslomnet.h"
#include "utils/Histograms.h"

/**
 * Initialize each node as its own module.
 */
void oslomnet_louvain::module_initializing() {
	for (int i = 0; i < dim; i++) {
		vertex_label.push_back(i);
		vertex_order.push_back(i);
		vertex_to_check.push_back(true);
		vertex_to_check_next.push_back(false);
		oslom_module newm(vertices[i]->stub_number);
		label_module.insert(std::make_pair(i, newm));
	}
}

/**
 * Find the "favourite" module of a node. The best module is the one with the lowest fitness score,
 *
 * @param node The node we examine.
 * @param new_label Label of the new module
 * @param stubs_into_new Number of edges from the node to the new module
 * @param stubs_into_old Number of edges from the node to the old module
 */
void
oslomnet_louvain::unweighted_favorite_of(const int &node, int &new_label, int &stubs_into_new,
                                         int &stubs_into_old) {
	double min_fitness = 10;
	int cur_label = vertex_label[node];
	new_label = cur_label;
	stubs_into_new = 0;
	stubs_into_old = 0;
	std::map<int, int> M; // M is a map module_label -> kin (internal stubs)
	for (int j = 0; j < vertices[node]->links->size(); j++) // for each neighbor of the node
		int_histogram(M, vertex_label[vertices[node]->links->l[j]],
		              vertices[node]->links->w[j].first); // (insert into M), neighbor, number of edges
	for (auto &itM : M) {
		int module_label = itM.first;
		int num_edges = itM.second;
		const oslom_module &module = label_module.find(module_label)->second;
		double fitness;
		int degree = vertices[node]->stub_number;
		if (module_label != cur_label) {
			int extStubs = total_stubs - module.ktot;
			fitness = Stochastics::topological_05(num_edges,
			                                      module.kout,
			                                      extStubs,
			                                      degree);
		} else {
			// Calc fitness of current module
			stubs_into_old = num_edges; // Current internal degree
			int kout_prime = module.kout - degree + 2 * num_edges;
			int extStubs = total_stubs - module.ktot + degree;
			fitness = Stochastics::topological_05(
					num_edges,
					kout_prime,
					extStubs,
					degree);
			fitness *= 0.999;        // to break possible ties
		}
		if (fitness < min_fitness) {
			stubs_into_new = num_edges;
			min_fitness = fitness;
			new_label = module_label;
		}
	}
}

void oslomnet_louvain::weighted_favorite_of(const int &node, int &fi, int &kp, int &kop) {
	double min_fitness = 10;
	fi = vertex_label[node];
	kp = 0;
	kop = 0;
	mapip M;        // M is a map module_label -> kin - win (internal stubs, internal weight)
	for (int j = 0; j < vertices[node]->links->size(); j++)
		int_histogram(vertex_label[vertices[node]->links->l[j]], M,
		              vertices[node]->links->w[j].first,
		              vertices[node]->links->w[j].second);

	for (auto &itM : M) {
		auto itOM = label_module.find(itM.first);
		double to_fit;
		if (itM.first != vertex_label[node]) {
			to_fit = Stochastics::topological_05(itM.second.first, itOM->second.kout,
			                                     total_stubs - itOM->second.ktot,
			                                     vertices[node]->stub_number);
			to_fit *= double(dim - itOM->second.nc + 1) / (
					std::min(dim - itOM->second.nc, itOM->second.kout / itM.second.first + 1) + 1);
			if (to_fit > 1)
				to_fit = 1;
		} else {
			kop = itM.second.first;
			int kout_prime = itOM->second.kout - vertices[node]->stub_number + 2 * kop;
			to_fit = Stochastics::topological_05(itM.second.first, kout_prime,
			                                     total_stubs - itOM->second.ktot +
			                                     vertices[node]->stub_number,
			                                     vertices[node]->stub_number);
			to_fit *= double(dim - itOM->second.nc + 2) /
			          (std::min(dim - itOM->second.nc + 1, kout_prime / kop + 1) + 1);

			if (to_fit > 1)
				to_fit = 1;
			to_fit *= 0.999;        // to break possible ties

		}

		double weight_fit = Stochastics::log_together(itM.second.second, itM.second.first);
		double fitness = Stochastics::log_together(-log(to_fit) - log(weight_fit), 2);

		if (fitness < min_fitness) {

			kp = itM.second.first;
			min_fitness = fitness;
			fi = itM.first;
		}
	}
}

/**
 * Move node from current module to the new module and update the degree of the modules.
 * @param node The node to move
 * @param new_label Label of the new module
 * @param stubs_into_new Number of edges from the node to the new module
 * @param stubs_into_old Number of edges from the node to the old module
 */
inline void
oslomnet_louvain::update_modules(const int &node, const int &new_label, const int &stubs_into_new,
                                 const int &stubs_into_old) {
	if (new_label != vertex_label[node]) {
		nodes_changed++;
		for (int j = 0; j < vertices[node]->links->size(); j++)
			vertex_to_check_next[vertices[node]->links->l[j]] = true;

		auto itm = label_module.find(vertex_label[node]);
		--(itm->second.nc);
		if (itm->second.nc == 0)
			label_module.erase(itm);
		else {
			itm->second.kout -= vertices[node]->stub_number - 2 * stubs_into_old;
			itm->second.ktot -= vertices[node]->stub_number;
		}

		itm = label_module.find(new_label);
		++(itm->second.nc);
		itm->second.kout += vertices[node]->stub_number - 2 * stubs_into_new;
		itm->second.ktot += vertices[node]->stub_number;

		vertex_label[node] = new_label;
	}
}

void oslomnet_louvain::single_pass_unweighted() {
	int new_label, stubs_into_new, stubs_into_old; // new_label is the label node i likes most, stubs_into_new is the number of internal stubs in module new_label, stubs_into_old is the same for vertex_label[i]
	// new_label = best module, stubs_into_new = stubs into best module, stubs_into_old = stubs into current module
	for (int &node : vertex_order) {
		if (vertex_to_check[node]) {
			unweighted_favorite_of(node, new_label, stubs_into_new, stubs_into_old);
			update_modules(node, new_label, stubs_into_new, stubs_into_old);
		}
	}
}

void oslomnet_louvain::single_pass_weighted() {
	int fi, kp, kop;    // fi is the label node i likes most, kp is the number od internal stubs in module fi, kop is the same for vertex_label[i]
	for (int &itd : vertex_order) {
		if (vertex_to_check[itd]) {
			weighted_favorite_of(itd, fi, kp, kop);
			update_modules(itd, fi, kp, kop);
		}
	}
}

void oslomnet_louvain::set_partition_collected(std::deque<std::deque<int>> &ten2) {
	ten2.clear();
	std::deque<std::deque<int>> M;

	// take partition from vertex_label  //******************************
	std::map<int, int> mems;
	for (int i = 0; i < dim; i++) {
		std::pair<std::map<int, int>::iterator, bool> itm_bool = mems.insert(
				std::make_pair(vertex_label[i], mems.size()));
		if (itm_bool.second) {
			std::deque<int> first;
			M.push_back(first);
		}

		M[itm_bool.first->second].push_back(i);
	}

	// check if subgraphs are connected  //******************************
	for (auto &i : M) {
		std::deque<std::deque<int>> link_per_node;
		std::deque<std::deque<std::pair<int, double> > > weights_per_node;
		set_subgraph(i, link_per_node, weights_per_node);
		StaticNetwork giovanni;
		giovanni.set_graph(link_per_node, weights_per_node, i);
		std::deque<std::deque<int>> gM;
		giovanni.set_connected_components(gM);

		if (gM.size() == 1)
			ten2.push_back(i);
		else {
			for (auto &j : gM) {
				giovanni.deque_id(j);
				ten2.push_back(j);
			}
		}
	}
}

/**
 * Try to merge modules, each module is a node in the current graph.
 * @param P
 * @return
 */
int oslomnet_louvain::collect_raw_groups_once(std::deque<std::deque<int>> &P) {
	module_initializing();
	int stopper = 0;
	int previous_nodes_changed = dim;
	int iteration = 0;
	while (true) {
		nodes_changed = 0;
		for (bool &check : vertex_to_check_next) // Vertexes to check in the next iteration
			check = false;
		shuffle_s(vertex_order);

		if (paras->weighted)
			single_pass_weighted();
		else
			single_pass_unweighted();

		if (paras->print_flag_subgraph && iteration % 10 == 0)
			std::cout << "iteration: " << iteration << " number of modules: "
			          << label_module.size() << std::endl;
		++iteration;
		/* this conditions means that at least max_iteration_convergence iterations are done and, after that, the number of nodes changed has to decrease (up to a facto 1e-3) */
		if (iteration > paras->max_iteration_convergence &&
		    double(nodes_changed - previous_nodes_changed) > 1e-3 * previous_nodes_changed)
			stopper++;

		if (stopper == paras->max_iteration_convergence || nodes_changed == 0 || iteration == dim)
			break;

		vertex_to_check = vertex_to_check_next;
		previous_nodes_changed = nodes_changed;
	}

	set_partition_collected(P);

	if (paras->print_flag_subgraph)
		std::cout << "collection done " << std::endl << std::endl;

	label_module.clear();
	vertex_label.clear();
	vertex_order.clear();
	vertex_to_check.clear();
	vertex_to_check_next.clear();
	nodes_changed = 0;

	return 0;
}

void prints(map_int_om &M) {
	for (auto &itm : M)
		std::cout << "module: " << itm.first << "\t\t\t\tnc= " << itm.second.nc
		          << "\t ktot= " << itm.second.ktot << "\t kout= " << itm.second.kout
		          << std::endl;
}
